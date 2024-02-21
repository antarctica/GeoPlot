import folium
import os
import json
import copy
import logging
import pandas as pd
import numpy as np
import geopandas as gpd
import re

from folium import plugins
from branca.colormap import linear
from branca.element import MacroElement
from shapely import wkt
from shapely.geometry import Polygon
from jinja2 import Template
from pyproj import Geod
from bas_geoplot.utils import convert_decimal_days


def params_object(layer, predefined=None, **kwargs):
    # Loading config of standards
    p_file = os.path.join(os.path.dirname(__file__),'config','interactive.json')
    with open(p_file, 'r') as f:
        p = json.load(f)[layer]

    # If a standard type for the specific layer is used replace the global standard values
    if not (predefined is None):
        for key, value in p[predefined].items():
            p[key] = value

    # Replace the global standard values with user defined key arguments
    for key, value in kwargs.items():
        p[key] = value

    return p


def build_sectors(cell_poly, dcx, dcy, cx, cy):
    """
        Split a rectangular polygon into 8 triangular sectors corresponding to the 8 angular ranges used for array data
    """

    coords = list(cell_poly.exterior.coords)
    x_gap = dcx - (np.sqrt(2) - 1) * dcy
    y_gap = dcy - (np.sqrt(2) - 1) * dcx
    centre = (cx, cy)

    # Find coordinates of sector vertices on the cell edge
    edge_coords = [0, 0, 0, 0, 0, 0, 0, 0]

    edge_coords[0] = (coords[0][0], coords[0][1] + y_gap)
    edge_coords[1] = (coords[0][0], coords[0][1] + 2*dcy - y_gap)
    edge_coords[2] = (coords[1][0] + x_gap, coords[1][1])
    edge_coords[3] = (coords[1][0] + 2*dcx - x_gap, coords[1][1])
    edge_coords[4] = (coords[2][0], coords[2][1] - y_gap)
    edge_coords[5] = (coords[2][0], coords[2][1] - 2*dcy + y_gap)
    edge_coords[6] = (coords[3][0] - x_gap, coords[3][1])
    edge_coords[7] = (coords[3][0] - 2*dcx + x_gap, coords[3][1])

    # Construct polygons in the order of the cases 1,2...-4
    sec_polys = [0, 0, 0, 0, 0, 0, 0, 0]

    sec_polys[0] = Polygon([centre, edge_coords[3], coords[2], edge_coords[4]])
    sec_polys[1] = Polygon([centre, edge_coords[4], edge_coords[5], centre])
    sec_polys[2] = Polygon([centre, edge_coords[5], coords[3], edge_coords[6]])
    sec_polys[3] = Polygon([centre, edge_coords[6], edge_coords[7], centre])
    sec_polys[4] = Polygon([centre, edge_coords[7], coords[0], edge_coords[0]])
    sec_polys[5] = Polygon([centre, edge_coords[0], edge_coords[1]])
    sec_polys[6] = Polygon([centre, edge_coords[1], coords[1], edge_coords[2]])
    sec_polys[7] = Polygon([centre, edge_coords[2], edge_coords[3]])

    return sec_polys


def sectorise_df(df,dn):
    """
        Split all cells in mesh into 8 sectors corresponding to the 8 angular ranges used for array data and assign the
        correct values to each sector
    """

    sec_df = pd.DataFrame(columns=[dn,'geometry'])

    for i, row in df.iterrows():
        sec_polys = build_sectors(row['geometry'], row['dcx'], row['dcy'], row['cx'], row['cy'])
        for j, poly in enumerate(sec_polys):
            idx = i*8 + j
            sec_df.loc[idx] = {dn:row[dn][j], 'geometry':poly}
    return sec_df


def scale_list_columns(column, scaling_factor):
    """
        Takes in a column that contains at least some list values. Scales and then rounds any floats within the lists
        to 3 decimal places.
    """
    scaled_column = []
    for c in column:
        if type(c) == float:
            scaled_column.append(c)
        else:
            scaled_column.append([round(x*scaling_factor, 3) for x in c])

    return scaled_column


def round_list_columns(column):
    """
        Takes in a column that contains at least some list values and rounds any floats within the lists
         to 3 decimal places.
    """
    rounded_column = []
    for c in column:
        if type(c) == float:
            rounded_column.append(c)
        else:
            rounded_column.append([round(x, 3) for x in c])

    return rounded_column


class BindColormap(MacroElement):
    """Binds a colormap to a given layer.

    Parameters
    ----------
    colormap : branca.colormap.ColorMap
        The colormap to bind.
    """
    def __init__(self, layer, colormap):
        super(BindColormap, self).__init__()
        self.layer = layer
        self.colormap = colormap
        self._template = Template(u"""
        {% macro script(this, kwargs) %}
            {{this.colormap.get_name()}}.svg[0][0].style.display = 'block';
            {{this._parent.get_name()}}.on('layeradd', function (eventLayer) {
                if (eventLayer.layer == {{this.layer.get_name()}}) {
                    {{this.colormap.get_name()}}.svg[0][0].style.display = 'block';
                }});
            {{this._parent.get_name()}}.on('layerremove', function (eventLayer) {
                if (eventLayer.layer == {{this.layer.get_name()}}) {
                    {{this.colormap.get_name()}}.svg[0][0].style.display = 'none';
                }});
        {% endmacro %}
        """)  # noqa


class Map:
    def __init__(self,predefined=None,**kwargs):

        """
            ---

            Interactive Plotting Toolkit leveraging Folium software package. During the class call a basemap is generated and subsequent function calls from the class functions as separate layers. Once all layers are completed the plot can either be shown in an interactive notebook or saved for viewing later

            This plotting class includes:
                Maps: Plotting Geopandas dataframes of polygon datatypes e.g. Environmental Mesh
                Paths: Plotting GeoJSONs of line information e.g. Route Paths
                Points: Plotting padnas dataframe of small point information e.g. waypoints
                Vectors: Plotting pandas dataframe of small vector information e.g. Currents
                MeshInfo: Plotting Geopandas dataframe labelling the data with all its attributes.

            ---

            Attributes:
                predefined (opt=None, string) - Predefined plotting formats given in
                    config/interactive.json of the package files
                **kwargs - Can be used to change information within the configuration files. See manual for more information.
            ---

        """

        # === Initialising layer info
        self._layer_info = {}

        # ==== Loading standard configs
        p = params_object('Basemap', predefined=predefined, **kwargs)

        if p['title']:
            title='{} &ensp;|&ensp;'.format(p['title'])
        else:
            title=''

        title_html = '''
            <h1 style="color:#003b5c;font-size:16px">
            &ensp;<img src="https://i.ibb.co/XtZdzDt/BAS-colour-eps.png" alt="BAS-colour-eps" border="0" style="width:179px;height:40px;"> 
            &ensp; |&ensp; {}
            </h1>
            </body>
            '''.format(title)

        if 'map_centre' in p.keys():
            map_centre = p['map_centre']
        else:
            map_centre = None

        if 'size' in p.keys():
            self.map = folium.Map(location=map_centre,zoom_start=p['zoom_start'],tiles=None,width=p['size'][0],height=p['size'][1])
        else:
            self.map = folium.Map(location=map_centre,zoom_start=p['zoom_start'],tiles=None)

        if p['offline_filepath']:
            self._offline_mode = True
            self._offline_mode_path = p['offline_filepath']
        else:
            self._offline_mode = False

        if p['basemap']:
            bsmap = folium.FeatureGroup(name='BaseMap')
            folium.TileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png',attr="toner-bcg", name='Basemap').add_to(bsmap)
            bsmap.add_to(self.map)

        if p['offline_coastlines']:
            bsmap = folium.FeatureGroup(name='Coastlines',show=True)
            antarctica = gpd.read_file(p['offline_coastlines'])
            folium.GeoJson(antarctica,
                    style_function=lambda feature: {
                        'color': 'black',
                        'weight': 0.5,
                    }, name="geojson").add_to(bsmap)
            bsmap.add_to(self.map)

        if (p['plot_title']) and (p['title'] is not None):
            self.map.get_root().html.add_child(folium.Element(title_html))

    def _layer(self,name,show=False):
        if name not in self._layer_info:
            lyr = folium.FeatureGroup(name='{}'.format(name),show=show)
            self._layer_info[name] = lyr
        return self._layer_info[name]

    def _add_plots_map(self):
        for key in self._layer_info.keys():
            self._layer_info[key].add_to(self.map)

    def show(self):
        """
            For use case in interactive notebooks, showing the plot.
        """

        self._add_plots_map()
        folium.LayerControl('topleft', collapsed=True).add_to(self.map)
        return self.map

    def save(self,file):
        """
            Saving the interactive plot to file

            Attributes:
                file (str): File path for output
        """
        map = self.show()

        html = map.get_root().render()
        if self._offline_mode:
            html = re.sub(r'https://cdn.jsdelivr.net/npm/leaflet@\d+\.\d+\.\d+/dist/leaflet.js',
                          f'{os.path.join(self._offline_mode_path, "leaflet.js")}', html)
            html = re.sub(r'https://cdn.jsdelivr.net/npm/bootstrap@\d+\.\d+\.\d+/dist/js/bootstrap.bundle.min.js',
                          f'{os.path.join(self._offline_mode_path, "bootstrap.bundle.min.js")}', html)
            html = re.sub(r'https://cdn.jsdelivr.net/npm/bootstrap@\d+\.\d+\.\d+/dist/css/bootstrap.min.css',
                          f'{os.path.join(self._offline_mode_path, "bootstrap.min.css")}', html)
            html = re.sub(r'https://netdna.bootstrapcdn.com/bootstrap/\d+\.\d+\.\d+/css/bootstrap.min.css',
                          f'{os.path.join(self._offline_mode_path, "bootstrap_netdna.min.css")}', html)
            html = re.sub(r'https://cdn.jsdelivr.net/npm/@fortawesome/fontawesome-free@\d+\.\d+\.\d+/css/all.min.css',
                          f'{os.path.join(self._offline_mode_path, "all.min.css")}', html)
            html = html.replace('https://code.jquery.com/jquery-1.12.4.min.js',
                                f'{os.path.join(self._offline_mode_path, "jquery-1.12.4.min.js")}')
            html = re.sub(r'https://maxcdn.bootstrapcdn.com/bootstrap/\d+\.\d+\.\d+/js/bootstrap.min.js',
                          f'{os.path.join(self._offline_mode_path, "bootstrap.min.js")}', html)
            html = re.sub(r'https://cdnjs.cloudflare.com/ajax/libs/Leaflet.awesome-markers/\d+\.\d+\.\d+/leaflet.awesome-markers.js',
                          f'{os.path.join(self._offline_mode_path, "leaflet.awesome-markers.js")}', html)
            html = re.sub(r'https://cdn.jsdelivr.net/npm/leaflet@\d+\.\d+\.\d+/dist/leaflet.css',
                          f'{os.path.join(self._offline_mode_path, "leaflet.css")}', html)
            html = re.sub(r'https://maxcdn.bootstrapcdn.com/bootstrap/\d+\.\d+\.\d+/css/bootstrap.min.css',
                          f'{os.path.join(self._offline_mode_path, "bootstrap.min.css")}', html)
            html = re.sub(r'https://maxcdn.bootstrapcdn.com/bootstrap/\d+\.\d+\.\d+/css/bootstrap-theme.min.css',
                          f'{os.path.join(self._offline_mode_path, "bootstrap-theme.min.css")}', html)
            html = re.sub(r'https://maxcdn.bootstrapcdn.com/font-awesome/\d+\.\d+\.\d+/css/font-awesome.min.css',
                          f'{os.path.join(self._offline_mode_path, "font-awesome.min.css")}', html)
            html = re.sub(r'https://cdnjs.cloudflare.com/ajax/libs/Leaflet.awesome-markers/\d+\.\d+\.\d+/leaflet.awesome-markers.css',
                          f'{os.path.join(self._offline_mode_path, "leaflet.awesome-markers.css")}', html)
            html = html.replace("https://cdn.jsdelivr.net/gh/python-visualization/folium/folium/templates/leaflet.awesome.rotate.min.css",
                                f'{os.path.join(self._offline_mode_path, "leaflet.awesome.rotate.min.css")}')

            html = re.sub(r'https://cdnjs.cloudflare.com/ajax/libs/d3/\d+\.\d+\.\d+/d3.min.js',
                          f'{os.path.join(self._offline_mode_path, "d3.min.js")}', html)
            html = re.sub(r'https://cdnjs.cloudflare.com/ajax/libs/leaflet-dvf/\d+\.\d+\.\d+/leaflet-dvf.markers.min.js',
                          f'{os.path.join(self._offline_mode_path, "leaflet-dvf.markers.min.js")}', html)
            html = html.replace('https://i.ibb.co/XtZdzDt/BAS-colour-eps.png',
                                f'{os.path.join(self._offline_mode_path, "BAS-colour-eps.png")}')

            html = html.replace('https://cdn.jsdelivr.net/gh/marslan390/BeautifyMarker/leaflet-beautify-marker-icon.min.js',
                                f'{os.path.join(self._offline_mode_path, "leaflet-beautify-marker-icon.min.js")}')
            html = html.replace("https://cdn.jsdelivr.net/gh/marslan390/BeautifyMarker/leaflet-beautify-marker-icon.min.css",
                                f'{os.path.join(self._offline_mode_path, "leaflet-beautify-marker-icon.min.css")}')

        with open(file,'w') as fp:
            fp.write(html)
            fp.close()

    def fit_to_bounds(self, bounds=None):
        """
            Change zoom level to match plotted data
        """
        # Use input bounds if provided
        if bounds:
            self.map.fit_bounds(bounds)
        # Otherwise try to retrieve bounds automatically from map
        else:
            bounds = self.map.get_bounds()
            if bounds == [[None, None], [None, None]]:
                logging.info("No bounds returned, can't fit to layers")
            else:
                self.map.fit_bounds(bounds)

    def Paths(self,geojson,name,show=True,predefined=None,arrows=False,**kwargs):
        """
            Overlays paths on the interactive plot with layer defined by `name`

            Attributes:
                geojson (dict): A geojson file with several features representing all
                    the separate paths
                name (string): Layer name to add to the interactive plot
                show (opt=True, boolean) - Show the layer on loading of plot
                predefined (opt=None, string) - Predefined plotting formats given in
                    config/interactive.json of the package files
        """

        p = params_object('Paths', predefined=predefined, **kwargs)

        # Defining the feature groups to add
        pths = self._layer(name,show=show)
        paths = geojson['features']

        no_path_name = True
        for path in copy.deepcopy(paths):
            if p['data_name'] in path['properties'].keys():
                no_path_name = False
        if no_path_name:
            return

        # Determining min-max values of all paths if colormap being used
        if type(p['line_color']) is dict:
            max_val = -np.inf
            min_val = np.inf
            for path in copy.deepcopy(paths):
                    if (np.array(path['properties'][p['data_name']])*p['scaling_factor']).max() > max_val:
                        max_val = (np.array(path['properties'][p['data_name']])*p['scaling_factor']).max()

                    if (np.array(path['properties'][p['data_name']])*p['scaling_factor']).min() < min_val:
                        min_val = (np.array(path['properties'][p['data_name']])*p['scaling_factor']).min()
        # Determining max travel-times of all paths
        for path in paths:
            points   = np.array(path['geometry']['coordinates'])

            start_wpt = path['properties']['from']
            end_wpt   = path['properties']['to']

            if p['data_name']:
                try:
                    data_val = np.array(path['properties'][p['data_name']])*p['scaling_factor']
                except:
                    continue
            else:
                data_val = np.array(len(points))

            # Find the max value for this path for display in pop-up
            path_max = np.max(data_val)

            points[:,0] = points[:,0]
            points = points[:,::-1]

            if type(p['line_color']) is dict:
                if "cmin" in p["line_color"].keys():
                    min_val = p["line_color"]['cmin']
                if "cmax" in p["line_color"].keys():
                    max_val = p["line_color"]['cmax']

                colormap = linear._colormaps[p["line_color"]['color']].scale(min_val,max_val)
                if p['unit'] == 'Days':
                    colormap.caption = '{} ({}, Max Value: {})'.format(name, p['unit'], convert_decimal_days(max_val))
                else:
                    colormap.caption = '{} ({}, Max Value: {:.3f})'.format(name,p['unit'],max_val)
                folium.ColorLine(points,data_val, colormap=colormap,nb_steps=50, weight=p['line_width'],
                                 opacity=p['line_opacity']).add_to(pths)

                if p['unit'] == 'Days':
                    folium.PolyLine(points, color='black', weight=p['line_width'], opacity=0.0,
                                    popup = "Path - {} to {}\n{} = {}".format(start_wpt, end_wpt, p['data_name'],
                                                                              convert_decimal_days(path_max))
                                    ).add_to(pths)
                else:
                    folium.PolyLine(points, color='black', weight=p['line_width'], opacity=0.0,
                                    popup="Path - {} to {}\n{} = {:.3f} {}".format(start_wpt, end_wpt, p['data_name'],
                                                                                   path_max, p['unit'])).add_to(pths)

                if arrows:
                    # Every 10 segments, draw a triangle as an arrow head
                    for idx in range(1, len(points), 10):
                        lon_diff = points[idx, 0] - points[idx-1, 0]
                        lat_diff = points[idx, 1] - points[idx-1, 1]
                        loc = [points[idx,0],points[idx,1]]
                        heading = -np.degrees(np.arctan2(lon_diff, lat_diff))
                        folium.RegularPolygonMarker(location=loc, rotation=heading,
                                                    color=colormap(data_val[idx]), fill=True,
                                                    number_of_sides=3, radius=10).add_to(pths)

            else:
                folium.PolyLine(points, color=p['line_color'], weight=p['line_width'], opacity=p['line_opacity'],
                                popup = "Path - {} to {}".format(start_wpt,end_wpt)).add_to(pths)

                if arrows:
                    # Every 10 segments, draw a triangle as an arrow head
                    for idx in range(1, len(points), 10):
                        lon_diff = points[idx, 0] - points[idx-1, 0]
                        lat_diff = points[idx, 1] - points[idx-1, 1]
                        loc = [points[idx,0],points[idx,1]]
                        heading = -np.degrees(np.arctan2(lon_diff, lat_diff))
                        folium.RegularPolygonMarker(location=loc, rotation=heading,
                                                    color=p['line_color'], fill=True,
                                                    number_of_sides=3, radius=10).add_to(pths)

            if p['path_points']:
                for idx in range(len(points)):
                    loc = [points[idx,0],points[idx,1]]
                    folium.Marker(
                        location=loc,
                        icon=folium.plugins.BeautifyIcon(icon='circle',
                                                    border_color='transparent',
                                                    background_color='transparent',
                                                    border_width=1,
                                                    text_color='black',
                                                    inner_icon_style='margin:0px;font-size:0.8em')
                    ).add_to(pths)
        
        if type(p['line_color']) is dict:
            self.map.add_child(colormap)
            self.map.add_child(BindColormap(pths,colormap))

    def Points(self,dataframe_points,name,show=True,predefined=None,**kwargs):
        """
            Overlays small collection of Points, such as waypoints and sites of interest.

            Attributes:
                dataframe_points (Pandas DataFrame): A Dataframe requiring at least columns of Latitude ('Lat') Longitude ('Long'), Name ('Name').
                name (string): Layer name to add to the interactive plot
                show (opt=True, boolean) - Show the layer on loading of plot
                predefined (opt=None, string) - Predefined plotting formats given in
                    config/interactive.json of the package files
        """

        p = params_object('Points', predefined=predefined, **kwargs)

        wpts      = self._layer(name,show=show)
 
        for id,wpt in dataframe_points.iterrows():
            loc = [wpt['Lat'], wpt['Long']]
            folium.Marker(
                location=loc,
                icon=plugins.BeautifyIcon(icon=p['icon'],
                                                border_color='transparent',
                                                background_color='transparent',
                                                border_width=p['line_width'],
                                                text_color=p['color'],
                                                inner_icon_style='margin:0px;font-size:{}em'.format(p["marker_size"])),
                popup="Name = {}\n Long = {:4f}\n Lat = {:4f}\n".format(wpt['Name'],loc[0],loc[1]),
            ).add_to(wpts)    

            if p['names']:
                folium.Marker(
                            location=loc,
                                icon=folium.features.DivIcon(
                                    icon_size=(250,36),
                                    icon_anchor=(0,0),
                                    html='<div style="font-size: {}pt">{}</div>'.format(p['names']['font_size'],wpt['Name']),
                                    ),
                ).add_to(wpts)

        wpts.add_to(self.map)

    def Maps(self,dataframe_pandas,name,show=True,predefined=None,plot_sectors=False,**kwargs):
        """
            Overlays a layer of scalars or booleans such as sea ice concentration or land.

            Attributes:
                dataframe_pandas (GeoPandas DataFrame): A Dataframe in Geopandas format. This requires at least a column
                    with 'geometry'. Additional plotting of specific columns is done in the configuration files.
                name (string): Layer name to add to the interactive plot
                show (opt=True, boolean): Show the layer on loading of plot
                predefined (opt=None, string): Predefined plotting formats given in config/interactive.json of the
                    package files
                plot_sectors (boolean): Display sectorised list values as eight individual polygons
        """
        p = params_object('Maps', predefined=predefined, **kwargs)

        dataframe_pandas = copy.copy(dataframe_pandas)
        dataframe_pandas['geometry'] = dataframe_pandas['geometry'].apply(wkt.loads)

        # Don't plot anything in the land cells unless, of course, we are plotting the land mask
        if 'land' in dataframe_pandas.keys() and p['data_name'] != 'land':
            dataframe_pandas = dataframe_pandas[dataframe_pandas['land']==False].reset_index(drop=True)

        # For array values we either plot each value as a separate polygon or just the average of the values in the list
        if any(type(d) is list for d in dataframe_pandas[p['data_name']]):
            if plot_sectors:
                dataframe_pandas = sectorise_df(dataframe_pandas, p['data_name'])
            else:
                dataframe_pandas[p['data_name']] = [np.mean(dn) for dn in dataframe_pandas[p['data_name']]]
        dataframe_geo = gpd.GeoDataFrame(dataframe_pandas,crs='EPSG:4326', geometry='geometry')

        feature_info = self._layer(name,show=show)

        if p['data_name']:
            try:
                if 'scaling_factor' in p.keys():
                    dataframe_geo[p['data_name']] = dataframe_geo[p['data_name']]*p['scaling_factor']
                else:
                    dataframe_geo[p['data_name']] = dataframe_geo[p['data_name']]
            except:
                raise print('Data name not in variables')

        if p['data_name'] and p['trim_min']:
            dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']] > p['trim_min']]
        if p['data_name'] and p['trim_max']:
            dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']] < p['trim_max']]

        if (type(p['fill_color']) is dict) and (p['data_name']):
            dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']].notna() & ~np.isinf(abs(dataframe_geo[p['data_name']]))]
            if 'cmin' in p['fill_color'].keys():
                cmin = p['fill_color']['cmin']
            else:
                cmin = dataframe_geo[p['data_name']].min()
            if 'cmax' in p['fill_color'].keys():
                cmax = p['fill_color']['cmax']
            else:
                cmax = dataframe_geo[p['data_name']].max()
            colormap = linear._colormaps[p['fill_color']['colormap']].scale(cmin,cmax)

            if len(dataframe_geo[p['data_name']].unique()) == 1:
                colormap.caption = '{} ({}, Singular Value = {:.3f})'.format(name,p['units'],dataframe_geo[p['data_name']].unique()[0])
            else:
                colormap.caption = '{} ({})'.format(name,p['units'])

            folium.GeoJson(
                dataframe_geo,
                style_function=lambda x: {
                        'fillColor': colormap(x['properties'][p["data_name"]]),
                        'color': p["line_color"],
                        'weight': p["line_width"],
                        'fillOpacity': p["fill_opacity"]
                    }
            ).add_to(feature_info)
            self.map.add_child(colormap)
            self.map.add_child(BindColormap(feature_info,colormap))
        else:
            if p['data_name'] and (type(dataframe_geo[p['data_name']].iloc[0].item()) is bool):
                dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']] == True]
            if not dataframe_geo.empty:
                folium.GeoJson(
                    dataframe_geo,
                    style_function=lambda x: {
                        'fillColor': p['fill_color'],
                        'color': p['line_color'],
                        'weight': p['line_width'],
                        'fillOpacity': p['fill_opacity']
                        }
                ).add_to(feature_info)

    def Vectors(self,mesh,name,show=True,predefined=None,**kwargs):
        """
            Overlays a layer of vectors such as currents.

            Attributes:
                dataframe_points (Pandas DataFrame): A Dataframe requiring at least columns of Longitude ('cx') Latitude ('cy'), Displacement in X ('uC') and Displacement in Y ('vC').
                name (string): Layer name to add to the interactive plot
                show (opt=True, boolean) - Show the layer on loading of plot
                predefined (opt=None, string) - Predefined plotting formats given in
                    config/interactive.json of the package files
        """

        p = params_object('Vectors', predefined=predefined, **kwargs)

        Vectors = mesh
        # Filter out empty vectors but allow single component vectors
        Vectors = Vectors[(Vectors[p['V']]!=0.0)|(Vectors[p['U']]!=0.0)].reset_index(drop=True)
        Vectors = Vectors.dropna(subset=[p['U'], p['V']]).reset_index(drop=True)

        if 'inaccessible' in Vectors.keys():
            Vectors = Vectors[Vectors['inaccessible']==False].reset_index(drop=True)

        vcts = self._layer(name,show=show)
        for idx,vec in Vectors.iterrows():
            loc =[[vec[p['Lat']],vec[p['Long']]],[vec[p['Lat']]+vec[p['V']]*p['scale'],vec[p['Long']]+vec[p['U']]*p['scale']]]
            mag = np.sqrt(vec[p['V']]**2 + vec[p['U']]**2)
            folium.PolyLine(loc, color=p['color'],weight=1.4, popup='{} (m/s)'.format(mag)).add_to(vcts)
            # get pieces of the line
            pairs = [(loc[idx], loc[idx-1]) for idx, val in enumerate(loc) if idx != 0]
            # get rotations from forward azimuth of the line pieces and add an offset of 90°
            geodesic = Geod(ellps='WGS84')
            rotations = [geodesic.inv(pair[0][1], pair[0][0], pair[1][1], pair[1][0])[0]+90 for pair in pairs]
            # create your arrow
            for pair, rot in zip(pairs, rotations):
                folium.RegularPolygonMarker(location=pair[0], color=p['color'], fill=True, fill_color=p['color'], fill_opacity=1,
                                            number_of_sides=3, rotation=rot,radius=4,weight=p['line_width']).add_to(vcts)

    def MeshInfo(self,mesh,name,predefined='PolarRoute',show=True,**kwargs):
        """
            Overlays a layer that gives information on all polygon boxes.

            Attributes:
                cellboxes (GeodataFrame): A Dataframe in Geopandas format. This requires at least a column with 'geometry'.
                name (string): Layer name to add to the interactive plot
                show (opt=True, boolean) - Show the layer on loading of plot
                predefined (opt=None, string) - Predefined plotting formats given in
                    config/interactive.json of the package files
        """

        p = params_object('MeshInfo', predefined=predefined, **kwargs)

        dataframe_pandas = copy.copy(pd.DataFrame(mesh))
        dataframe_pandas['geometry'] = dataframe_pandas['geometry'].apply(wkt.loads)
        dataframe_geo = gpd.GeoDataFrame(dataframe_pandas,crs='EPSG:4326', geometry='geometry')

        p = p['fields']

        column_names = ['geometry']
        for col_info in p:
            try:
                column_name = col_info['Name']
                data = col_info['data']
                dataframe_geo[column_name] = dataframe_geo[data]
                # Apply scaling factor
                if 'scaling_factor' in col_info.keys():
                    if pd.api.types.is_numeric_dtype(dataframe_geo[column_name]):
                        dataframe_geo[column_name] = dataframe_geo[column_name]*col_info['scaling_factor']
                    else:
                        dataframe_geo[column_name] = scale_list_columns(dataframe_geo[column_name], col_info['scaling_factor'])
                # Round values within array data types
                elif any(type(d) is list for d in dataframe_geo[column_name]):
                    dataframe_geo[column_name] = round_list_columns(dataframe_geo[column_name])
                column_names.append(column_name)
            # Catch KeyError when one of the default Mesh Info fields is not present
            except KeyError:
                continue
        dataframe_geo = dataframe_geo[column_names]

        feature_info = self._layer(name,show=show)
        folium.GeoJson(dataframe_geo,
            style_function= lambda x: {
                    'fillColor': 'black',
                    'color': 'black',
                    'weight': 0.0,
                    'fillOpacity': 0.0
                },
            tooltip=folium.GeoJsonTooltip(
                fields=list(dataframe_geo.columns[1:]),
                aliases=list(dataframe_geo.columns[1:]),
                localize=True
            )
        ).add_to(feature_info)


def plot_mesh(mesh_filename, basemap=False, mesh_info=False, **kwargs):
    """
        Takes in a mesh and returns a Map object plotting a selection of common data interactively.

        Args:
            mesh_filename (str or dict): Input mesh as either a file path or dictionary
            basemap (bool): Determine whether to plot the openstreetmap basemap layer
            mesh_info (bool): Determine whether to include the mesh info pop-up
        Returns:
            mp (Map): Interactive Map object with all relevant data plotted on it
    """

    if isinstance(mesh_filename, str):
        with open(mesh_filename, 'r') as fl:
            info = json.load(fl)
    elif isinstance(mesh_filename, dict):
        info = mesh_filename
    else:
        raise TypeError('Input to plot_mesh should be either a filename as a string or a dict object!')

    mesh = pd.DataFrame(info['cellboxes'])
    region = info['config']['mesh_info']['region']
    split_level = info['config']['mesh_info']['splitting']['split_depth']
    if 'rm_title' in kwargs and kwargs['rm_title'] == True:
        output=None
    else:
        output = 'Start Date: {}, End Date: {} | Split level: {}'.format(region['start_time'],
                                                                            region['end_time'], split_level)

    # Put mesh bounds in format required by fit_to_bounds
    mesh_bounds = [[region["lat_min"], region["long_min"]], [region["lat_max"], region["long_max"]]]

    mp = Map(title=output,basemap=basemap, **kwargs)

    # Plot maps of mesh info
    if 'cx' in mesh.columns:
        logging.debug('Plotting mesh grid')
        mp.Maps(mesh, 'MeshGrid', predefined='cx')
    if 'SIC' in mesh.columns:
        logging.debug("Plotting Sea Ice Concentration")
        mp.Maps(mesh, 'SIC', predefined='SIC')
    if 'ext_ice' in mesh.columns:
        logging.debug("Plotting Extreme Ice areas")
        mp.Maps(mesh, 'Extreme Ice', predefined='Extreme Sea Ice Conc')
    if 'land' in mesh.columns:
        logging.debug("Plotting Land Mask")
        mp.Maps(mesh, 'Land Mask', predefined='Land Mask')
    if 'shallow' in mesh.columns:
        logging.debug("Plotting shallow areas")
        mp.Maps(mesh, 'Shallows', predefined='Shallows')
    if 'elevation' in mesh.columns:
        logging.debug("Plotting elevation")
        mp.Maps(mesh, 'Elevation', predefined='Elev', show=False)
    if 'fuel' in mesh.columns:
        logging.debug('Plotting Fuel usage per day and tCO2e')
        mp.Maps(mesh, 'Fuel', predefined='Fuel (Tonnes/Day)', show=False)
        mp.Maps(mesh, 'tCO2e', predefined='tCO2e', show=False)
    if 'battery' in mesh.columns:
        logging.debug('Plotting battery usage')
        mp.Maps(mesh, 'Battery Usage', predefined='Battery Usage', show=False)
    if 'speed' in mesh.columns:
        logging.debug('Plotting vessel maximum speed')
        mp.Maps(mesh, 'Max Speed', predefined='Max Speed (knots)', show=False)
    if ('uC' in mesh.columns) and ('vC' in mesh.columns):
        mesh['mC'] = np.sqrt(mesh['uC']**2 + mesh['vC']**2)
        logging.debug('Plotting currents')
        mp.Vectors(mesh,'Currents', show=False, predefined='Currents')
    if ('u10' in mesh.columns) and ('v10' in mesh.columns):
        mesh['m10'] = np.sqrt(mesh['u10'] ** 2 + mesh['v10'] ** 2)
        mp.Vectors(mesh, 'Winds', predefined='Winds', show=False)
        mp.Maps(mesh, 'Wind Speed', predefined='Wind Speed', show=False)
        logging.debug('Plotting winds')
    if 'swh' in mesh.columns:
        logging.debug("Plotting Significant Wave Height")
        mp.Maps(mesh, 'Significant Wave Height', predefined='Significant Wave Height')
    if 'hmax' in mesh.columns:
        logging.debug("Plotting Max Wave Height")
        mp.Maps(mesh, 'Max Wave Height', predefined='Max Wave Height', show=False)
    if 'mwd' in mesh.columns:
        logging.debug("Plotting Mean Wave Direction")
        mp.Maps(mesh, 'Mean Wave Direction', predefined='Wave Direction', show=False)
    if 'mwp' in mesh.columns:
        logging.debug("Plotting Mean Wave Period")
        mp.Maps(mesh, 'Mean Wave Period', predefined='Wave Period', show=False)
    if 'wind_dir' in mesh.columns:
        logging.debug("Plotting Wind Direction")
        mp.Maps(mesh, 'Wind Direction', predefined='Wind Direction', show=False)
    if 'wind_mag' in mesh.columns:
        logging.debug("Plotting Wind Magnitude")
        mp.Maps(mesh, 'Wind Magnitude', predefined='Wind Magnitude', show=False)
    if ('uW' in mesh.columns) and ('vW' in mesh.columns):
        mp.Vectors(mesh, 'Wave Direction', predefined='Wave Direction', show=False)
        logging.debug('Plotting wave direction')
    if 'ext_waves' in mesh.columns:
        logging.debug("Plotting Extreme Wave areas")
        mp.Maps(mesh, 'Extreme Waves', predefined='Extreme Waves')
    if 'offshore_platforms' in mesh.columns:
        logging.debug("Plotting Offshore Platforms")
        mp.Maps(mesh, 'Offshore Platforms', predefined='Offshore Platforms')

    # Plot routes and waypoints
    if 'paths' in info.keys():
        logging.debug('Plotting paths')
        paths = info['paths']
        mp.Paths(paths, 'Route - Traveltime', predefined='Traveltime (Days)')
        mp.Paths(paths, 'Route - Distance', predefined='Distance (Nautical miles)', show=False)
        mp.Paths(paths, 'Route - Max Speed', predefined='Max Speed (knots)', show=False)
        if 'fuel' in mesh.columns:
            mp.Paths(paths, 'Route - Fuel', predefined='Fuel', show=False)
            mp.Paths(paths, 'Route - tCO2e', predefined='tCO2e', show=False)
        if 'battery' in mesh.columns:
            mp.Paths(paths, 'Route - Battery', predefined='Battery', show=False)
    if 'waypoints' in info.keys():
        logging.debug('Plotting waypoints')
        waypoints = pd.DataFrame(info['waypoints'])
        mp.Points(waypoints, 'Waypoints', names={"font_size":10.0})

    if mesh_info:
        mp.MeshInfo(mesh, 'Mesh Info', show=False)

    mp.fit_to_bounds(mesh_bounds)

    return mp
