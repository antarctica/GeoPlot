import folium
import os
import json
import copy
import pandas as pd

import numpy as np
from folium import plugins
from folium.plugins import TimestampedGeoJson

from branca.colormap import linear
from branca.element import MacroElement
from shapely import wkt
from shapely.geometry import Polygon
import geopandas as gpd
from jinja2 import Template
from pyproj import Geod




def paramsObject(layer,predefined=None,**kwargs):
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

    sec_df = pd.DataFrame(columns=[dn,'geometry'])

    for i, row in df.iterrows():
        sec_polys = build_sectors(row['geometry'], row['dcx'], row['dcy'], row['cx'], row['cy'])
        for j, poly in enumerate(sec_polys):
            idx = i*8 + j
            sec_df.loc[idx] = {dn:row[dn][j], 'geometry':poly}
    return sec_df


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
        p = paramsObject('Basemap',predefined=predefined,**kwargs)

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

        if p['plot_title']:
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
        folium.LayerControl(collapsed=True).add_to(self.map)
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
            html = html.replace('https://cdn.jsdelivr.net/npm/leaflet@1.6.0/dist/leaflet.js',f'{os.path.join(self._offline_mode_path, "leaflet.js")}')
            html = html.replace('https://code.jquery.com/jquery-1.12.4.min.js',f'{os.path.join(self._offline_mode_path, "jquery-1.12.4.min.js")}')
            html = html.replace('https://maxcdn.bootstrapcdn.com/bootstrap/3.2.0/js/bootstrap.min.js',f'{os.path.join(self._offline_mode_path, "bootstrap.min.js")}')
            html = html.replace('https://cdnjs.cloudflare.com/ajax/libs/Leaflet.awesome-markers/2.0.2/leaflet.awesome-markers.js',f'{os.path.join(self._offline_mode_path, "leaflet.awesome-markers.js")}')

            html = html.replace('https://cdn.jsdelivr.net/npm/leaflet@1.6.0/dist/leaflet.css',f'{os.path.join(self._offline_mode_path, "leaflet.css")}')
            html = html.replace('https://maxcdn.bootstrapcdn.com/bootstrap/3.2.0/css/bootstrap.min.css',f'{os.path.join(self._offline_mode_path, "bootstrap.min.css")}')
            html = html.replace("https://maxcdn.bootstrapcdn.com/bootstrap/3.2.0/css/bootstrap-theme.min.css",f'{os.path.join(self._offline_mode_path, "bootstrap-theme.min.css")}')
            html = html.replace("https://maxcdn.bootstrapcdn.com/font-awesome/4.6.3/css/font-awesome.min.css",f'{os.path.join(self._offline_mode_path, "font-awesome.min.css")}')
            html = html.replace("https://cdnjs.cloudflare.com/ajax/libs/Leaflet.awesome-markers/2.0.2/leaflet.awesome-markers.css",f'{os.path.join(self._offline_mode_path, "leaflet.awesome-markers.css")}')
            html = html.replace("https://cdn.jsdelivr.net/gh/python-visualization/folium/folium/templates/leaflet.awesome.rotate.min.css",f'{os.path.join(self._offline_mode_path, "leaflet.awesome.rotate.min.css")}')

            html = html.replace('https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.5/d3.min.js',f'{os.path.join(self._offline_mode_path, "d3.min.js")}')
            html = html.replace('https://cdnjs.cloudflare.com/ajax/libs/leaflet-dvf/0.3.0/leaflet-dvf.markers.min.js',f'{os.path.join(self._offline_mode_path, "leaflet-dvf.markers.min.js")}')
            html = html.replace('https://i.ibb.co/XtZdzDt/BAS-colour-eps.png',f'{os.path.join(self._offline_mode_path, "BAS-colour-eps.png")}')

            html = html.replace('https://cdn.jsdelivr.net/gh/marslan390/BeautifyMarker/leaflet-beautify-marker-icon.min.js',f'{os.path.join(self._offline_mode_path, "leaflet-beautify-marker-icon.min.js")}')
            html = html.replace("https://cdn.jsdelivr.net/gh/marslan390/BeautifyMarker/leaflet-beautify-marker-icon.min.css",f'{os.path.join(self._offline_mode_path, "leaflet-beautify-marker-icon.min.css")}')

            
            
        with open(file,'w') as fp:
            fp.write(html)
            fp.close()


    def Paths(self,geojson,name,show=True,predefined=None,**kwargs):
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
        p = paramsObject('Paths',predefined=predefined,**kwargs)

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

            points[:,0] = points[:,0]
            points = points[:,::-1]

            if type(p['line_color']) is dict:
                if "cmin" in p["line_color"].keys():
                    min_val = p["line_color"]['cmin']
                if "cmax" in p["line_color"].keys():
                    max_val = p["line_color"]['cmax']

                colormap = linear._colormaps[p["line_color"]['color']].scale(min_val,max_val)
                colormap.caption = '{} ({},Max Value={:.3f})'.format(name,p['unit'],max_val)
                folium.ColorLine(points,data_val,colormap=colormap,nb_steps=50, weight=p['line_width'], opacity=p['line_opacity']).add_to(pths)

                folium.PolyLine(points,color='black', weight=p['line_width'],opacity=0.0,popup = "Path - {} to {}\n{} = {:.3f} {}".format(start_wpt,end_wpt,p['data_name'],max_val,p['unit'])).add_to(pths)

            else:
                folium.PolyLine(points,color=p['line_color'], weight=p['line_width'], opacity=p['line_opacity'],popup = "Path - {} to {}".format(start_wpt,end_wpt)).add_to(pths)


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

        p = paramsObject('Points',predefined=predefined,**kwargs)

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
        p = paramsObject('Maps',predefined=predefined,**kwargs)

        dataframe_pandas = copy.copy(dataframe_pandas)
        dataframe_pandas['geometry'] = dataframe_pandas['geometry'].apply(wkt.loads)

        # Don't plot anything in the land cells unless, of course, we are plotting the land mask
        if 'land' in dataframe_pandas.keys() and p['data_name'] != 'land':
            dataframe_pandas = dataframe_pandas[dataframe_pandas['land']==False].reset_index(drop=True)

        # For array values we either plot each value as a separate polygon or just the average of the values in the list
        if type(dataframe_pandas[p['data_name']][0]) is list:
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


        p = paramsObject('Vectors',predefined=predefined,**kwargs)

        Vectors = mesh
        Vectors = Vectors[(Vectors[p['V']]!=0.0)&(Vectors[p['U']]!=0.0)].reset_index(drop=True)
        Vectors = Vectors.dropna(subset=[p['U'], p['V']]).reset_index(drop=True)

        if 'land' in Vectors.keys():
            Vectors = Vectors[Vectors['land']==False].reset_index(drop=True)

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

        p = paramsObject('MeshInfo',predefined=predefined,**kwargs)



        dataframe_pandas = pd.DataFrame(mesh)
        dataframe_pandas['geometry'] = dataframe_pandas['geometry'].apply(wkt.loads)
        dataframe_geo = gpd.GeoDataFrame(dataframe_pandas,crs='EPSG:4326', geometry='geometry')

        p = p['fields']

        column_names = []
        column_names.append('geometry')
        for col_info in p:
            try:
                column_name = col_info['Name']
                data = col_info['data']
                dataframe_geo[column_name] = dataframe_geo[data]
                if 'scaling_factor' in col_info.keys():
                    dataframe_geo[column_name] = dataframe_geo[column_name]*col_info['scaling_factor']
                column_names.append(column_name)
            except:
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

    def _MapArray(self,array,bounds,name,show=True):
        feature_info = self._layer('{}'.format(name),show=show)
        colormap = linear._colormaps['BuPu_09'].scale(0,100)
        colormap.caption = '{} (%)'.format(name)
        Zc = np.zeros((array.shape[0],array.shape[1],4))
        for ii in range(Zc.shape[0]):
            for jj in range(Zc.shape[1]):
                z = array[ii,jj]
                if np.isnan(z):
                    Zc[ii,jj,:] = 0.0
                else:
                    Zc[ii,jj,:] = list(colormap.rgba_floats_tuple(z))


        folium.raster_layers.ImageOverlay(
            image=Zc[::-1],
            bounds=bounds,
            mercator_project=True,
            opacity=0.6
        ).add_to(feature_info)
        self.map.add_child(colormap)
        self.map.add_child(BindColormap(feature_info,colormap))




    def _Geotiff(self,path,name,show=True):
        import rasterio
        src = rasterio.open(path)
        indx=1
        trying=True
        while trying:
            try:
                Z = src.read(indx)[::10,::10]
                feature_info = self._layer('{} - Band {}'.format(name,indx),show=show)
                opc = 0.8
                colormap = linear._colormaps['viridis'].scale(Z.min(),Z.max())
                colormap.caption = '{} - Band {}'.format(name,indx)
                Zc = np.zeros((Z.shape[0],Z.shape[1],4))
                for ii in range(Z.shape[0]):
                    for jj in range(Z.shape[1]):
                        z = Z[ii,jj]
                        if np.isnan(z):
                            Zc[ii,jj,:] = 0.0
                        else:
                            Zc[ii,jj,:] = list(colormap.rgba_floats_tuple(z))


                img = folium.raster_layers.ImageOverlay(
                    name="Band {}".format(indx),
                    bounds = [[src.bounds[1], src.bounds[0]], [src.bounds[3], src.bounds[2]]],
                    image=Zc,
                    opacity=opc,
                    mercator_project=True,
                    pixelated=False,
                    show = show
                )
                img.add_to(feature_info)
                
                self.map.add_child(colormap)
                self.map.add_child(BindColormap(feature_info,colormap))
            except:
                trying=False
            indx+=1




    def _TimeData(self,geojson,predefined=None,**kwargs):
        p = paramsObject('TimeData',predefined=predefined,**kwargs)

        for ii in range(len(geojson['features'])):
            geojson['features'][ii]['properties']['style'] = {}
            geojson['features'][ii]['properties']['style']["weight"]  = p['weight']
            geojson['features'][ii]['properties']['style']["color"]   = p['color']
            geojson['features'][ii]['properties']['style']["opacity"] = p['opacity']

            geojson['features'][ii]['properties']['icon'] = p["icon"],
            geojson['features'][ii]['properties']['iconstyle'] = {'color': p['color'],'iconSize': [0.1,0.1]}

        TimestampedGeoJson(
            geojson,
            period=p['period'],
            duration=p['duration'],
            auto_play=False,
            loop=False,
            loop_button=True,
            date_options='DD/MM/YYYY HH:mm:ss',
            add_last_point=p['point']
        ).add_to(self.map)
    





# def TimeMapSDA(PATH,map):
#     #'/Users/jsmith/Documents/Research/Researcher_BAS/RoutePlanning/SDADT-Positions'
#     Info = SDAPosition(PATH)
#     lines=[]
#     Points = Info[['Long','Lat']].to_numpy()
#     Points[:,0] = Points[:,0]-360
#     entry = {}
#     entry['coordinates'] = Points.tolist()
#     entry['dates'] = Info['Time'].dt.strftime('%Y-%m-%dT%H:%M:%S').to_list()
#     entry['color'] = 'blue'
#     lines.append(entry)


#     TMS = Info['Time'].dt.strftime('%Y-%m-%dT%H:%M:%S').to_list()
#     pointfeatures = [
#         {
#             "type": "Feature",
#             "geometry": {
#                 "type": "Point",
#                 "coordinates": pt.tolist(),
#             },
#             'properties': {
#                 'time': TMS[idx],
#                 'style': {'color': ''},
#                 'icon': 'circle',
#                 'iconstyle': {
#                     'fillColor': '#black',
#                     'fillOpacity': 0.8,
#                     'stroke': 'true',
#                     'radius': 2
#                 }
#     },

#         }
#         for idx,pt in enumerate(Points)
#     ]


#     features = [
#         {
#             "type": "Feature",
#             "geometry": {
#                 "type": "LineString",
#                 "coordinates": line["coordinates"],
#             },
#             "properties": {
#                 "times": line["dates"],
#                 "style": {
#                     "weight": line["weight"] if "weight" in line else 3,
#                     'color': 'blue',
#                     "line-dasharray": [0.1, 1.8]
#                 },
#                 'icon': 'circle',
#                 'iconstyle': {'color': 'blue','iconSize': [1,1]}
#             },
#         }
#         for line in lines
#     ]

#     features =  features + pointfeatures

#     TimestampedGeoJson(
#         {
#             "type": "FeatureCollection",
#             "features": features,
#         },
#         period="PT1H",
#         duration="P7D",
#         auto_play=False,
#         add_last_point=True,
#         max_speed=50
#     ).add_to(map)
#     return map


# def MapMesh(cellGrid,map,threshold=0.8):
#     DF = MeshDF(cellGrid)
#     LandDF = DF[DF['Land'] == True]
#     IceDF  = DF[DF['Land'] == False]
#     ThickIceDF = IceDF[IceDF['Ice Area'] >= threshold*100]
#     ThinIceDF  = IceDF[IceDF['Ice Area'] < threshold*100]

#     # ==== Plotting Ice ==== 
#     iceInfo = folium.FeatureGroup(name='Ice Mesh')
#     folium.GeoJson(
#         IceDF,
#         style_function=lambda x: {
#                 'fillColor': 'white',
#                 'color': 'gray',
#                 'weight': 0.5,
#                 'fillOpacity': x['properties']['Ice Area']/100
#             }
#     ).add_to(iceInfo)
#     folium.GeoJson(
#         ThickIceDF,
#         style_function=lambda x: {
#                 'color': 'red',
#                 'weight': 0.5,
#                 'fillOpacity': 0.0
#             }
#     ).add_to(iceInfo)
#     iceInfo.add_to(map)

#     # ===== Plotting Land =====
#     landInfo = folium.FeatureGroup(name='Land Mesh')
#     folium.GeoJson(
#         LandDF,
#         style_function= lambda x: {
#                 'fillColor': 'green',
#                 'color': 'gray',
#                 'weight': 0.5,
#                 'fillOpacity': 0.3
#             }
#     ).add_to(landInfo)
#     landInfo.add_to(map)

#     # ===== Plotting Mesh Info =====
#     bathInfo = folium.FeatureGroup(name='Bathymetry Mesh',show=False)
#     colormap = linear.Reds_09.scale(min(ThinIceDF['Depth']),max(ThinIceDF['Depth']))
#     folium.GeoJson(
#         IceDF,
#         style_function=lambda x: {
#                 'fillColor': colormap(x['properties']['Depth']),
#                 'color': 'gray',
#                 'weight': 0.5,
#                 'fillOpacity': 0.3
#             }
#     ).add_to(bathInfo)
#     bathInfo.add_to(map)
#     # ===== Plotting Mesh Info =====
#     meshInfo = folium.FeatureGroup(name='Mesh Information',show=False)
#     folium.GeoJson(
#         DF,
#         style_function= lambda x: {
#                 'fillColor': 'black',
#                 'color': 'gray',
#                 'weight': 0.5,
#                 'fillOpacity': 0.
#             },
#         tooltip=folium.GeoJsonTooltip(
#             fields=['Ice Area', 'Land','Cx','Cy','Depth','Vector','Index'],
#             aliases=['Ice Area (%)', 'Land','Centroid Cx [Long]','Centroid Cy [Lat]','Depth(m)','Vector (m/s)','Cell Index'],
#             localize=True
#         ),
#         name='Land Grid'
#     ).add_to(meshInfo)
#     meshInfo.add_to(map)
#     return map