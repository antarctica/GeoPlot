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
    def __init__(self,Title,predefined=None,**kwargs):


        '''
            Initialising a map instance
        '''

        # ==== Loading standard configs
        p = paramsObject('Basemap',predefined=predefined,**kwargs)

        title_html = '''
            <h1 style="color:#003b5c;font-size:16px">
            &ensp;<img src='https://i.ibb.co/JH2zknX/Small-Logo.png' alt="BAS-colour-eps" border="0" style="width:40px;height:40px;"> 
            <img src="https://i.ibb.co/XtZdzDt/BAS-colour-eps.png" alt="BAS-colour-eps" border="0" style="width:179px;height:40px;"> 
            &ensp;|&ensp; RoutePlanner &ensp;|&ensp; {}
            </h1>
            </body>
            '''.format(Title)   
        self.map = folium.Map(location=p['map_centre'],zoom_start=p['zoom_start'],tiles=None)
        
        bsmap = folium.FeatureGroup(name='BaseMap')
        folium.TileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/NatGeo_World_Map/MapServer/tile/{z}/{y}/{x}.png',attr="toner-bcg", name='Basemap').add_to(bsmap)
        bsmap.add_to(self.map)
        bsmap = folium.FeatureGroup(name='Dark BaseMap',show=False)
        folium.TileLayer(tiles="https://{s}.basemaps.cartocdn.com/rastertiles/dark_all/{z}/{x}/{y}.png",attr='&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>',name='darkmatter').add_to(bsmap)
        bsmap.add_to(self.map)

        if p['plot_title']:
            self.map.get_root().html.add_child(folium.Element(title_html))

    def Paths(self,geojson,name,show=True,predefined=None,**kwargs):
        

        p = paramsObject('Paths',predefined=predefined,**kwargs)

        # Defining the feature groups to add
        pths        = folium.FeatureGroup(name='{}'.format(name),show=show)
        if p['path_points']:
            pths_points = folium.FeatureGroup(name='{} - Path Points'.format(name),show=show)
            
        paths = geojson['features']

        # Determining min-max values of all paths if colormap being used
        if type(p['line_color']) is dict:
            max_val = -np.inf
            min_val = np.inf
            for path in copy.deepcopy(paths):
                if np.array(path['properties'][p['data_name']]).max() > max_val:
                    max_val = np.array(path['properties'][p['data_name']]).max()

                if np.array(path['properties'][p['data_name']]).min() < min_val:
                    min_val = np.array(path['properties'][p['data_name']]).min()

        # Determining max travel-times of all paths
        for path in paths:
            points   = np.array(path['geometry']['coordinates'])
            if p['data_name']:
                data_val = np.array(path['properties'][p['data_name']])
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

            else:
                folium.PolyLine(points,color=p['line_color'], weight=p['line_width'], opacity=p['line_opacity']).add_to(pths)


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
                    ).add_to(pths_points)
        
        if p['path_points']:
            pths_points.add_to(self.map)
        if type(p['line_color']) is dict:
            self.map.add_child(pths)
            self.map.add_child(colormap)
            self.map.add_child(BindColormap(pths,colormap))
        else:
            pths.add_to(self.map)


    def Points(self,dataframe_points,name,show=True,predefined=None,**kwargs):

        p = paramsObject('Points',predefined=predefined,**kwargs)

        wpts      = folium.FeatureGroup(name='{}'.format(name),show=show)

        if p['names']:
            wpts_name = folium.FeatureGroup(name='{} - Names'.format(name),show=show)


        for id,wpt in dataframe_points.iterrows():
            loc = [wpt['Lat'], wpt['Long']]
            folium.Marker(
                location=loc,
                icon=plugins.BeautifyIcon(icon='circle',
                                                border_color='transparent',
                                                background_color='transparent',
                                                border_width=p['line_width'],
                                                text_color=p['color'],
                                                inner_icon_style='margin:0px;font-size:0.8em'),
                popup="<b>{} [{:4f},{:4f}]</b>".format(wpt['Name'],loc[0],loc[1]),
            ).add_to(wpts)    

            if p['names']:

                folium.Marker(
                            location=loc,
                                icon=folium.features.DivIcon(
                                    icon_size=(250,36),
                                    icon_anchor=(0,0),
                                    html='<div style="font-size: {}pt">{}</div>'.format(p['names']['font_size'],wpt['Name']),
                                    ),
                ).add_to(wpts_name)


        wpts.add_to(self.map)
        if p['names']:
            wpts_name.add_to(self.map)


    def Maps(self,dataframe_pandas,name,show=True,predefined=None,**kwargs):
        '''
            Plotting a map type object
        '''

        p = paramsObject('Maps',predefined=predefined,**kwargs)
        self.p = p
        dataframe_pandas['geometry'] = dataframe_pandas['geometry'].apply(wkt.loads)
        dataframe_geo = gpd.GeoDataFrame(dataframe_pandas,crs='EPSG:4326', geometry='geometry')

        feature_info = folium.FeatureGroup(name='{}'.format(name),show=show)


        if p['data_name']:
            try:
                dataframe_geo[p['data_name']]
            except:
                print('Dataname not in variables')



        if p['data_name'] and p['trim_min']:
            dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']] > p['trim_min']]
        if p['data_name'] and p['trim_max']:
            dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']] < p['trim_max']]
            


        if (type(p['fill_color']) is dict) and (p['data_name']):
            dataframe_geo
            dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']].notna()]


            cmin = dataframe_geo[p['data_name']].min()
            cmax = dataframe_geo[p['data_name']].max()

            colormap = linear._colormaps[p['fill_color']['colormap']].scale(cmin,cmax)
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
            self.map.add_child(feature_info)
            self.map.add_child(colormap)
            self.map.add_child(BindColormap(feature_info,colormap))
        else:
            if p['data_name'] and (type(dataframe_geo[p['data_name']].iloc[0].item()) is bool):
                dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']] == True]
                
            folium.GeoJson(
                dataframe_geo,
                style_function=lambda x: {
                    'fillColor': p['fill_color'],
                    'color': p['line_color'],
                    'weight': p['line_width'],
                    'fillOpacity': p['fill_opacity']
                    }
            ).add_to(feature_info)
            self.map.add_child(feature_info)



    def Vectors(self,Currents,name,show=True,predefined=None,**kwargs):



        vcts = folium.FeatureGroup(name='Currents',show=show)
        for idx,vec in Currents.iterrows():
            loc =[[vec['Y'],vec['X']],[vec['Y']+vec['V']*scale,vec['X']+vec['U']*scale]]
            folium.PolyLine(loc, color="gray",weight=1.4).add_to(vcts)
            # get pieces of the line
            pairs = [(loc[idx], loc[idx-1]) for idx, val in enumerate(loc) if idx != 0]
            # get rotations from forward azimuth of the line pieces and add an offset of 90°
            geodesic = Geod(ellps='WGS84')
            rotations = [geodesic.inv(pair[0][1], pair[0][0], pair[1][1], pair[1][0])[0]+90 for pair in pairs]
            # create your arrow
            for pair, rot in zip(pairs, rotations):
                folium.RegularPolygonMarker(location=pair[0], color='gray', fill=True, fill_color='gray', fill_opacity=1,
                                            number_of_sides=3, rotation=rot,radius=2,weight=0.8).add_to(vectors)

        vcts.add_to(self.map)



    def MeshInfo(self,cellboxes,name,show=True):
        dataframe_pandas = pd.DataFrame(cellboxes)
        dataframe_pandas['geometry'] = dataframe_pandas['geometry'].apply(wkt.loads)
        dataframe_geo = gpd.GeoDataFrame(dataframe_pandas,crs='EPSG:4326', geometry='geometry')


        feature_info = folium.FeatureGroup(name='{}'.format(name),show=show)
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
        feature_info.add_to(self.map)

    def TimeData(self,geojson,predefined=None,**kwargs):
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
    


    def Show(self):
        folium.LayerControl(collapsed=True).add_to(self.map)
        return self.map


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


# def MapCurrents(cellGrid,map,show=False,scale=15):
#     import folium
#     from pyproj import Geod
#     def bearing(st,en):
#         import numpy as np
#         long1,lat1 = st
#         long2,lat2 = en
#         dlong = long2-long1
#         dlat  = lat2-lat1
#         vector_1 = [0, 1]
#         vector_2 = [dlong, dlat]
#         if np.linalg.norm(vector_2) == 0:
#             return np.nan
#         unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
#         unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
#         dot_product = np.dot(unit_vector_1, unit_vector_2)
#         angle = np.arccos(dot_product)/(np.pi/180)*np.sign(vector_2[0])
#         if (angle==0) & (np.sign(dlat)==-1):
#             angle=180
#         if angle < 0:
#             angle = angle +360
#         angle
#         return angle

#     cellGrid
#     X=[];Y=[];U=[];V=[];
#     for ii in range(len(cellGrid.cellBoxes)):
#         cellbox = cellGrid.cellBoxes[ii]
#         if not isinstance(cellbox, CellBox):
#             continue

#         X.append(cellbox.cx)
#         Y.append(cellbox.cy)
#         U.append(cellbox.getuC())
#         V.append(cellbox.getvC())
#     Currents = pd.DataFrame({'X':X,'Y':Y,'U':U,'V':V})
#     Currents = Currents.dropna()
#     Currents['X'] = Currents['X']


#     vectors = folium.FeatureGroup(name='Currents',show=show)
#     for idx,vec in Currents.iterrows():
#         loc =[[vec['Y'],vec['X']],[vec['Y']+vec['V']*scale,vec['X']+vec['U']*scale]]
#         folium.PolyLine(loc, color="gray",weight=1.4).add_to(vectors)
#         # get pieces of the line
#         pairs = [(loc[idx], loc[idx-1]) for idx, val in enumerate(loc) if idx != 0]
#         # get rotations from forward azimuth of the line pieces and add an offset of 90°
#         geodesic = Geod(ellps='WGS84')
#         rotations = [geodesic.inv(pair[0][1], pair[0][0], pair[1][1], pair[1][0])[0]+90 for pair in pairs]
#         # create your arrow
#         for pair, rot in zip(pairs, rotations):
#             folium.RegularPolygonMarker(location=pair[0], color='gray', fill=True, fill_color='gray', fill_opacity=1,
#                                         number_of_sides=3, rotation=rot,radius=2,weight=0.8).add_to(vectors)

#     vectors.add_to(map)
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


# def BaseMap(location=[-58,-63.7],logo=True,logoPos=[5,88]):
#     map = folium.Map(location=location,zoom_start=2.6,tiles=None)
#     bsmap = folium.FeatureGroup(name='BaseMap')
#     folium.TileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/NatGeo_World_Map/MapServer/tile/{z}/{y}/{x}.png',attr="toner-bcg", name='Basemap').add_to(bsmap)
#     bsmap.add_to(map)
#     if logo:
#         folium.plugins.FloatImage('https://i.ibb.co/JH2zknX/Small-Logo.png',bottom=logoPos[1],left=logoPos[0]).add_to(map)
#     return map

# def LayerControl(map,collapsed=True):
#     folium.LayerControl(collapsed=collapsed).add_to(map)
#     return map