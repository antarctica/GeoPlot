import matplotlib
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.collections import LineCollection


import os,json, copy
import pandas as pd
import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt

from shapely import wkt
import geopandas as gpd
from branca.colormap import linear



def paramsObject(layer,predefined=None,**kwargs):
    # Loading config of standards
    p_file = os.path.join(os.path.dirname(__file__),'config','static.json')
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

class Map:
    def __init__(self, predefined=None,**kwargs):
    
        self.zorder = 0 

        # ==== Loading standard configs
        p = paramsObject('Basemap',predefined=predefined,**kwargs)

        if p['title']:
            title='{}'.format(p['title'])
        else:
            title=''


        #Initialising the basemap
        if 'CRS' in p.keys():
            if p['CRS'] == 'Mercartor':
                self.ccrs = ccrs.Mercator()
            if p['CRS'] == 'Orthographic':
                self.ccrs = ccrs.Orthographic(central_longitude=p['pole_centre'][0],central_latitude=p['pole_centre'][1])
        else:
           self.ccrs = ccrs.Mercator()
  
        self.fig = plt.figure(figsize=(p['figure_size'][0],p['figure_size'][1]))
        matplotlib.rcParams.update({'font.size': 16})

        self.ax = plt.axes(projection=self.ccrs)
        self._bounds = p["bounds"]

        self.ax.set_extent([p["bounds"][0][0]+1e-6,p["bounds"][1][0]-1e-6,
                            p["bounds"][0][1]+1e-6,p["bounds"][1][1]-1e-6], crs=ccrs.PlateCarree())
        self.ax.add_image(cimgt.GoogleTiles(), 3)
        self.ax.coastlines(resolution='50m')
        self.ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,linewidth=0.5,linestyle='--')
        self.ax.add_feature(cfeature.BORDERS)
        plt.title(r'{}'.format(title),fontsize=14,loc='left',color='blue')


    def Points(self,dataframe_points,predefined=None,**kwargs):

        p = paramsObject('Points',predefined=predefined,**kwargs)
        self.zorder+=1

        if type(p['color']) is str:
            self.ax.scatter(dataframe_points['Long'],dataframe_points['Lat'],p["marker_size"],marker=p['icon'],transform=ccrs.PlateCarree(),color=p['color'],zorder=self.zorder)

        if type(p['color']) is dict:
            colormap = linear._colormaps[p['color']['colormap']].scale(dataframe_points[p['color']['data_name']].min(), 
                                                                            dataframe_points[p['color']['data_name']].max())
            colours = [colormap.rgba_floats_tuple(ii) for ii in np.array(dataframe_points[p['color']['data_name']])]
            self.ax.scatter(dataframe_points['Long'],dataframe_points['Lat'],p["marker_size"],marker=p['icon'],transform=ccrs.PlateCarree(),color=colours,zorder=self.zorder)

        if type(p['names']) is dict:
            transform = ccrs.PlateCarree()._as_mpl_transform(self.ax)


            dataframe_points = dataframe_points[(dataframe_points['Long'] >= self._bounds[0][0]+1e-6) &\
                                                (dataframe_points['Long'] <= self._bounds[1][0]+1e-6) &\
                                                (dataframe_points['Lat'] >= self._bounds[0][1]+1e-6) &\
                                                (dataframe_points['Lat'] <= self._bounds[1][1]+1e-6)]

            for idx,row in dataframe_points.iterrows():
                self.ax.annotate(row['Name'], xy=(row['Long']+p['names']['dX'], row['Lat']+p['names']['dY']),xycoords=transform,color=p['names']['color'], size=p['names']['Size'])


    def Paths(self,geojson,predefined=None,**kwargs):
        '''
            Plotting a paths type object
        '''

        p = paramsObject('Paths',predefined=predefined,**kwargs)
        self.zorder+=1

        # Defining the feature groups to add
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

            start_wpt = path['properties']['from']
            end_wpt   = path['properties']['to']

            if p['data_name']:
                data_val = np.array(path['properties'][p['data_name']])
            else:
                data_val = np.array(len(points))


            if type(p['line_color']) is dict:
                # Add ColorLines
                x = self.ccrs.transform_points(x=points[:,0], y=points[:,1],
                                                src_crs=ccrs.PlateCarree())
                xcs = np.array([x[:,0],x[:,1]]).T.reshape(-1,1,2) 
                segments = np.concatenate([xcs[:-1], xcs[1:]], axis=1)   
                lc = LineCollection(segments, cmap=p["line_color"]['color'], linewidth=p['line_width'],norm=plt.Normalize(vmin=min_val, vmax=max_val),zorder=self.zorder) 
                lc.set_array(data_val)                                           
                self.ax.add_collection(lc) 
                
            
            else:
                self.ax.plot(points[:,0],points[:,1],transform=ccrs.PlateCarree(),linewidth=p['line_width'],color=p["line_color"],alpha=p["line_opacity"],zorder=self.zorder)

            if p['path_points']:
                self.ax.scatter(points[:,0],points[:,1],15,color=p["line_color"])
        
        if type(p['line_color']) is dict:
            cbaxes = inset_axes(self.ax, '25%', '3%', loc =1)
            cbaxes.set_facecolor([1,1,1,0.7])
            cb=self.fig.colorbar(lc,cax=cbaxes,orientation='horizontal',label='{} ({},Max Value={:.3f})'.format(p['data_name'],p['unit'],max_val)) #make colorbar


    def Maps(self,dataframe_pandas,predefined=None,**kwargs):
        '''
            Plotting a map type object
        '''

        p = paramsObject('Maps',predefined=predefined,**kwargs)
        self.zorder+=1

        dataframe_pandas = copy.copy(dataframe_pandas)
        dataframe_pandas['geometry'] = dataframe_pandas['geometry'].apply(wkt.loads)
        dataframe_geo = gpd.GeoDataFrame(dataframe_pandas,crs='EPSG:4326', geometry='geometry')

        if p['data_name'] and p['trim_min']:
            dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']] > p['trim_min']]
        if p['data_name'] and p['trim_max']:
            dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']] < p['trim_max']]
        if p['data_name'] and (type(dataframe_geo[p['data_name']].iloc[0].item()) is bool):
            dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']] == True]

        if (type(p['fill_color']) is dict) and (p['data_name']):
            dataframe_geo = dataframe_geo[dataframe_geo[p['data_name']].notna() & ~np.isinf(abs(dataframe_geo[p['data_name']]))]
            cmin = dataframe_geo[p['data_name']].min()
            cmax = dataframe_geo[p['data_name']].max()

            colormap = linear._colormaps[p['fill_color']['colormap']].scale(cmin,cmax)

            for _,poly in dataframe_geo.iterrows():
                data = poly[p['data_name']]
                if not np.isnan(data):
                    colour = colormap.rgba_floats_tuple(data)
                    self.ax.add_geometries([poly['geometry']], crs=ccrs.PlateCarree(), edgecolor=p["line_color"], facecolor=colour,alpha=p["fill_opacity"],lw=p["line_width"],zorder=self.zorder)
            


        else:
            if (p["fill_opacity"] == 0.0):
                for _,poly in dataframe_geo.iterrows():
                    x,y = poly['geometry'].exterior.coords.xy
                    self.ax.plot(np.array(x),np.array(y),color=p["line_color"],linewidth=p["line_width"],transform=ccrs.PlateCarree(),zorder=self.zorder)
            else:
                for _,poly in dataframe_geo.iterrows():
                    self.ax.add_geometries([poly['geometry']], crs=ccrs.PlateCarree(), edgecolor=p["line_color"], facecolor=p['fill_color'],alpha=p["fill_opacity"],lw=p["line_width"],zorder=self.zorder)   


    #     if type(p['line_color']) is dict:

 
    def Vectors(self,dataframe_pandas,predefined=None,**kwargs):
        '''
        
        '''
        p = paramsObject('Vectors',predefined=predefined,**kwargs)
        self.zorder+=1

        dataframe_pandas = copy.copy(dataframe_pandas)

        self.ax.quiver(dataframe_pandas[p['Long']].to_numpy(),dataframe_pandas[p['Lat']].to_numpy(),dataframe_pandas[p['U']].to_numpy(),dataframe_pandas[p['V']].to_numpy(),zorder=self.zorder,transform=ccrs.PlateCarree(),width=p['line_width'],scale=p['scale'])

    def savefig(self,filename,**kwargs):
        plt.savefig(filename,**kwargs)
    
    def show(self):
        plt.show()