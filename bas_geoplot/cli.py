import argparse
import json
import inspect
import logging
import numpy as np

from bas_geoplot import __version__ as version
from bas_geoplot.utils import setup_logging, timed_call

@setup_logging
def get_args(default_output: str):
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output",default=default_output,help="Output file")
    ap.add_argument("-v", "--verbose",default=False,action="store_true",help="Turn on DEBUG level logging")
    ap.add_argument("-s", "--static",default=False,action="store_true",help="Save the plot as a static .PNG")           
    ap.add_argument("-c", "--currents_paths",default='',help="Path to currents file")
    ap.add_argument("-l", "--coastlines",default='',help="Loading Offline Coastlines")
    ap.add_argument("-j", "--offline_filepath",default='',help="Location of Offline File Information")
    ap.add_argument("-p", "--plot_sectors",default=False,action="store_true",help="Plot array values as separate polygons")
    ap.add_argument("mesh", type=argparse.FileType('r'),help="file location of mesh to be plot")
    return ap.parse_args()

@timed_call
def plot_mesh_cli():
    from bas_geoplot.interactive import Map
    import pandas as pd

    args = get_args("interactive_plot.html")
    logging.info("{} {}".format(inspect.stack()[0][3][:-4], version))
    info = json.load(args.mesh)
    mesh = pd.DataFrame(info['cellboxes'])

    output = ' '.join(args.output.split('/')[-1].split('.')[:-1])
    output = '{} | Start Date: {}, End Date: {}'.format(output, info['config']['Mesh_info']['Region']['startTime'], info['config']['Mesh_info']['Region']['endTime'])


    if args.offline_filepath != '':
        logging.debug("offline .js & .css datastore - {}".format(args.offline_filepath))
        if args.coastlines != '':
            mp = Map(title=output,offline_coastlines=args.coastlines,offline_filepath=args.offline_filepath)
        else:
            mp = Map(title=output,offline_filepath=args.offline_filepath)
    else:
        if args.coastlines != '':
            mp = Map(title=output,offline_coastlines=args.coastlines)
        else:
            mp = Map(title=output)


    if 'SIC' in mesh.columns:
        logging.debug("plotting Sea Ice Concentration")
        mp.Maps(mesh,'SIC',predefined='SIC')
    if 'ext_ice' in mesh.columns:
        logging.debug("plotting Extreme Ice areas")
        mp.Maps(mesh,'Extreme Ice',predefined='Extreme Sea Ice Conc')
    if 'land' in mesh.columns:
        logging.debug("plotting Land Mask")
        mp.Maps(mesh,'Land Mask',predefined='Land Mask')
    if 'fuel' in mesh.columns:
        logging.debug('plotting Fuel usage per day and tCO2e')
        mp.Maps(mesh,'Fuel',predefined='Fuel (Tonnes/Day)',show=False,plot_sectors=args.plot_sectors)
        mp.Maps(mesh,'tCO2e',predefined='tCO2e',show=False,plot_sectors=args.plot_sectors)
    if 'speed' in mesh.columns:
        logging.debug('plotting vessel maximum speed')
        mp.Maps(mesh,'Max Speed',predefined='Max Speed (knots)',show=False,plot_sectors=args.plot_sectors)
    if ('uC' in mesh.columns) and ('vC' in mesh.columns):
        mesh['mC'] = np.sqrt(mesh['uC']**2 + mesh['vC']**2)
        logging.debug('plotting currents')
        if args.currents_paths != '':
            logging.debug('plotting currents from file')
            currents = pd.read_csv(args.currents_paths)
            currents = currents[(currents['cx'] >=  info['config']['Mesh_info']['Region']['longMin']) & (currents['cx'] <=  info['config']['Mesh_info']['Region']['longMax']) & (currents['cy'] >=  info['config']['Mesh_info']['Region']['latMin']) & (currents['cy'] <=  info['config']['Mesh_info']['Region']['latMax'] )].reset_index(drop=True)
            mp.Vectors(currents,'Currents - Raw Data', show=False, predefined='Currents')
        mp.Vectors(mesh,'Currents - Mesh', show=False, predefined='Currents')

    if 'paths' in info.keys():
        logging.debug('plotting paths')
        paths = info['paths']
        mp.Paths(paths,'Routes - Traveltimes',predefined='Traveltime (Days)')
        mp.Paths(paths,'Routes - Distance',predefined='Distance (Nautical miles)',show=False)
        mp.Paths(paths,'Routes - Max Speed',predefined='Max Speed (knots)',show=False)
        mp.Paths(paths,'Routes - Fuel',predefined='Fuel',show=False)
        mp.Paths(paths,'Routes - tCO2e',predefined='tCO2e',show=False)
    if 'waypoints' in info.keys():
        logging.debug('plotting waypoints')
        waypoints = pd.DataFrame(info['waypoints'])
        mp.Points(waypoints,'Waypoints',names={"font_size":10.0})
    mp.MeshInfo(mesh,'Mesh Info',show=False)
    logging.info('Saving plot to {}'.format(args.output))
    mp.save(args.output)




