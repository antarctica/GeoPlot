import argparse
import json
import inspect
import logging
import numpy as np
import pandas as pd

from bas_geoplot import __version__ as version
from bas_geoplot.utils import setup_logging, timed_call
from bas_geoplot.interactive import Map


@setup_logging
def get_args(default_output: str):
    """
        Add required command line arguments for all CLI entry points.
    """

    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default=default_output, help="Output file")
    ap.add_argument("-v", "--verbose", default=False, action="store_true", help="Turn on DEBUG level logging")
    ap.add_argument("-s", "--static", default=False, action="store_true", help="Save the plot as a static .PNG")
    ap.add_argument("-c", "--currents_paths", default='', help="Path to currents file")
    ap.add_argument("-l", "--coastlines", default='', help="Loading Offline Coastlines")
    ap.add_argument("-j", "--offline_filepath", default='', help="Location of Offline File Information")
    ap.add_argument("-p", "--plot_sectors", default=False, action="store_true",
                    help="Plot array values as separate polygons")
    ap.add_argument("mesh", type=argparse.FileType('r'), help="file location of mesh to be plot")
    ap.add_argument('--version', action='version',
                    version='%(prog)s {version}'.format(version=version))
    ap.add_argument("-t", "--rm_titlebar", default=False, action="store_true", help="Remove titlebar from html")
    ap.add_argument("-r", "--route", default=None, help="Plot additional route on mesh")
    ap.add_argument("-a", "--arrows", default=False, action="store_true", help="Add directional arrows to routes")

    return ap.parse_args()


@timed_call
def plot_mesh_cli():
    """
        CLI entry point to plot an environmental mesh and associated routes/waypoints
    """

    # Set output location and load mesh info
    args = get_args("interactive_plot.html")
    logging.info("{} {}".format(inspect.stack()[0][3][:-4], version))
    info = json.load(args.mesh)
    mesh = pd.DataFrame(info['cellboxes'])
    region = info['config']['Mesh_info']['Region']

    # Set-up title bar
    if args.rm_titlebar:
        output = None
    else:
        output = ' '.join(args.output.split('/')[-1].split('.')[:-1])
        output = '{} | Start Date: {}, End Date: {}'.format(output, region['startTime'], region['endTime'])

    # Put mesh bounds in format required by fit_to_bounds
    mesh_bounds = [[region["latMin"], region["longMin"]], [region["latMax"], region["longMax"]]]

    # Initialise Map object
    if args.offline_filepath != '':
        logging.debug("Offline .js & .css datastore - {}".format(args.offline_filepath))
        if args.coastlines != '':
            mp = Map(title=output, offline_coastlines=args.coastlines, offline_filepath=args.offline_filepath)
        else:
            mp = Map(title=output, offline_filepath=args.offline_filepath)
    else:
        if args.coastlines != '':
            mp = Map(title=output, offline_coastlines=args.coastlines)
        else:
            mp = Map(title=output)

    # Plot maps of mesh info
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
        mp.Maps(mesh, 'Fuel', predefined='Fuel (Tonnes/Day)', show=False, plot_sectors=args.plot_sectors)
        mp.Maps(mesh, 'tCO2e', predefined='tCO2e', show=False, plot_sectors=args.plot_sectors)
    if 'battery' in mesh.columns:
        logging.debug('Plotting battery usage')
        mp.Maps(mesh, 'Battery Usage', predefined='Battery Usage', show=False, plot_sectors=args.plot_sectors)
    if 'speed' in mesh.columns:
        logging.debug('Plotting vessel maximum speed')
        mp.Maps(mesh, 'Max Speed', predefined='Max Speed (knots)', show=False,plot_sectors=args.plot_sectors)
    if ('uC' in mesh.columns) and ('vC' in mesh.columns):
        mesh['mC'] = np.sqrt(mesh['uC']**2 + mesh['vC']**2)
        logging.debug('Plotting currents')
        if args.currents_paths != '':
            logging.debug('Plotting currents from file')
            currents = pd.read_csv(args.currents_paths)
            currents = currents[(currents['cx'] >=  info['config']['Mesh_info']['Region']['longMin']) &
                                (currents['cx'] <=  info['config']['Mesh_info']['Region']['longMax']) &
                                (currents['cy'] >=  info['config']['Mesh_info']['Region']['latMin']) &
                                (currents['cy'] <=  info['config']['Mesh_info']['Region']['latMax'] )
            ].reset_index(drop=True)
            mp.Vectors(currents,'Currents - Raw Data', show=False, predefined='Currents')
        mp.Vectors(mesh,'Currents - Mesh', show=False, predefined='Currents')
    if ('u10' in mesh.columns) and ('v10' in mesh.columns):
        mesh['m10'] = np.sqrt(mesh['u10'] ** 2 + mesh['v10'] ** 2)
        mp.Vectors(mesh, 'Winds', predefined='Winds')
        logging.debug('Plotting winds')
    if 'swh' in mesh.columns:
        logging.debug("Plotting Wave Height")
        mp.Maps(mesh, 'Wave Height', predefined='Wave Height')
    if ('uW' in mesh.columns) and ('vW' in mesh.columns):
        mp.Vectors(mesh, 'Wave Direction', predefined='Wave Direction', show=False)
        logging.debug('Plotting wave direction')
    if 'ext_waves' in mesh.columns:
        logging.debug("Plotting Extreme Wave areas")
        mp.Maps(mesh, 'Extreme Waves', predefined='Extreme Waves')

    # Plot routes and waypoints
    if 'paths' in info.keys():
        logging.debug('Plotting paths')
        paths = info['paths']
        mp.Paths(paths, 'Route - Traveltime', predefined='Traveltime (Days)', arrows=args.arrows)
        mp.Paths(paths, 'Route - Distance', predefined='Distance (Nautical miles)', show=False, arrows=args.arrows)
        mp.Paths(paths, 'Route - Max Speed', predefined='Max Speed (knots)', show=False, arrows=args.arrows)
        if 'fuel' in mesh.columns:
            mp.Paths(paths, 'Route - Fuel', predefined='Fuel', show=False, arrows=args.arrows)
            mp.Paths(paths, 'Route - tCO2e', predefined='tCO2e', show=False, arrows=args.arrows)
        if 'battery' in mesh.columns:
            mp.Paths(paths, 'Route - Battery', predefined='Battery', show=False, arrows=args.arrows)
    if 'waypoints' in info.keys():
        logging.debug('Plotting waypoints')
        waypoints = pd.DataFrame(info['waypoints'])
        mp.Points(waypoints, 'Waypoints', names={"font_size":10.0})
    if args.route:
        logging.debug('Plotting user defined route')
        with open(args.route, "r") as f:
            route_json = json.load(f)
        mp.Paths(route_json, 'User Route - Traveltime', predefined='Traveltime (Days)', arrows=args.arrows)
        mp.Paths(route_json, 'User Route - Fuel', predefined='Fuel', show=False, arrows=args.arrows)

    # Set-up mesh info and save map to html file
    mp.MeshInfo(mesh, 'Mesh Info', show=False)
    mp.fit_to_bounds(mesh_bounds)
    logging.info('Saving plot to {}'.format(args.output))
    mp.save(args.output)
