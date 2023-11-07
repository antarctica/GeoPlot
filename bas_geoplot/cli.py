import argparse
import json
import inspect
import logging
import numpy as np
import pandas as pd

from bas_geoplot import __version__ as version
from bas_geoplot.utils import setup_logging, timed_call, gpx_route_import
from bas_geoplot.interactive import Map


@setup_logging
def get_args(default_output: str):
    """
        Add required command line arguments for all CLI entry points.
    """

    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default=default_output, help="Output file")
    ap.add_argument("-v", "--verbose", default=False, action="store_true",
                    help="Turn on DEBUG level logging")
    ap.add_argument("-s", "--static", default=False, action="store_true",
                    help="Save the plot as a static .PNG")
    ap.add_argument("-c", "--currents_paths", default='', help="Path to currents file")
    ap.add_argument("-l", "--coastlines", default='', help="Loading Offline Coastlines")
    ap.add_argument("-j", "--offline_filepath", default='', help="Location of Offline File Information")
    ap.add_argument("-p", "--plot_sectors", default=False, action="store_true",
                    help="Plot array values as separate polygons")
    ap.add_argument("mesh", type=argparse.FileType('r'), help="file location of mesh to be plot")
    ap.add_argument('--version', action='version',
                    version='%(prog)s {version}'.format(version=version))
    ap.add_argument("-t", "--rm_titlebar", default=False, action="store_true",
                    help="Remove titlebar from html")
    ap.add_argument("-r", "--route", default=None, help="Plot additional route on mesh")
    ap.add_argument("-a", "--arrows", default=False, action="store_true",
                    help="Add directional arrows to routes")
    ap.add_argument("--custom_title", default="", help="Add a custom title to the plot")

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
    region = info['config']['mesh_info']['region']
    split_level = info['config']['mesh_info']['splitting']['split_depth']

    # Set-up title bar
    if args.rm_titlebar:
        output = None
    else:
        if args.custom_title:
            output = '{} | Start Date: {}, End Date: {} | Split level: {}'.format(args.custom_title,
                                                                                  region['start_time'],
                                                                                  region['end_time'], split_level)
        else:
            output = ' '.join(args.output.split('/')[-1].split('.')[:-1])
            output = '{} | Start Date: {}, End Date: {} | Split level: {}'.format(output, region['start_time'],
                                                                                  region['end_time'], split_level)

    # Put mesh bounds in format required by fit_to_bounds
    mesh_bounds = [[region["lat_min"], region["long_min"]], [region["lat_max"], region["long_max"]]]

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
            currents = currents[(currents['cx'] >=  info['config']['mesh_info']['region']['long_min']) &
                                (currents['cx'] <=  info['config']['mesh_info']['region']['long_max']) &
                                (currents['cy'] >=  info['config']['mesh_info']['region']['lat_min']) &
                                (currents['cy'] <=  info['config']['mesh_info']['region']['lat_max'] )
            ].reset_index(drop=True)
            mp.Vectors(currents,'Currents - Raw Data', show=False, predefined='Currents')
        mp.Vectors(mesh,'Currents - Mesh', show=False, predefined='Currents')
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
        mp.Paths(paths, 'Route - Traveltime', predefined='Traveltime (Days)', arrows=args.arrows)
        mp.Paths(paths, 'Route - Distance', predefined='Distance (Nautical miles)', show=False,
                 arrows=args.arrows)
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
        # Read in as either GPX or GeoJSON
        filetype = args.route.split('.')[-1]
        if filetype == 'gpx':
            route_json = gpx_route_import(args.route)
            
            mp.Paths(route_json, 'User Route - GPX', predefined='black', arrows=args.arrows)
        elif filetype in ['json', 'geojson']:
            route_json = json.load(open(args.route))
            mp.Paths(route_json, 'User Route - Traveltime', predefined='black', arrows=args.arrows)
            mp.Paths(route_json, 'User Route - Fuel', predefined='green', show=False, arrows=args.arrows)
        else:
            raise NameError("User defined route needs to be GPX or GeoJSON file!")


    # Set-up mesh info and save map to html file
    mp.MeshInfo(mesh, 'Mesh Info', show=False)
    mp.fit_to_bounds(mesh_bounds)
    logging.info('Saving plot to {}'.format(args.output))
    mp.save(args.output)
