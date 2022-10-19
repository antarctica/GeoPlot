import argparse
import json
import inspect
import logging

from bas_geoplot import __version__ as version
from bas_geoplot.utils import setup_logging, timed_call

@setup_logging
def get_args(default_output: str):
    ap = argparse.ArgumentParser()

    ap.add_argument("-o", "--output",
                    default=default_output,
                    help="Output file")
    ap.add_argument("-v", "--verbose",
                    default=False,
                    action="store_true",
                    help="Turn on DEBUG level logging")
    ap.add_argument("-s", "--static",
                    default=False,
                    action="store_true",
                    help="Save the plot as a static .PNG")           

    ap.add_argument("mesh", type=argparse.FileType('r'),
                    help="file location of mesh to be plot")

    return ap.parse_args()

@timed_call
def plot_mesh_cli():
    from bas_geoplot.interactive import Map
    import pandas as pd

    args = get_args("interactive_plot.html")
    logging.info("{} {}".format(inspect.stack()[0][3][:-4], version))

    info = json.load(args.mesh)

    mesh = pd.DataFrame(info['cellboxes'])

    mp = Map(title=args.output)
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
        logging.debug('plotting fuel requirements')
        mp.Maps(mesh,'Fuel',predefined='Fuel',show=False)
    if 'speed' in mesh.columns:
        logging.debug('plotting vessel maximum speed')
        mp.Maps(mesh, 'speed', predefined = 'Speed', show = False)

    if 'paths' in info.keys():
        logging.debug('plotting paths')
        paths = info['paths']
        mp.Paths(paths,'Routes',predefined='Route Traveltime Paths')
    if 'waypoints' in info.keys():
        logging.debug('plotting waypoints')
        waypoints = pd.DataFrame(info['waypoints'])
        mp.Points(waypoints,'Waypoints',names={"font_size":10.0})

    mp.MeshInfo(mesh,'Mesh Info',show=False)

    logging.info('Saving plot to {}'.format(args.output))
    mp.save(args.output)




