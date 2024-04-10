import logging
import time
import tracemalloc
import numpy as np
import geopandas as gpd
import json
from functools import wraps
from shapely.geometry import Polygon, MultiPolygon, Point
import shapely

def memory_trace(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        tracemalloc.start(20)
        res = func(*args, **kwargs)
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('traceback')

        stat = top_stats[0]
        logging.info("{} memory blocks: {.1f} KiB".
                     format(stat.count, stat.size / 1024))
        logging.info("\n".join(stat.traceback.format()))
        return res
    return wrapper


def timed_call(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        res = func(*args, **kwargs)
        end = time.perf_counter()
        logging.info("Timed call to {} took {:02f} seconds".
                     format(func.__name__, end - start))
        return res
    return wrapper


def setup_logging(func,
                  log_format="[%(asctime)-17s :%(levelname)-8s] - %(message)s"):
    """Wraps a CLI endpoint and sets up logging for it

    This is probably not the smoothest implementation, but it's an educational
    one for people who aren't aware of decorators and how they're implemented.
    In addition, it supports a nice pattern for CLI endpoints

    TODO: start handling level configuration from logging yaml config

    :param func:
    :param log_format:
    :return:
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        parsed_args = func(*args, **kwargs)
        level = logging.INFO

        if hasattr(parsed_args, "verbose") and parsed_args.verbose:
            level = logging.DEBUG

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt="%d-%m-%y %T",
        )

        logging.getLogger("cdsapi").setLevel(logging.WARNING)
        logging.getLogger("matplotlib").setLevel(logging.WARNING)
        logging.getLogger("matplotlib.pyplot").setLevel(logging.WARNING)
        logging.getLogger("requests").setLevel(logging.WARNING)
        logging.getLogger("tensorflow").setLevel(logging.WARNING)
        logging.getLogger("urllib3").setLevel(logging.WARNING)
        return parsed_args
    return wrapper


def convert_decimal_days(decimal_days, mins=False):
    """
    Convert decimal days to more readable Days, Hours and (optionally) Minutes
    Args:
        decimal_days (float): Number of days as a decimal
        mins (bool): Determines whether to return minutes or decimal hours
    Returns:
        new_time (str): The time in the new format
    """
    frac_d, days = np.modf(decimal_days)
    hours = frac_d * 24.0

    if mins:
        frac_h, hours = np.modf(hours)
        minutes = round(frac_h * 60.0)
        if days:
            if round(days) == 1:
                new_time = f"{round(days)} day {round(hours)} hours {minutes} minutes"
            else:
                new_time = f"{round(days)} days {round(hours)} hours {minutes} minutes"
        elif hours:
            new_time = f"{round(hours)} hours {minutes} minutes"
        else:
            new_time = f"{minutes} minutes"
    else:
        hours = round(hours, 2)
        if days:
            if round(days) == 1:
                new_time = f"{round(days)} day {hours} hours"
            else:
                new_time = f"{round(days)} days {hours} hours"
        else:
            new_time = f"{hours} hours"

    return new_time

def gpx_route_import(f_name):
    """
        Function to import a route in gpx format and convert it to geojson format

        Args:
            f_name: Filename of gpx route file

        Returns:
            geojson: Route in geojson format
    """
    gdf_r = gpd.read_file(f_name, layer="routes")
    gdf_p = gpd.read_file(f_name, layer="route_points")

    # Drop empty fields from original gpx file
    gdf_r = gdf_r.dropna(how='all', axis=1)
    # Convert route to geojson linestring
    geojson = json.loads(gdf_r.to_json())

    # Extract start and end waypoints and add to geojson properties
    geojson['features'][0]['properties']['from'] = gdf_p['name'].iloc[0]
    geojson['features'][0]['properties']['to'] = gdf_p['name'].iloc[-1]
    # Spoof traveltime so that it plots successfully
    geojson['features'][0]['properties']['traveltime'] = [-1]*len(
                            geojson['features'][0]['geometry']['coordinates'])

    return geojson

def split_at_antimeridian(gdf):
    """
    Splits a cellbox at the antimeridian into two seperate polygons, and 
    combines them into a MultiPolygon. Also adds a wavy edge to the antimeridian
    side in order to denote a continuous cellbox on each side of the antimeridian
    """
    def create_ellipse(centre=(0,0), semi_x_y=(1,1)):
        """
        Creates a shapely object of an ellipse centred at 'centre', 
        and with semi major/minor axes of 'semi_x_y'
        """
        circ = Point(centre).buffer(0.5)
        ellipse = shapely.affinity.scale(circ, semi_x_y[0], semi_x_y[1])
        return ellipse
    
    # Split cellboxes that cross the antimeridian into multipolygons on each side
    for idx, row in gdf.iterrows():
        shape = row['geometry']
        # If it crosses antimeridian, it will be a multipolygon
        if shape.geom_type == 'MultiPolygon':
            polygons = list(shape.geoms)
            

            poly_a = polygons[0]
            poly_b = polygons[1]
            
            # Create wavy line to represent a cellbox that goes over antimeridian
            xx, yy = poly_a.exterior.coords.xy
            cb_height = yy[1] - yy[0]
            cb_width = 180 - xx[0]
            # Create two scaled ellipses along the cellbox boundary on the antimeridian
            ellipse_a_1 = create_ellipse(centre = (180, yy[0] + 0.75*cb_height), 
                                         semi_x_y=(cb_width*0.25, cb_height*0.5))
            ellipse_a_2 = create_ellipse(centre = (180, yy[0] + 0.25*cb_height), 
                                         semi_x_y=(cb_width*0.25, cb_height*0.5))
            # Do same on negative boundary
            xx, yy = poly_b.exterior.coords.xy
            cb_width = 180 + xx[2]
            ellipse_b_1 = create_ellipse(centre = (-180, yy[0] + 0.75*cb_height), 
                                         semi_x_y=(cb_width*0.25, cb_height*0.5))
            ellipse_b_2 = create_ellipse(centre = (-180, yy[0] + 0.25*cb_height), 
                                         semi_x_y=(cb_width*0.25, cb_height*0.5))
            # Add/subtract ellipses from antimeridian line
            poly_a = poly_a.union(ellipse_a_1).difference(ellipse_a_2)
            poly_b = poly_b.union(ellipse_b_2).difference(ellipse_b_1)

            # Combine into a multipolygon and replace original
            new_polygon = MultiPolygon([poly_a, poly_b])
            # Replace the cellbox in place in the original dataframe
            gdf.loc[idx, 'geometry'] = new_polygon
    
    return gdf