import os
import sys
import geopandas as gpd
from pyproj import CRS
from shapely.geometry import LineString, MultiLineString

from paths import *
from unit_conversions import *

'''
This script will take a shapefile as input and perform the following operations:
1. Reproject to new EPSG code
2. Simplify the geometry
3. Dissolve the geometry, export to a new shapefile
3. Buffer the geometry by a user-defined distance, export to a new shapefile
4. Buffer the geometry by 1000ft, export to new shapefile
'''
#using pyproj, this function will check that an input EPSG code is valid
try:
    in_shape = sys.argv[1]
    if not os.path.exists(in_shape):
        raise ValueError("Shapefile path does not exist.")
except IndexError:
    print("Please provide a shapefile path as the first argument.")
    sys.exit(1)

try:
    out_crs = int(sys.argv[2])
    out_crs = str(out_crs)
except ValueError:
    print("EPSG codes must only contain numbers.")
    sys.exit(1)

try:
    buffer_distance = int(sys.argv[3])
    if buffer_distance <= 0:
        raise ValueError("Buffer distance must be a positive integer.")
except ValueError:
    print("Buffer distance must be a positive integer.")
    sys.exit(1)

try:
    buffer_units = sys.argv[4]
    if buffer_units not in ['m', 'ft']:
        raise ValueError("Buffer units must be either 'm' or 'ft'.")
except IndexError:
    print("Please provide buffer units as the fourth argument.")
    sys.exit(1)

# Define the desired coordinate system and reproject
def reproject(out_crs, gdf):
    '''reproject the shapefile to a new coordinate system'''
    target_crs = CRS.from_epsg(out_crs)
    reprojected_gdf = gdf.to_crs(target_crs)
    return reprojected_gdf

# Dissolve the shapefile 
def dissolve(reprojected_gdf):
    '''dissolve the shapefile into a single geometry'''
    dissolved_geom = reprojected_gdf['geometry'].unary_union
    dissolved_gdf = gpd.GeoDataFrame(geometry=[dissolved_geom])
    return dissolved_gdf

def simplify_gdf(dissolved_gdf):
    '''simplify the geometry by a user-defined tolerance'''
    tolerance = 0.5 # Define the tolerance for simplification 
    simplified_geoms = []
    for geom in dissolved_gdf['geometry']:
        simplified_geom = geom.simplify(tolerance, preserve_topology=True)
        if isinstance(simplified_geom, MultiLineString):
            for line in simplified_geom:
                if isinstance(line, LineString):
                    simplified_geoms.append(line)
    # Convert the simplified geometries back to a GeoDataFrame
    simplified_gdf = gpd.GeoDataFrame(geometry=simplified_geoms)
    return simplified_gdf

def CRS_units(buffer_units, buffer_distance):
    '''parse the user-defined buffer units and calculate the correct buffer distance'''
    #read infile CRS, parse units
    crs = gpd.read_file(in_shape).crs
    unit_name = crs.axis_info[0].unit_name
    #options = metre, degree, US survey foot, foot
    if buffer_units == 'ft':
        if unit_name == "US survey foot" or unit_name == "foot":
            return buffer_distance
        elif unit_name == "metre":
            buffer_distance = us_survey_feet_to_meters(buffer_distance)
            return buffer_distance
        else:
            print("Something broke!")
            exit()
    elif buffer_units == 'm':
        if unit_name == "metre":
            return buffer_distance
        elif unit_name == "US survey foot" or unit_name == "foot":
            buffer_distance = meters_to_us_survey_feet(buffer_distance)  
            return buffer_distance
        else:
            print("Something broke!!!")
            exit()
    else:
        print("EPSG code provided is in degrees or is invalid. Please provide an EPSG code in US ft/ Intl Ft/ meters. Exiting...")
        exit()

def buffer_gdf(simplified_gdf, buffer_distance):
    '''buffer the simplified geometry by a user-defined distance'''
    buffered_geoms = [geom.buffer(buffer_distance/2) for geom in simplified_gdf['geometry']]
    buffered_gdf = gpd.GeoDataFrame(geometry=buffered_geoms) # Convert the buffered geometries back to a GeoDataFrame
    return buffered_gdf


def main(in_shape, out_crs, buffer_distance, buffer_units):

    # Load the shapefile
    gdf = gpd.read_file(in_shape)

    # Create Export File Paths
    parent_dir_with_filename = parse_filepath(in_shape)[1]
    ## "dissolved"
    dissolved_file_path = os.path.join(parent_dir_with_filename + "_simplified_dissolved.shp")
    ## "1000ft_buffer"  
    RR_buffer_path = os.path.join(parent_dir_with_filename + "_1000ft_buffer.shp")
    ## user-defined buffer distance + units 
    custom_buffer_path = os.path.join(parent_dir_with_filename + "_" + str(buffer_distance) + buffer_units + "_buffer.shp")

    # reproject, dissolve, simplify, and export
    reprojected_gdf = reproject(out_crs, gdf)
    dissolved_gdf = dissolve(reprojected_gdf)
    simplified_gdf = simplify_gdf(dissolved_gdf)
    simplified_gdf.to_file(dissolved_file_path, driver='ESRI Shapefile')

    # Custom Buffer
    buffer_distance = CRS_units(buffer_units, buffer_distance)
    buffered_gdf = buffer_gdf(simplified_gdf, buffer_distance)
    buffered_gdf = dissolve(buffered_gdf)

    # RR Buffer
    RR_buffer_distance = CRS_units("ft", 1000)
    RR_buffered_gdf = buffer_gdf(simplified_gdf, RR_buffer_distance)
    RR_buffered_gdf = dissolve(RR_buffered_gdf)

    # Apply CRS to the buffered geometries
    formatted_epsg = "EPSG:" + out_crs
    buffered_gdf.crs = formatted_epsg
    RR_buffered_gdf.crs = formatted_epsg

    # Export buffered shapes
    buffered_gdf.to_file(custom_buffer_path, driver='ESRI Shapefile')
    RR_buffered_gdf.to_file(RR_buffer_path, driver='ESRI Shapefile')

if __name__ == "__main__":
    main(in_shape, out_crs, buffer_distance, buffer_units)