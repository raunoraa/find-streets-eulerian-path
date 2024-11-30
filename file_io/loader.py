r"""
This submodule deals with loading the files.
"""

import os
import xml.etree.ElementTree as ET
import geojson


def load_osm_to_dict(file_path):
    """
    Load and parse an OSM XML file and convert it into a dictionary of ways and their associated tags.

    Args:
        file_path (str): The path to the OSM XML file to be parsed.

    Returns:
        dict: A dictionary where the keys are `way_id` and the values are dictionaries of tags for each way.
            The structure is {way_id: {tag1: value1, tag2: value2, ...}, ...}

    Example:

    ::

        osm_dict = load_osm_to_dict("path/to/osm_file.xml")
    """
    osm_dict = {}

    # Use `iterparse` to load the XML incrementally
    context = ET.iterparse(file_path, events=("start", "end"))
    for event, elem in context:
        if event == "end" and elem.tag == "way":
            osm_id = elem.get("id")
            if osm_id:
                # Collect tags into a dictionary for the current way
                tags = {tag.get("k"): tag.get("v") for tag in elem.findall("tag")}
                osm_dict[osm_id] = tags
            # Clear the element from memory to save space
            elem.clear()

    return osm_dict


def load_geojson(file_path):
    """
    Load and parse a GeoJSON file.

    Args:
        file_path (str): The path to the GeoJSON file to be loaded.

    Returns:
        dict: The parsed GeoJSON data as a dictionary.

    Example:

    ::

        geojson_data = load_geojson("path/to/geojson/file.geojson")
    """
    with open(file_path, "r") as f:
        return geojson.load(f)


def check_file_exists(file_path):
    """
    Check if a file exists at the specified path.

    Args:
        file_path (str): The path to the file to check.

    Returns:
        bool: True if the file exists, False otherwise.

    Example:

    ::

        if check_file_exists("path/to/file.txt"):
            print("File exists!")
        else:
            print("File not found!")
    """
    return os.path.isfile(file_path)
