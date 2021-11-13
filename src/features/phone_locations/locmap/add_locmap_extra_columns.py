import warnings
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, Point
from fastkml import kml
import os

class polygon_set:
    def __init__(self, kml_filepath):
        self.polygons = []
        self.polygon_index = 0 # a pointer to save time for searching
        self.kml_filepath = kml_filepath
        self.empty_flag = False
        self.load_polygons()
    
    def get_current_polygon(self):
        if self.empty_flag:
            return None
        return self.polygons[self.polygon_index]

    def get_next_polygon(self):
        self.polygon_index += 1
        if (self.polygon_index >= self.num_polygons):
            self.polygon_index = 0
        return self.get_current_polygon()

    def load_polygons(self):
        listOfPolygons = []
        if not os.path.exists(self.kml_filepath):
            self.empty_flag = True
        else:
            with open(self.kml_filepath, "rb") as f:
                doc = f.read()
            k = kml.KML()
            k.from_string(doc)
            self.polygons = [g.geometry for f in k.features() for g in f.features() if g is not None]
        self.num_polygons = len(self.polygons)
    
    def whether_contains_point(self,point):
        # Leverage the fact of the continuity of lications,
        # so that a point could be in the same polygon if the previous point is
        if (self.get_current_polygon().contains(point)):
            return True
        count = 1
        while count < self.num_polygons:
            p = self.get_next_polygon()
            if (p.contains(point)):
                break
            count += 1
        return count < self.num_polygons

location_data = pd.read_csv(snakemake.input["sensor_input"])
provider = snakemake.params["provider"]
provider_doryab = snakemake.params["provider_doryab"]

kml_folder_root = provider["KML_FILE_ROOT"]
identifier = provider["IDENTIFIER"]
locmap_types = provider["TYPES"]
radius_from_home = provider_doryab["RADIUS_FOR_HOME"]

base_locmap_types = ["study","exercise","greens","living"]
locmaps_to_compute = list(set(locmap_types) & set(base_locmap_types))

for locmap_type in locmaps_to_compute:
    kml_filepath = os.path.join(kml_folder_root, identifier, locmap_type + ".kml")
    ps = polygon_set(kml_filepath)
    if os.path.exists(kml_filepath):
        if (location_data.empty):
            location_data[f"locmap_isin_{locmap_type}"] = np.nan
        else:
            location_data[f"locmap_isin_{locmap_type}"] = location_data.apply(
                    lambda row: int(ps.whether_contains_point(Point(row['double_longitude'], row['double_latitude']))),
                axis=1)
    else:
        raise FileNotFoundError(f"Error: missing files for {locmap_type}")

location_data[f"locmap_isin_home"] = location_data.apply(lambda row: 1 if row["distance_from_home"] <= radius_from_home else 0, axis=1)

location_data.to_csv(snakemake.output[0], index=False)