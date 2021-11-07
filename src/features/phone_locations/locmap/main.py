from datetime import time
import numpy as np
import pandas as pd
import os

def find_range(index_list):
    adjacent_pairs = []
    l = len(index_list)
    if (l > 1):
        start = index_list[0]
        for i, idx in enumerate(index_list[:-1]):
            if (idx + 1 == index_list[i + 1]):
                continue
            else:
                end = idx
                adjacent_pairs.append((start,end))
                start = index_list[i + 1]
        end = index_list[-1]
        adjacent_pairs.append((start,end))
    elif(l==1):
        adjacent_pairs = [(index_list[0], index_list[0])]
    else:
        adjacent_pairs = []
    return adjacent_pairs

def duration_in_locmap_grpby_old(df_location_data, polygons, locmap_str):
    df_location_data["in_map"] = df_location_data.apply(
        lambda row: polygons.whether_contains_point(Point(row['double_longitude'], row['double_latitude'])),
        axis=1)
    index_in_locmap = np.where(df_location_data["in_map"])[0]
    index_pairs = find_range(index_in_locmap)
    duration_list = []
    timestmap_series = df_location_data["timestamp"]
    for (index_start, index_end) in index_pairs:
        if (index_start == 0):
            start_timestamp = timestmap_series.iloc[0]
        else:
            start_timestamp = (timestmap_series.iloc[index_start-1] + timestmap_series.iloc[index_start]) / 2 # get the middle timestamp
        if (index_end == (len(timestmap_series) - 1)):
            end_timestamp = timestmap_series.iloc[-1]
        else:
            end_timestamp = (timestmap_series.iloc[index_end] + timestmap_series.iloc[index_end + 1]) / 2 # get the middle timestamp
        duration_list.append(end_timestamp - start_timestamp)
    minutes_in_locmap = sum(duration_list) / 60000.0
    total_duration = (df_location_data["timestamp"].iloc[-1] - df_location_data["timestamp"].iloc[0]) / 60000.0
    percent_in_locmap = 100 * minutes_in_locmap / total_duration
    return pd.Series({f"duration_in_locmap_{locmap_str}": minutes_in_locmap, f"percent_in_locmap_{locmap_str}": percent_in_locmap})

def duration_in_locmap_grpby(df_location_data, locmap_str):
    minutes_in_locmap = df_location_data[df_location_data[f"locmap_isin_{locmap_str}"] == 1]["duration"].sum()
    total_duration = df_location_data["duration"].sum()
    percent_in_locmap = 100.0 * minutes_in_locmap / total_duration
    return pd.Series({f"duration_in_locmap_{locmap_str}": minutes_in_locmap, f"percent_in_locmap_{locmap_str}": percent_in_locmap})

def locmap_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    location_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_features = provider["FEATURES"]
    kml_folder_root = provider["KML_FILE_ROOT"]
    identifier = provider["IDENTIFIER"]
    locmap_types = provider["TYPES"]

    base_features = ["duration_in_locmap", "percent_in_locmap"]
    base_locmap_types = ["study","exercise","greens"]

    features_to_compute = list(set(requested_features) & set(base_features))
    locmaps_to_compute = list(set(locmap_types) & set(base_locmap_types))
    features_to_compute_final = [f"{feat}_{locmap}" for feat in features_to_compute for locmap in locmaps_to_compute]

    locmap_features = pd.DataFrame(columns=["local_segment"])

    # if not location_data.empty:
    #     location_data_filt = filter_data_by_segment(location_data, time_segment)
    #     for locmap_type in locmaps_to_compute:
    #         kml_filepath = os.path.join(kml_folder_root, identifier, locmap_type + ".kml")
    #         ps = polygon_set(kml_filepath)
    #         if os.path.exists(kml_filepath):
    #             location_data_duration_locmap = location_data_filt.groupby("local_segment").apply(
    #                 lambda x : duration_in_locmap_grpby(x, polygons = ps, locmap_str = locmap_type)).reset_index()
    #             if (locmap_features.shape[0] == 0):
    #                 locmap_features = location_data_duration_locmap.copy()
    #             else:
    #                 locmap_features = locmap_features.merge(location_data_duration_locmap,
    #                     left_on = "local_segment", right_on = "local_segment")
    
    if not location_data.empty:
        location_data_filt = filter_data_by_segment(location_data, time_segment)
        for locmap_type in locmaps_to_compute:
                location_data_duration_locmap = location_data_filt.groupby("local_segment").apply(
                    lambda x : duration_in_locmap_grpby(x, locmap_str = locmap_type)).reset_index()
                if (locmap_features.shape[0] == 0):
                    locmap_features = location_data_duration_locmap.copy()
                else:
                    locmap_features = locmap_features.merge(location_data_duration_locmap,
                        left_on = "local_segment", right_on = "local_segment")
    assert set(locmap_features.columns[1:]).issubset(set(features_to_compute_final))
    return locmap_features
