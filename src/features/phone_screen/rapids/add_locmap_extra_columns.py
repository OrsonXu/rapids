import warnings
import numpy as np
import pandas as pd
import os


screen_data = pd.read_csv(snakemake.input["sensor_input_screen"])

if (not os.path.exists(snakemake.input["sensor_input_locmap"])):
    raise FileNotFoundError(f"Error: missing files for {locmap_type}. Please turn on Doryab and LocMap Location Providers")
location_data = pd.read_csv(snakemake.input["sensor_input_locmap"])

def find_intersections(screen_rages, locmap_ranges):
    idx_screen, idx_locmap = 0, 0
    judge_list = []
    while idx_screen < len(screen_ranges):
    # for idx_screen, (screen_start, screen_end) in enumerate(screen_ranges):
        screen_start, screen_end = screen_ranges[idx_screen]
        if (idx_locmap >= len(locmap_ranges)):
            judge_list.append(0)
            idx_screen += 1
            continue
        locmap_start, locmap_end = locmap_ranges[idx_locmap]
        if (screen_end <= locmap_start):
            judge_list.append(0)
            idx_screen += 1
        elif (screen_end <= locmap_end):
            judge_list.append(1)
            idx_screen += 1
        else: # screen_end > locmap_end
            if (screen_start < locmap_end):
                judge_list.append(1)
                idx_screen += 1
            else:
                idx_locmap += 1
    return judge_list

def find_in_locmap(row, location_data_locmap):
    start_timestamp = row["start_timestamp"]
    end_timestamp = row["end_timestamp"]
    return sum((location_data_locmap["start_timestamp"] < end_timestamp) & (location_data_locmap["end_timestamp"] > start_timestamp)) > 0

    if (location_data_locmap.empty):
        return False
    else:
        # A loose condition
        # An episode will be considered as long as there is any overlap
        
        # # Strategy 1 - This seems to be slow
        # query = location_data_locmap[(location_data_locmap["start_timestamp"] < end_timestamp) & (location_data_locmap["end_timestamp"] > start_timestamp)]
        # return not query.empty
        
        # Strategy 2 - Can be more efficient
        for idx, row_locmap in location_data_locmap.iterrows(): 
            if (start_timestamp < row_locmap["end_timestamp"]) \
                and (end_timestamp > row_locmap["start_timestamp"]): 
                return True
        return False

screen_ranges = screen_data[["start_timestamp", "end_timestamp"]].values
locmap_columns = [c for c in location_data.columns if c.startswith("locmap_isin_")]
for locmap_col in locmap_columns:
    locmap_ranges = location_data[location_data[locmap_col] == 1][["start_timestamp", "end_timestamp"]].values
    screen_data[locmap_col] = find_intersections(screen_ranges, locmap_ranges)

screen_data.to_csv(snakemake.output[0], index=False)