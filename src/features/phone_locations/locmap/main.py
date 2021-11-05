from datetime import time
import numpy as np
import pandas as pd


def locmap_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    location_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_features = provider["FEATURES"]
    print(sensor_data_files)
    print(requested_features)
