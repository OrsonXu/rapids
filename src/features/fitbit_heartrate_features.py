import pandas as pd
import numpy as np
from scipy.stats import entropy
import json


heartrate_data = pd.read_csv(snakemake.input[0], parse_dates=["local_date_time", "local_date"])
day_segment = snakemake.params["day_segment"]
features = snakemake.params["features"]


heartrate_features = pd.DataFrame(columns=["local_date"] + ["heartrate_" + day_segment + "_" + x for x in features])
if not heartrate_data.empty:
    device_id = heartrate_data["device_id"][0]
    num_rows_per_minute = heartrate_data.groupby(["local_date", "local_hour", "local_minute"]).count().mean()["device_id"]
    if day_segment != "daily":
        heartrate_data =heartrate_data[heartrate_data["local_day_segment"] == day_segment]
    
    if not heartrate_data.empty:
        heartrate_features = pd.DataFrame()
     
        # get stats of heartrate
        if "maxhr" in features:
            heartrate_features["heartrate_" + day_segment + "_maxhr"] = heartrate_data.groupby(["local_date"])["heartrate"].max()
        if "minhr" in features:
            heartrate_features["heartrate_" + day_segment + "_minhr"] = heartrate_data.groupby(["local_date"])["heartrate"].min()
        if "avghr" in features:
            heartrate_features["heartrate_" + day_segment + "_avghr"] = heartrate_data.groupby(["local_date"])["heartrate"].mean()
        if "medianhr" in features:
            heartrate_features["heartrate_" + day_segment + "_medianhr"] = heartrate_data.groupby(["local_date"])["heartrate"].median()
        if "modehr" in features:
            heartrate_features["heartrate_" + day_segment + "_modehr"] = heartrate_data.groupby(["local_date"])["heartrate"].agg(lambda x: pd.Series.mode(x)[0])
        if "stdhr" in features:
            heartrate_features["heartrate_" + day_segment + "_stdhr"] = heartrate_data.groupby(["local_date"])["heartrate"].std()
        if "diffmaxmodehr" in features:
            heartrate_features["heartrate_" + day_segment + "_diffmaxmodehr"] = heartrate_data.groupby(["local_date"])["heartrate"].max() - heartrate_data.groupby(["local_date"])["heartrate"].agg(lambda x: pd.Series.mode(x)[0])
        if "diffminmodehr" in features:
            heartrate_features["heartrate_" + day_segment + "_diffminmodehr"] = heartrate_data.groupby(["local_date"])["heartrate"].agg(lambda x: pd.Series.mode(x)[0]) - heartrate_data.groupby(["local_date"])["heartrate"].min()
        if "entropyhr" in features:
            heartrate_features["heartrate_" + day_segment + "_entropyhr"] = heartrate_data.groupby(["local_date"])["heartrate"].agg(entropy)

        # get number of minutes in each heart rate zone
        for feature_name in list(set(["lengthoutofrange", "lengthfatburn", "lengthcardio", "lengthpeak"]) & set(features)):
            heartrate_zone = heartrate_data[heartrate_data["heartrate_zone"] == feature_name[6:]]
            heartrate_features["heartrate_" + day_segment + "_" + feature_name] = heartrate_zone.groupby(["local_date"])["device_id"].count() / num_rows_per_minute
            heartrate_features.fillna(value={"heartrate_" + day_segment + "_" + feature_name: 0}, inplace=True)

        heartrate_features = heartrate_features.reset_index()

heartrate_features.to_csv(snakemake.output[0], index=False)
