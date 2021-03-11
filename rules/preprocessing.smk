# This rull will dissapear
rule restore_sql_file:
    input:
        sql_file = "data/external/rapids_example.sql",
        db_credentials = ".env"
    params:
        group = "No_GROUP" #config["DATABASE_GROUP"]
    output:
        touch("data/interim/restore_sql_file.done")
    script:
        "../src/data/restore_sql_file.py"

rule create_example_participant_files:
    output:
        expand("data/external/participant_files/{pid}.yaml", pid = ["example01", "example02"])
    shell:
        "echo 'PHONE:\n  DEVICE_IDS: [a748ee1a-1d0b-4ae9-9074-279a2b6ba524]\n  PLATFORMS: [android]\n  LABEL: test-01\n  START_DATE: 2020-04-23\n  END_DATE: 2020-05-04\nFITBIT:\n  DEVICE_IDS: [a748ee1a-1d0b-4ae9-9074-279a2b6ba524]\n  LABEL: test-01\n  START_DATE: 2020-04-23\n  END_DATE: 2020-05-04\n' >> ./data/external/participant_files/example01.yaml && echo 'PHONE:\n  DEVICE_IDS: [13dbc8a3-dae3-4834-823a-4bc96a7d459d]\n  PLATFORMS: [ios]\n  LABEL: test-02\n  START_DATE: 2020-04-23\n  END_DATE: 2020-05-04\nFITBIT:\n  DEVICE_IDS: [13dbc8a3-dae3-4834-823a-4bc96a7d459d]\n  LABEL: test-02\n  START_DATE: 2020-04-23\n  END_DATE: 2020-05-04\n' >> ./data/external/participant_files/example02.yaml"

rule create_participants_files:
    input:
        participants_file = config["CREATE_PARTICIPANT_FILES"]["CSV_FILE_PATH"] 
    params:
        config = config["CREATE_PARTICIPANT_FILES"]
    script:
        "../src/data/create_participants_files.R"

rule pull_phone_data:
    input: unpack(pull_phone_data_input_with_mutation_scripts)
    params:
        data_configuration = config["PHONE_DATA_STREAMS"][config["PHONE_DATA_STREAMS"]["USE"]],
        sensor = "phone_" + "{sensor}",
        tables = lambda wildcards: config["PHONE_" + str(wildcards.sensor).upper()]["TABLE"],
    output:
        "data/raw/{pid}/phone_{sensor}_raw.csv"
    script:
        "../src/data/streams/pull_phone_data.R"

rule compute_time_segments:
    input: 
        config["TIME_SEGMENTS"]["FILE"],
        "data/external/participant_files/{pid}.yaml"
    params:
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"],
        pid = "{pid}"
    output:
        segments_file = "data/interim/time_segments/{pid}_time_segments.csv",
        segments_labels_file = "data/interim/time_segments/{pid}_time_segments_labels.csv",
    script:
        "../src/data/compute_time_segments.py"

rule phone_readable_datetime:
    input:
        sensor_input = "data/raw/{pid}/phone_{sensor}_raw.csv",
        time_segments = "data/interim/time_segments/{pid}_time_segments.csv",
        pid_file = "data/external/participant_files/{pid}.yaml",
        tzcodes_file = input_tzcodes_file,
    params:
        device_type = "phone",
        timezone_parameters = config["TIMEZONE"],
        pid = "{pid}",
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"],
        include_past_periodic_segments = config["TIME_SEGMENTS"]["INCLUDE_PAST_PERIODIC_SEGMENTS"]
    output:
        "data/raw/{pid}/phone_{sensor}_with_datetime.csv"
    script:
        "../src/data/datetime/readable_datetime.R"

rule phone_yielded_timestamps:
    input:
        all_sensors = expand("data/raw/{{pid}}/{sensor}_raw.csv", sensor = map(str.lower, config["PHONE_DATA_YIELD"]["SENSORS"]))
    params:
        sensors = config["PHONE_DATA_YIELD"]["SENSORS"] # not used but needed so the rule is triggered if this array changes
    output:
        "data/interim/{pid}/phone_yielded_timestamps.csv"
    script:
        "../src/data/phone_yielded_timestamps.R"

rule phone_yielded_timestamps_with_datetime:
    input:
        sensor_input = "data/interim/{pid}/phone_yielded_timestamps.csv",
        time_segments = "data/interim/time_segments/{pid}_time_segments.csv",
        pid_file = "data/external/participant_files/{pid}.yaml",
        tzcodes_file = input_tzcodes_file,
    params:
        device_type = "phone",
        timezone_parameters = config["TIMEZONE"],
        pid = "{pid}",
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"],
        include_past_periodic_segments = config["TIME_SEGMENTS"]["INCLUDE_PAST_PERIODIC_SEGMENTS"]
    output:
        "data/interim/{pid}/phone_yielded_timestamps_with_datetime.csv"
    script:
        "../src/data/datetime/readable_datetime.R"

rule unify_ios_android:
    input:
        sensor_data = "data/raw/{pid}/{sensor}_with_datetime.csv",
        participant_info = "data/external/participant_files/{pid}.yaml"
    params:
        sensor = "{sensor}",
    output:
        "data/raw/{pid}/{sensor}_with_datetime_unified.csv"
    script:
        "../src/data/unify_ios_android.R"

rule process_phone_locations_types:
    input:
        locations = "data/raw/{pid}/phone_locations_raw.csv",
        phone_sensed_timestamps = "data/interim/{pid}/phone_yielded_timestamps.csv",
    params:
        consecutive_threshold = config["PHONE_LOCATIONS"]["FUSED_RESAMPLED_CONSECUTIVE_THRESHOLD"],
        time_since_valid_location = config["PHONE_LOCATIONS"]["FUSED_RESAMPLED_TIME_SINCE_VALID_LOCATION"],
        locations_to_use = config["PHONE_LOCATIONS"]["LOCATIONS_TO_USE"]
    output:
        "data/interim/{pid}/phone_locations_processed.csv"
    script:
        "../src/data/process_location_types.R"

rule phone_locations_processed_with_datetime:
    input:
        sensor_input = "data/interim/{pid}/phone_locations_processed.csv",
        time_segments = "data/interim/time_segments/{pid}_time_segments.csv",
        pid_file = "data/external/participant_files/{pid}.yaml",
        tzcodes_file = input_tzcodes_file,
    params:
        device_type = "phone",
        timezone_parameters = config["TIMEZONE"],
        pid = "{pid}",
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"],
        include_past_periodic_segments = config["TIME_SEGMENTS"]["INCLUDE_PAST_PERIODIC_SEGMENTS"]
    output:
        "data/interim/{pid}/phone_locations_processed_with_datetime.csv"
    script:
        "../src/data/datetime/readable_datetime.R"

rule phone_locations_processed_with_datetime_with_home:
    input:
        sensor_input = "data/interim/{pid}/phone_locations_processed_with_datetime.csv"
    params:
        dbscan_eps = config["PHONE_LOCATIONS"]["HOME_INFERENCE"]["DBSCAN_EPS"],
        dbscan_minsamples = config["PHONE_LOCATIONS"]["HOME_INFERENCE"]["DBSCAN_MINSAMPLES"],
        threshold_static = config["PHONE_LOCATIONS"]["HOME_INFERENCE"]["THRESHOLD_STATIC"],
        clustering_algorithm = config["PHONE_LOCATIONS"]["HOME_INFERENCE"]["CLUSTERING_ALGORITHM"]
    output: 
        "data/interim/{pid}/phone_locations_processed_with_datetime_with_home.csv"
    script:
        "../src/data/infer_home_location.py"

rule resample_episodes:
    input:
        "data/interim/{pid}/{sensor}_episodes.csv"
    output:
        "data/interim/{pid}/{sensor}_episodes_resampled.csv"
    script:
        "../src/features/utils/resample_episodes.R"

rule resample_episodes_with_datetime:
    input:
        sensor_input = "data/interim/{pid}/{sensor}_episodes_resampled.csv",
        time_segments = "data/interim/time_segments/{pid}_time_segments.csv",
        pid_file = "data/external/participant_files/{pid}.yaml",
        tzcodes_file = input_tzcodes_file,
    params:
        device_type = lambda wildcards: wildcards.sensor.split("_")[0],
        timezone_parameters = config["TIMEZONE"],
        pid = "{pid}",
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"],
        include_past_periodic_segments = config["TIME_SEGMENTS"]["INCLUDE_PAST_PERIODIC_SEGMENTS"]
    output:
        "data/interim/{pid}/{sensor}_episodes_resampled_with_datetime.csv"
    script:
        "../src/data/datetime/readable_datetime.R"


rule phone_application_categories:
    input:
        "data/raw/{pid}/phone_applications_{type}_with_datetime.csv"
    params:
        catalogue_source = lambda wildcards: config["PHONE_APPLICATIONS_" + str(wildcards.type).upper()]["APPLICATION_CATEGORIES"]["CATALOGUE_SOURCE"],
        catalogue_file = lambda wildcards: config["PHONE_APPLICATIONS_" + str(wildcards.type).upper()]["APPLICATION_CATEGORIES"]["CATALOGUE_FILE"],
        update_catalogue_file = lambda wildcards: config["PHONE_APPLICATIONS_" + str(wildcards.type).upper()]["APPLICATION_CATEGORIES"]["UPDATE_CATALOGUE_FILE"],
        scrape_missing_genres = lambda wildcards: config["PHONE_APPLICATIONS_" + str(wildcards.type).upper()]["APPLICATION_CATEGORIES"]["SCRAPE_MISSING_CATEGORIES"]
    output:
        "data/raw/{pid}/phone_applications_{type}_with_datetime_with_categories.csv"
    script:
        "../src/data/application_categories.R"

rule pull_wearable_data:
    input: unpack(pull_wearable_data_input_with_mutation_scripts)
    params:
        data_configuration = lambda wildcards: config[wildcards.device_type.upper() +"_DATA_STREAMS"][config[wildcards.device_type.upper() +"_DATA_STREAMS"]["USE"]],
        device_type = "{device_type}",
        sensor = "{device_type}" + "_" + "{sensor}",
        pid = "{pid}",
        tables = lambda wildcards: config[wildcards.device_type.upper() + "_" + str(wildcards.sensor).upper()]["TABLE"],
    wildcard_constraints:
        device_type="(empatica|fitbit)"
    output:
        "data/raw/{pid}/{device_type}_{sensor}_raw.csv"
    script:
        "../src/data/streams/pull_wearable_data.R"

rule fitbit_readable_datetime:
    input:
        sensor_input = "data/raw/{pid}/fitbit_{sensor}_raw.csv",
        time_segments = "data/interim/time_segments/{pid}_time_segments.csv",
        pid_file = "data/external/participant_files/{pid}.yaml",
        tzcodes_file = input_tzcodes_file,
    params:
        device_type = "fitbit",
        timezone_parameters = config["TIMEZONE"],
        pid = "{pid}",
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"],
        include_past_periodic_segments = config["TIME_SEGMENTS"]["INCLUDE_PAST_PERIODIC_SEGMENTS"]
    output:
        "data/raw/{pid}/fitbit_{sensor}_with_datetime.csv"
    script:
        "../src/data/datetime/readable_datetime.R"

rule empatica_readable_datetime:
    input:
        sensor_input = "data/raw/{pid}/empatica_{sensor}_raw.csv",
        time_segments = "data/interim/time_segments/{pid}_time_segments.csv",
        pid_file = "data/external/participant_files/{pid}.yaml",
        tzcodes_file = input_tzcodes_file,
    params:
        device_type = "empatica",
        timezone_parameters = config["TIMEZONE"],
        pid = "{pid}",
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"],
        include_past_periodic_segments = config["TIME_SEGMENTS"]["INCLUDE_PAST_PERIODIC_SEGMENTS"]
    output:
        "data/raw/{pid}/empatica_{sensor}_with_datetime.csv"
    script:
        "../src/data/datetime/readable_datetime.R"
