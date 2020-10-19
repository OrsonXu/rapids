rule restore_sql_file:
    input:
        sql_file = "data/external/rapids_example.sql",
        db_credentials = ".env"
    params:
        group = config["DOWNLOAD_PARTICIPANTS"]["GROUP"]
    output:
        touch("data/interim/restore_sql_file.done")
    script:
        "../src/data/restore_sql_file.py"

rule create_example_participant_files:
    output:
        expand("data/external/{pid}", pid = ["example01", "example02"])
    shell:
        "echo 'a748ee1a-1d0b-4ae9-9074-279a2b6ba524\nandroid\ntest01\n2020/04/23,2020/05/04\n' >> ./data/external/example01 && echo '13dbc8a3-dae3-4834-823a-4bc96a7d459d\nios\ntest02\n2020/04/23,2020/05/04\n' >> ./data/external/example02"

rule download_participants:
    params:
        group = config["DOWNLOAD_PARTICIPANTS"]["GROUP"],
        ignored_device_ids = config["DOWNLOAD_PARTICIPANTS"]["IGNORED_DEVICE_IDS"],
        timezone = config["TIMEZONE"]
    priority: 1
    script:
        "../src/data/download_participants.R"

rule download_dataset:
    input:
        "data/external/{pid}"
    params:
        group = config["DOWNLOAD_DATASET"]["GROUP"],
        sensor = "{sensor}",
        table = lambda wildcards: config[str(wildcards.sensor).upper()]["DB_TABLE"],
        timezone = config["TIMEZONE"],
        aware_multiplatform_tables = config["PHONE_ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"] + "," + config["PHONE_ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"] + "," + config["PHONE_CONVERSATION"]["DB_TABLE"]["ANDROID"] + "," + config["PHONE_CONVERSATION"]["DB_TABLE"]["IOS"],
    output:
        "data/raw/{pid}/{sensor}_raw.csv"
    script:
        "../src/data/download_dataset.R"

rule compute_day_segments:
    input: 
        config["DAY_SEGMENTS"]["FILE"]
    params:
        day_segments_type = config["DAY_SEGMENTS"]["TYPE"],
        pid = "{pid}"
    output:
        segments_file = "data/interim/day_segments/{pid}_day_segments.csv",
        segments_labels_file = "data/interim/day_segments/{pid}_day_segments_labels.csv",
    script:
        "../src/data/compute_day_segments.py"

rule phone_readable_datetime:
    input:
        sensor_input = "data/raw/{pid}/phone_{sensor}_raw.csv",
        day_segments = "data/interim/day_segments/{pid}_day_segments.csv"
    params:
        timezones = None,
        fixed_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        day_segments_type = config["DAY_SEGMENTS"]["TYPE"],
        include_past_periodic_segments = config["DAY_SEGMENTS"]["INCLUDE_PAST_PERIODIC_SEGMENTS"]
    output:
        "data/raw/{pid}/phone_{sensor}_with_datetime.csv"
    script:
        "../src/data/readable_datetime.R"

rule phone_sensed_bins:
    input:
        all_sensors = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor = map(str.lower, config["PHONE_VALID_SENSED_BINS"]["PHONE_SENSORS"]))
    params:
        bin_size = config["PHONE_VALID_SENSED_BINS"]["BIN_SIZE"]
    output:
        "data/interim/{pid}/phone_sensed_bins.csv"
    script:
        "../src/data/phone_sensed_bins.R"

rule phone_sensed_timestamps:
    input:
        all_sensors = expand("data/raw/{{pid}}/{sensor}_raw.csv", sensor = map(str.lower, config["PHONE_VALID_SENSED_BINS"]["PHONE_SENSORS"]))
    output:
        "data/interim/{pid}/phone_sensed_timestamps.csv"
    script:
        "../src/data/phone_sensed_timestamps.R"

rule phone_valid_sensed_days:
    input:
        phone_sensed_bins =  "data/interim/{pid}/phone_sensed_bins.csv"
    params:
        min_valid_hours_per_day = "{min_valid_hours_per_day}",
        min_valid_bins_per_hour = "{min_valid_bins_per_hour}"
    output:
        "data/interim/{pid}/phone_valid_sensed_days_{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins.csv"
    script:
        "../src/data/phone_valid_sensed_days.R"


rule unify_ios_android:
    input:
        sensor_data = "data/raw/{pid}/{sensor}_with_datetime.csv",
        participant_info = "data/external/{pid}"
    params:
        sensor = "{sensor}",
    output:
        "data/raw/{pid}/{sensor}_with_datetime_unified.csv"
    script:
        "../src/data/unify_ios_android.R"

rule process_phone_location_types:
    input:
        locations = "data/raw/{pid}/phone_locations_raw.csv",
        phone_sensed_timestamps = "data/interim/{pid}/phone_sensed_timestamps.csv",
    params:
        consecutive_threshold = config["PHONE_LOCATIONS"]["FUSED_RESAMPLED_CONSECUTIVE_THRESHOLD"],
        time_since_valid_location = config["PHONE_LOCATIONS"]["FUSED_RESAMPLED_TIME_SINCE_VALID_LOCATION"],
        locations_to_use = config["PHONE_LOCATIONS"]["LOCATIONS_TO_USE"]
    output:
        "data/interim/{pid}/phone_locations_processed.csv"
    script:
        "../src/data/process_location_types.R"

rule readable_datetime_location_processed:
    input:
        sensor_input = "data/interim/{pid}/phone_locations_processed.csv",
        day_segments = "data/interim/day_segments/{pid}_day_segments.csv"
    params:
        timezones = None,
        fixed_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        day_segments_type = config["DAY_SEGMENTS"]["TYPE"],
        include_past_periodic_segments = config["DAY_SEGMENTS"]["INCLUDE_PAST_PERIODIC_SEGMENTS"]
    output:
        "data/interim/{pid}/phone_locations_processed_with_datetime.csv"
    script:
        "../src/data/readable_datetime.R"

rule phone_application_categories:
    input:
        "data/raw/{pid}/phone_applications_foreground_with_datetime.csv"
    params:
        catalogue_source = config["PHONE_APPLICATIONS_FOREGROUND"]["APPLICATION_CATEGORIES"]["CATALOGUE_SOURCE"],
        catalogue_file = config["PHONE_APPLICATIONS_FOREGROUND"]["APPLICATION_CATEGORIES"]["CATALOGUE_FILE"],
        update_catalogue_file = config["PHONE_APPLICATIONS_FOREGROUND"]["APPLICATION_CATEGORIES"]["UPDATE_CATALOGUE_FILE"],
        scrape_missing_genres = config["PHONE_APPLICATIONS_FOREGROUND"]["APPLICATION_CATEGORIES"]["SCRAPE_MISSING_CATEGORIES"]
    output:
        "data/raw/{pid}/phone_applications_foreground_with_datetime_with_categories.csv"
    script:
        "../src/data/application_categories.R"

rule fitbit_heartrate_with_datetime:
    input:
        expand("data/raw/{{pid}}/{fitbit_table}_raw.csv", fitbit_table=config["HEARTRATE"]["DB_TABLE"])
    params:
        local_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        fitbit_sensor = "heartrate"
    output:
        summary_data = "data/raw/{pid}/fitbit_heartrate_summary_with_datetime.csv",
        intraday_data = "data/raw/{pid}/fitbit_heartrate_intraday_with_datetime.csv"
    script:
        "../src/data/fitbit_readable_datetime.py"

rule fitbit_step_with_datetime:
    input:
        expand("data/raw/{{pid}}/{fitbit_table}_raw.csv", fitbit_table=config["STEP"]["DB_TABLE"])
    params:
        local_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        fitbit_sensor = "steps"
    output:
        intraday_data = "data/raw/{pid}/fitbit_step_intraday_with_datetime.csv"
    script:
        "../src/data/fitbit_readable_datetime.py"

rule fitbit_sleep_with_datetime:
    input:
        expand("data/raw/{{pid}}/{fitbit_table}_raw.csv", fitbit_table=config["SLEEP"]["DB_TABLE"])
    params:
        local_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        fitbit_sensor = "sleep"
    output:
        summary_data = "data/raw/{pid}/fitbit_sleep_summary_with_datetime.csv",
        intraday_data = "data/raw/{pid}/fitbit_sleep_intraday_with_datetime.csv"
    script:
        "../src/data/fitbit_readable_datetime.py"
