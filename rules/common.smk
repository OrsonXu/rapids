def get_script_language(script_path):
    from pathlib import Path
    script_path = Path(script_path)
    if not script_path.exists():
        raise ValueError("The following provider feature script does not exist: " + str(script_path))

    if script_path.name.endswith(".py"):
        return "python"
    elif script_path.name.endswith(".R"):
        return "r"


# Features.smk #########################################################################################################
def find_features_files(wildcards):
    feature_files = []
    for provider_key, provider in config[(wildcards.sensor_key).upper()]["PROVIDERS"].items():
        if provider["COMPUTE"]:
            feature_files.extend(expand("data/interim/{{pid}}/{sensor_key}_features/{sensor_key}_{language}_{provider_key}.csv", sensor_key=wildcards.sensor_key.lower(), language=get_script_language(provider["SRC_SCRIPT"]), provider_key=provider_key.lower()))
    return(feature_files)

def optional_steps_sleep_input(wildcards):
    if config["STEP"]["EXCLUDE_SLEEP"]["EXCLUDE"] == True and config["STEP"]["EXCLUDE_SLEEP"]["TYPE"] == "FITBIT_BASED":
        return  "data/raw/{pid}/fitbit_sleep_summary_with_datetime.csv"
    else:
        return []

def input_merge_sensor_features_for_individual_participants(wildcards):
    feature_files = []
    for config_key in config.keys():
        if config_key.startswith(("PHONE", "FITBIT", "EMPATICA")) and "PROVIDERS" in config[config_key] and isinstance(config[config_key]["PROVIDERS"], dict):
            for provider_key, provider in config[config_key]["PROVIDERS"].items():
                if "COMPUTE" in provider.keys() and provider["COMPUTE"]:
                    feature_files.append("data/processed/features/{pid}/" + config_key.lower() + ".csv")
                    break
    return feature_files

def get_phone_sensor_names():
    phone_sensor_names = []
    for config_key in config.keys():
        if config_key.startswith(("PHONE")) and "PROVIDERS" in config[config_key]:
            if config_key != "PHONE_DATA_YIELD" and config_key not in phone_sensor_names:
                    phone_sensor_names.append(config_key)
    return phone_sensor_names

def pull_phone_data_input_with_mutation_scripts(wilcards):
    import yaml
    input = dict()
    phone_stream = config["PHONE_DATA_STREAMS"]["USE"]

    input["participant_file"] = "data/external/participant_files/{pid}.yaml"
    input["rapids_schema_file"] = "src/data/streams/rapids_columns.yaml"
    input["stream_format"] = "src/data/streams/" + phone_stream + "/format.yaml"

    if Path("src/data/streams/"+ phone_stream + "/container.R").exists():
        input["stream_container"] = "src/data/streams/"+ phone_stream + "/container.R"
    elif Path("src/data/streams/"+ phone_stream + "/container.py").exists():
        input["stream_container"] = "src/data/streams/"+ phone_stream + "/container.py"
    else:
        raise ValueError("The container script for {stream} is missing: src/data/streams/{stream}/container.[py|R]".format(stream=empatica_stream))

    schema = yaml.load(open(input.get("stream_format"), 'r'), Loader=yaml.FullLoader)
    sensor = ("phone_" + wilcards.sensor).upper()
    if sensor not in schema:
        raise ValueError("{sensor} is not defined in the schema {schema}".format(sensor=sensor, schema=input.get("stream_format")))

    for device_os in schema[sensor].keys():
        if "MUTATION" not in schema[sensor][device_os]:
            raise ValueError("MUTATION is missing from [{sensor}][{device_os}] of {schema}".format(sensor=sensor, device_os=device_os,schema=input.get("stream_format")))
        if "COLUMN_MAPPINGS" not in schema[sensor][device_os]["MUTATION"]:
            raise ValueError("COLUMN_MAPPINGS is missing from [{sensor}][{device_os}][MUTATION] of {schema}".format(sensor=sensor, device_os=device_os, schema=input.get("stream_format")))
        if "SCRIPTS" not in schema[sensor][device_os]["MUTATION"]:
            raise ValueError("SCRIPTS is missing from [{sensor}][{device_os}][MUTATION] of {schema}".format(sensor=sensor, device_os=device_os, schema=input.get("stream_format")))

        scripts = schema[sensor][device_os]["MUTATION"]["SCRIPTS"]
        if isinstance(scripts, list):
            for idx, script in enumerate(scripts):
                if not script.lower().endswith((".py", ".r")):
                    raise ValueError("Mutate scripts can only be Python or R scripts (.py, .R).\n   Instead we got {script} in \n   [{sensor}][{device_os}] of {schema}".format(script=script, sensor=sensor, device_os=device_os, schema=input.get("stream_format")))
                input["mutationscript"+str(idx)] = script
    return input

def input_tzcodes_file(wilcards):
    from pathlib import Path
    if config["TIMEZONE"]["TYPE"] == "MULTIPLE":
        if not config["TIMEZONE"]["MULTIPLE"]["TZCODES_FILE"].lower().endswith(".csv"):
            raise ValueError("[TIMEZONE][MULTIPLE][TZCODES_FILE] should point to a CSV file, instead you typed: " + config["TIMEZONE"]["MULTIPLE"]["TZCODES_FILE"])
        if not Path(config["TIMEZONE"]["MULTIPLE"]["TZCODES_FILE"]).exists():
            raise ValueError("[TIMEZONE][MULTIPLE][TZCODES_FILE] should point to a CSV file, the file in the path you typed does not exist: " + config["TIMEZONE"]["MULTIPLE"]["TZCODES_FILE"])
        return [config["TIMEZONE"]["MULTIPLE"]["TZCODES_FILE"]]
    return []

def pull_wearable_data_input_with_mutation_scripts(wilcards):
    import yaml
    from pathlib import Path
    input = dict()
    device = wilcards.device_type.upper()
    device_stream = config[device+"_DATA_STREAMS"]["USE"]

    input["participant_file"] = "data/external/participant_files/{pid}.yaml"
    input["rapids_schema_file"] = "src/data/streams/rapids_columns.yaml"
    input["stream_format"] = "src/data/streams/" + device_stream + "/format.yaml"

    if Path("src/data/streams/"+ device_stream + "/container.R").exists():
        input["stream_container"] = "src/data/streams/"+ device_stream + "/container.R"
    elif Path("src/data/streams/"+ device_stream + "/container.py").exists():
        input["stream_container"] = "src/data/streams/"+ device_stream + "/container.py"
    else:
        raise ValueError("The container script for {stream} is missing: src/data/streams/{stream}/container.[py|R]".format(stream=device_stream))

    schema = yaml.load(open(input.get("stream_format"), 'r'), Loader=yaml.FullLoader)
    sensor = (device + "_" + wilcards.sensor).upper()
    if sensor not in schema:
        raise ValueError("{sensor} is not defined in the schema {schema}".format(sensor=sensor, schema=input.get("stream_format")))
    
    if "MUTATION" not in schema[sensor]:
        raise ValueError("MUTATION is missing from [{sensor}] of {schema}".format(sensor=sensor, schema=input.get("stream_format")))
    if "COLUMN_MAPPINGS" not in schema[sensor]["MUTATION"]:
        raise ValueError("COLUMN_MAPPINGS is missing from [{sensor}][MUTATION] of {schema}".format(sensor=sensor, schema=input.get("stream_format")))
    if "SCRIPTS" not in schema[sensor]["MUTATION"]:
        raise ValueError("SCRIPTS is missing from [{sensor}][MUTATION] of {schema}".format(sensor=sensor, schema=input.get("stream_format")))

    scripts = schema[sensor]["MUTATION"]["SCRIPTS"]
    if isinstance(scripts, list):
        for idx, script in enumerate(scripts):
            if not script.lower().endswith((".py", ".r")):
                raise ValueError("Mutate scripts can only be Python or R scripts (.py, .R).\n   Instead we got {script} in [{sensor}] of {schema}".format(script=script, sensor=sensor, schema=input.get("stream_format")))
            input["mutationscript"+str(idx)] = script
    return input

