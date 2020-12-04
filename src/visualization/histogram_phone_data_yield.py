import pandas as pd
import plotly.express as px


phone_data_yield = pd.read_csv(snakemake.input[0])

# make sure the input file contains "phone_data_yield_rapids_ratiovalidyieldedminutes" and "phone_data_yield_rapids_ratiovalidyieldedhours" columns
if ("phone_data_yield_rapids_ratiovalidyieldedminutes" not in phone_data_yield.columns) or ("phone_data_yield_rapids_ratiovalidyieldedhours" not in phone_data_yield.columns):
    raise ValueError("Please make sure [PHONE_DATA_YIELD][RAPIDS][COMPUTE] is True AND [PHONE_DATA_YIELD][RAPIDS][FEATURES] contains [ratiovalidyieldedminutes, ratiovalidyieldedhours].")

html_file = open(snakemake.output[0], "a", encoding="utf-8")
if phone_data_yield.empty:
    html_file.write("There is no sensor data for the sensors in [PHONE_DATA_YIELD][SENSORS].")
else:
    # plot ratio valid yielded minutes histogram
    fig_ratiovalidyieldedminutes = px.histogram(phone_data_yield, x="phone_data_yield_rapids_ratiovalidyieldedminutes", color="local_segment_label")
    fig_ratiovalidyieldedminutes.update_layout(title="Histogram of valid yielded minutes ratio per time segment.")
    html_file.write(fig_ratiovalidyieldedminutes.to_html(full_html=False, include_plotlyjs="cdn"))

    # plot ratio valid yielded hours histogram
    fig_ratiovalidyieldedhours = px.histogram(phone_data_yield, x="phone_data_yield_rapids_ratiovalidyieldedhours", color="local_segment_label")
    fig_ratiovalidyieldedhours.update_layout(title="Histogram of valid yielded hours ratio per time segment.")
    html_file.write(fig_ratiovalidyieldedhours.to_html(full_html=False, include_plotlyjs="cdn"))

html_file.close()
