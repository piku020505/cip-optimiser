"""CIP optimiser explorer:  streamlit run app.py"""
import pandas as pd
import streamlit as st
from cipopt import Plant, Recipe, simulate, wash_time_needed, rinse_volume_to_endpoint

st.set_page_config(page_title="CIP optimiser", layout="wide")
st.title("CIP optimiser · clean to spec with less water, caustic and energy")
st.caption("Pick a caustic temperature and NaOH strength. The app sets the wash time needed for a 3-log soil reduction (+20% margin) and stops the rinse on conductivity, then compares with today's fixed-time recipe.")
with st.sidebar:
    temp = st.slider("Caustic temperature (°C)", 60.0, 85.0, 74.0, 0.5)
    conc = st.slider("NaOH concentration (%)", 0.75, 2.5, 1.15, 0.05)
    rec = st.checkbox("Caustic recovery tank", True)
    st.divider(); st.caption("Plant assumptions")
    flow = st.number_input("Flow (L/min)", 60.0, 400.0, 180.0, 10.0)
    hold = st.number_input("Circuit hold-up (L)", 50.0, 600.0, 150.0, 10.0)
p = Plant(flow_lpm=flow, holdup_l=hold)
base = simulate(p, Recipe())
r = Recipe(temp_c=temp, conc_pct=conc, wash_min=wash_time_needed(p, temp, conc), rinse1_l=rinse_volume_to_endpoint(p, conc), recovery=rec)
new = simulate(p, r)
cols = st.columns(5)
for c, (k, lab, unit) in zip(cols, [("water_l", "Water", "L"), ("naoh_kg", "NaOH", "kg"), ("fuel_kwh", "Boiler fuel", "kWh"), ("co2_kg", "CO₂", "kg"), ("duration_min", "Duration", "min")]):
    c.metric(lab, f"{new[k]:.1f} {unit}", f"{new[k] - base[k]:+.1f} vs today", delta_color="inverse")
st.write(f"Wash time set to **{r.wash_min:.1f} min** · soil reduction **{new['log_reduction']:.1f}-log** (target {p.log_target:.0f}) · intermediate rinse stops at **{r.rinse1_l:.0f} L** (residual NaOH {new['residual_naoh_ppm']:.1f} ppm) · cost **₹{new['cost_inr']:.0f}** vs ₹{base['cost_inr']:.0f} today")
if new["duration_min"] > base["duration_min"]:
    st.warning("This recipe makes CIP longer than today, so the line is down for longer.")
st.dataframe(pd.DataFrame([base, new], index=["Today", "This recipe"])[["duration_min", "water_l", "naoh_kg", "fuel_kwh", "co2_kg", "cost_inr", "log_reduction"]].round(2), use_container_width=True)
st.caption("Kinetic, thermal and cost parameters are assumptions (see README); calibrate them with swab/ATP results and conductivity traces from a real circuit.")
