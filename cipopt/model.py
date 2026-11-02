"""Clean-in-place (CIP) model for a filling line's product circuit.

Cycle: pre-rinse -> caustic (NaOH) wash -> intermediate rinse -> final rinse.

* Soil removal in the caustic step is first order:  residual = exp(-k t),
  k = k_ref * exp(-Ea/R * (1/T - 1/T_ref)) * (c / c_ref)^n
  (temperature by Arrhenius, concentration by a power law).
* Rinsing a well-mixed circuit of hold-up V_mix:  C(V) = C0 * exp(-V / V_mix),
  so the rinse water needed to reach an end-point concentration is V_mix * ln(C0 / C_end).
* Heating: the caustic charge is heated from ambient (m cp dT, boiler efficiency),
  then held at temperature against heat loss UA (T - T_amb).

Every parameter below is an assumption chosen to be in a typical range for a small
cosmetics filling circuit; replace them with measured values (flow, hold-up, swab or
ATP results, conductivity traces) before using the numbers.
"""
from dataclasses import dataclass, replace
import numpy as np

R = 8.314  # J/mol/K


@dataclass(frozen=True)
class Plant:
    flow_lpm: float = 180.0        # ~1.5 m/s in DN50 pipe
    holdup_l: float = 150.0        # well-mixed volume of the circuit (V_mix)
    caustic_charge_l: float = 400.0
    t_amb_c: float = 25.0
    ua_w_per_k: float = 250.0      # heat loss of tank + circuit while holding temperature
    boiler_eff: float = 0.85
    # cleaning kinetics (caustic step)
    k_ref_per_min: float = 0.20    # at 70 C and 1.0 % NaOH
    t_ref_c: float = 70.0
    c_ref_pct: float = 1.0
    ea_j_mol: float = 40_000.0
    n_conc: float = 0.5
    log_target: float = 3.0        # required soil reduction (3-log = 99.9 %)
    safety: float = 1.2            # margin on wash time and rinse volume
    endpoint_ppm: float = 10.0     # residual NaOH allowed at end of intermediate rinse
    # caustic recovery (re-use the wash solution from a recovery tank)
    reuse_cycles: int = 5            # a charge is dumped after this many cycles
    topup_frac: float = 0.10         # make-up for carry-over and dilution each cycle
    return_temp_c: float = 55.0      # recovered caustic returns warm
    # costs and emissions (assumptions, INR)
    water_inr_per_kl: float = 60.0     # supply + effluent treatment
    naoh_inr_per_kg: float = 45.0
    heat_inr_per_kwh: float = 3.5      # boiler fuel per kWh of heat delivered to fuel
    co2_kg_per_kwh_fuel: float = 0.20  # natural-gas-fired boiler


@dataclass(frozen=True)
class Recipe:
    prerinse_min: float = 5.0
    temp_c: float = 75.0
    conc_pct: float = 1.5
    wash_min: float = 30.0
    rinse1_l: float = None   # None = fixed-time rinse (rinse1_min); a number = end-point rinse volume
    rinse1_min: float = 10.0
    final_min: float = 10.0
    recovery: bool = False   # caustic recovery tank in use


def k_clean(p: Plant, temp_c, conc_pct):
    T, Tr = temp_c + 273.15, p.t_ref_c + 273.15
    return p.k_ref_per_min * np.exp(-p.ea_j_mol / R * (1 / T - 1 / Tr)) * (conc_pct / p.c_ref_pct) ** p.n_conc


def log_reduction(p: Plant, temp_c, conc_pct, wash_min):
    return k_clean(p, temp_c, conc_pct) * wash_min / np.log(10)


def wash_time_needed(p: Plant, temp_c, conc_pct):
    """Minutes of caustic contact for the target log reduction, with the safety margin."""
    return p.safety * p.log_target * np.log(10) / k_clean(p, temp_c, conc_pct)


def rinse_volume_to_endpoint(p: Plant, conc_pct):
    c0_ppm = conc_pct * 10_000
    return p.safety * p.holdup_l * np.log(c0_ppm / p.endpoint_ppm)


def simulate(p: Plant, r: Recipe):
    """Resources, time, cost and CO2 for one CIP cycle."""
    rinse1_l = r.rinse1_l if r.rinse1_l is not None else r.rinse1_min * p.flow_lpm
    rinse1_min = rinse1_l / p.flow_lpm
    if r.recovery:
        fresh_l = p.caustic_charge_l / p.reuse_cycles + p.topup_frac * p.caustic_charge_l
        t_start = p.return_temp_c
    else:
        fresh_l, t_start = p.caustic_charge_l, p.t_amb_c
    water_l = (r.prerinse_min + r.final_min) * p.flow_lpm + rinse1_l + fresh_l
    naoh_kg = fresh_l * r.conc_pct / 100  # ~1 kg/L solution
    m_kg = p.caustic_charge_l
    heat_kwh = (m_kg * 4.18 * max(r.temp_c - t_start, 0) / 3600 + p.ua_w_per_k * (r.temp_c - p.t_amb_c) * r.wash_min * 60 / 3.6e6)
    fuel_kwh = heat_kwh / p.boiler_eff
    residual_ppm = r.conc_pct * 10_000 * np.exp(-rinse1_l / p.holdup_l)
    out = dict(
        log_reduction=float(log_reduction(p, r.temp_c, r.conc_pct, r.wash_min)),
        residual_naoh_ppm=float(residual_ppm),
        duration_min=float(r.prerinse_min + r.wash_min + rinse1_min + r.final_min),
        water_l=float(water_l), naoh_kg=float(naoh_kg), heat_kwh=float(heat_kwh), fuel_kwh=float(fuel_kwh),
        co2_kg=float(fuel_kwh * p.co2_kg_per_kwh_fuel),
    )
    out["cost_inr"] = float(water_l / 1000 * p.water_inr_per_kl + naoh_kg * p.naoh_inr_per_kg + fuel_kwh * p.heat_inr_per_kwh)
    out["clean_ok"] = out["log_reduction"] >= p.log_target
    out["rinse_ok"] = out["residual_naoh_ppm"] <= p.endpoint_ppm
    return out
