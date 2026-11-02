"""Reproduce every number, table and figure in README.md:   python run_analysis.py

1. Evaluate a typical fixed-time CIP recipe (baseline).
2. Find the recipe (temperature, NaOH concentration, wash time) with the lowest cost per
   cycle that still achieves a 3-log soil reduction, without making CIP longer than today,
   and switch the intermediate rinse from fixed time to a conductivity end-point.
3. Map the trade-off (cost heat-map, water-vs-energy Pareto front) and scale to a year.
"""
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from cipopt import Plant, Recipe, simulate, wash_time_needed, rinse_volume_to_endpoint

ROOT = Path(__file__).parent; FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
G, G2, CORAL, GOLD = "#1D5C38", "#58B97F", "#C9473C", "#A8843A"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
P = Plant()
BASE = Recipe()
CIP_PER_YEAR = {"orders as they arrive": 324, "optimised sequence (ChangeoverIQ)": 132}  # family changes per line per year


T_BOUNDS, C_BOUNDS = (60.0, 85.0), (0.75, 2.5)   # practical operating window (assumption)


def recipe_for(temp_c, conc_pct, recovery=False):
    return Recipe(temp_c=temp_c, conc_pct=conc_pct, wash_min=wash_time_needed(P, temp_c, conc_pct),
                  rinse1_l=rinse_volume_to_endpoint(P, conc_pct), recovery=recovery)


base = simulate(P, BASE)
print("baseline", {k: round(v, 2) if isinstance(v, float) else v for k, v in base.items()})

# ------------------------------------------------ grid (for the map and to seed / check the optimiser)
T = np.linspace(*T_BOUNDS, 51); C = np.linspace(*C_BOUNDS, 71)
cost = np.full((len(C), len(T)), np.nan); dur = np.full_like(cost, np.nan); fuel = np.full_like(cost, np.nan)
for i, c in enumerate(C):
    for j, t in enumerate(T):
        s = simulate(P, recipe_for(t, c)); cost[i, j] = s["cost_inr"]; dur[i, j] = s["duration_min"]; fuel[i, j] = s["fuel_kwh"]
feasible = (dur <= base["duration_min"] + 1e-9) & (fuel <= base["fuel_kwh"] + 1e-9)
masked = np.where(feasible, cost, np.inf)
i0, j0 = np.unravel_index(np.argmin(masked), masked.shape)

# ------------------------------------------------ continuous optimisation (SLSQP)
def f(x): return simulate(P, recipe_for(*x))["cost_inr"]
cons = [{"type": "ineq", "fun": lambda x: base["duration_min"] - simulate(P, recipe_for(*x))["duration_min"]},
        {"type": "ineq", "fun": lambda x: base["fuel_kwh"] - simulate(P, recipe_for(*x))["fuel_kwh"]}]
opt = minimize(f, x0=[T[j0], C[i0]], bounds=[T_BOUNDS, C_BOUNDS], constraints=cons, method="SLSQP")
t_opt, c_opt = float(opt.x[0]), float(opt.x[1])
best_r = recipe_for(t_opt, c_opt); chem = simulate(P, best_r)
best_r = recipe_for(t_opt, c_opt, recovery=True); best = simulate(P, best_r)
grid_best = float(masked.min())
print(f"optimum T={t_opt:.1f} C c={c_opt:.2f}% wash={best_r.wash_min:.1f} min chem-cost={chem['cost_inr']:.0f} (grid {grid_best:.0f}) with recovery {best['cost_inr']:.0f}")

# rinse-only change (keep today's chemistry, end-point rinse) to separate the two levers
rinse_only = simulate(P, replace(BASE, rinse1_l=rinse_volume_to_endpoint(P, BASE.conc_pct)))

def pct(a, b): return 100 * (1 - a / b)
keys = ["duration_min", "water_l", "naoh_kg", "fuel_kwh", "co2_kg", "cost_inr", "log_reduction", "residual_naoh_ppm"]
table = {"Today (fixed-time recipe)": base, "+ end-point rinse": rinse_only, "+ optimised T and NaOH": chem, "+ caustic recovery": best}
annual = {}
for label, n in CIP_PER_YEAR.items():
    annual[label] = {"cips": n,
                     "today": {k: base[k] * n for k in ("water_l", "naoh_kg", "fuel_kwh", "co2_kg", "cost_inr")},
                     "optimised": {k: best[k] * n for k in ("water_l", "naoh_kg", "fuel_kwh", "co2_kg", "cost_inr")}}
combined = {k: base[k] * CIP_PER_YEAR["orders as they arrive"] - best[k] * CIP_PER_YEAR["optimised sequence (ChangeoverIQ)"] for k in ("water_l", "naoh_kg", "fuel_kwh", "co2_kg", "cost_inr")}

# ------------------------------------------------ Pareto: water vs fuel (no duration limit), by sweeping T
pareto = []
for t in np.linspace(*T_BOUNDS, 26):
    cs = np.linspace(*C_BOUNDS, 176)
    sims = [simulate(P, recipe_for(t, c)) for c in cs]
    k = int(np.argmin([s["cost_inr"] for s in sims]))
    pareto.append((t, cs[k], sims[k]["water_l"], sims[k]["fuel_kwh"], sims[k]["duration_min"]))

summary = dict(
    baseline_recipe=BASE.__dict__, optimised_recipe=dict(temp_c=round(t_opt, 1), conc_pct=round(c_opt, 2), wash_min=round(best_r.wash_min, 1), rinse1_l=round(best_r.rinse1_l, 0)),
    per_cycle={k: {kk: round(v[kk], 3) for kk in keys} for k, v in table.items()},
    reduction_pct={kk: round(pct(best[kk], base[kk]), 1) for kk in ("water_l", "naoh_kg", "fuel_kwh", "co2_kg", "cost_inr")},
    optimiser=dict(success=bool(opt.success), slsqp_cost=round(chem["cost_inr"], 2), grid_cost=round(grid_best, 2)),
    annual_per_line=annual, combined_with_sequencing_saving_per_line_year={k: round(v, 1) for k, v in combined.items()},
    plant=P.__dict__)
(ROOT / "results.json").write_text(json.dumps(summary, indent=2, default=float))
for k, v in table.items():
    print(f"{k:40s}", {kk: round(v[kk], 2) for kk in keys})
print("reduction %", summary["reduction_pct"]); print("combined", summary["combined_with_sequencing_saving_per_line_year"])

# ------------------------------------------------ figures
fig, ax = plt.subplots(figsize=(8.5, 5))
im = ax.pcolormesh(T, C, cost, shading="auto", cmap="Greens_r")
ax.contour(T, C, dur, levels=[base["duration_min"]], colors=CORAL, linewidths=2)
ax.contour(T, C, fuel, levels=[base["fuel_kwh"]], colors=GOLD, linewidths=2, linestyles="--")
ax.text(60.5, 2.4, "red: CIP time = today  ·  gold dashed: boiler fuel = today", color="k", fontsize=8.5)
ax.plot(BASE.temp_c, BASE.conc_pct, "o", ms=10, mfc="white", mec=CORAL, mew=2, label="today's recipe")
ax.plot(t_opt, c_opt, "*", ms=18, color=GOLD, mec="k", label="optimum (no longer, no more fuel)")
fig.colorbar(im, label="Cost per CIP cycle (INR)")
ax.set_xlabel("Caustic temperature (°C)"); ax.set_ylabel("NaOH concentration (%)"); ax.legend(loc="lower right", frameon=True)
ax.set_title("Cost of a validated (3-log) CIP cycle", loc="left", fontweight="bold")
fig.tight_layout(); fig.savefig(FIG / "cost_map.png", dpi=160); plt.close(fig)

fig, axes = plt.subplots(1, 5, figsize=(13, 3.6))
for ax, (k, lab) in zip(axes, [("water_l", "Water (L)"), ("naoh_kg", "NaOH (kg)"), ("fuel_kwh", "Boiler fuel (kWh)"), ("co2_kg", "CO₂ (kg)"), ("duration_min", "Duration (min)")]):
    vals = [base[k], rinse_only[k], chem[k], best[k]]
    b = ax.bar(["Today", "+rinse", "+T/NaOH", "+recovery"], vals, color=[CORAL, G2, G2, G]); ax.bar_label(b, fmt="%.0f" if max(vals) > 20 else "%.1f", fontsize=7.5)
    ax.set_title(lab, fontsize=10); ax.tick_params(axis="x", labelsize=8, rotation=35)
fig.suptitle("Resources per CIP cycle", x=0.01, ha="left", fontweight="bold")
fig.tight_layout(); fig.savefig(FIG / "per_cycle.png", dpi=160); plt.close(fig)

fig, ax = plt.subplots(figsize=(8.5, 4))
V = np.linspace(0, 2000, 400)
ax.semilogy(V, BASE.conc_pct * 1e4 * np.exp(-V / P.holdup_l), color=G)
ax.axhline(P.endpoint_ppm, color=CORAL, ls="--"); ax.text(30, P.endpoint_ppm * 1.4, f"end-point {P.endpoint_ppm:.0f} ppm NaOH (read by conductivity)", color=CORAL, fontsize=8.5)
ax.axvline(BASE.rinse1_min * P.flow_lpm, color="k", ls=":"); ax.text(BASE.rinse1_min * P.flow_lpm - 20, 2e3, "fixed 10-min rinse", ha="right", fontsize=8.5)
ax.axvline(rinse_volume_to_endpoint(P, BASE.conc_pct), color=G, lw=2); ax.text(rinse_volume_to_endpoint(P, BASE.conc_pct) + 20, 2e3, "stop at end-point\n(+20% margin)", color=G, fontsize=8.5)
ax.set_xlabel("Rinse water used (L)"); ax.set_ylabel("Residual NaOH (ppm)"); ax.set_ylim(1e-2, 3e4)
ax.set_title("Intermediate rinse: stop on conductivity, not on the clock", loc="left", fontweight="bold")
fig.tight_layout(); fig.savefig(FIG / "rinse_endpoint.png", dpi=160); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.5, 4.3))
pw = np.array([q[2] for q in pareto]); pf = np.array([q[3] for q in pareto]); pt = np.array([q[0] for q in pareto])
sc = ax.scatter(pw, pf, c=pt, cmap="RdYlGn_r", s=40); fig.colorbar(sc, label="Wash temperature (°C)")
ax.plot(base["water_l"], base["fuel_kwh"], "o", ms=10, mfc="white", mec=CORAL, mew=2, label="today")
ax.plot(chem["water_l"], chem["fuel_kwh"], "*", ms=16, color=GOLD, mec="k", label="optimum chemistry (before recovery)")
ax.set_xlabel("Water per cycle (L)"); ax.set_ylabel("Boiler fuel per cycle (kWh)"); ax.legend(frameon=False)
ax.set_title("Water vs energy trade-off (cost-optimal NaOH at each temperature)", loc="left", fontweight="bold", fontsize=10)
fig.tight_layout(); fig.savefig(FIG / "pareto_water_energy.png", dpi=160); plt.close(fig)
print("done")
