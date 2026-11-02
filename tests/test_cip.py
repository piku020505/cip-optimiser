import numpy as np
from cipopt import Plant, Recipe, simulate, log_reduction, wash_time_needed, rinse_volume_to_endpoint

P = Plant()


def test_hotter_and_stronger_clean_faster():
    assert log_reduction(P, 80, 1.0, 20) > log_reduction(P, 70, 1.0, 20)
    assert log_reduction(P, 70, 2.0, 20) > log_reduction(P, 70, 1.0, 20)


def test_wash_time_meets_target_with_margin():
    t = wash_time_needed(P, 72, 1.2)
    assert np.isclose(log_reduction(P, 72, 1.2, t), P.log_target * P.safety)


def test_endpoint_rinse_meets_endpoint():
    r = Recipe(conc_pct=1.5, rinse1_l=rinse_volume_to_endpoint(P, 1.5))
    s = simulate(P, r)
    assert s["rinse_ok"] and s["residual_naoh_ppm"] < P.endpoint_ppm


def test_recovery_cuts_naoh_water_and_fuel():
    a = simulate(P, Recipe()); b = simulate(P, Recipe(recovery=True))
    assert b["naoh_kg"] < a["naoh_kg"] and b["water_l"] < a["water_l"] and b["fuel_kwh"] < a["fuel_kwh"]


def test_baseline_recipe_is_clean():
    assert simulate(P, Recipe())["clean_ok"]
