"""Numerical baseline from the public nprobust illustration data."""

from pathlib import Path

import numpy as np
import pandas as pd

from nprobust import kdbwselect, kdrobust, lpbwselect, lprobust


DATA = Path(__file__).resolve().parents[2] / "nprobust_data.csv"


def illustration_data():
    data = pd.read_csv(DATA)
    t = data["t"].to_numpy()
    return {
        "t": t,
        "chol1": data["chol1"].to_numpy(),
        "cholf": data["cholf"].to_numpy(),
    }


def test_density_estimates_match_illustration_baseline():
    data = illustration_data()
    chol1 = data["chol1"]
    eval_grid = np.linspace(chol1.min(), chol1.max(), 5)
    result = kdrobust(chol1[data["t"] == 0], eval=eval_grid, bwselect="imse-dpi")

    np.testing.assert_allclose(
        result.Estimate[["tau.us", "tau.bc", "se.rb"]].to_numpy(dtype=float),
        np.array(
            [
                [0.0025859566367191712, 0.0017439604302553788, 0.0004959976725936008],
                [0.015422146044751598, 0.016838878505664023, 0.0011309750484081325],
                [0.0016492713323383573, 0.0014283335541471313, 0.00046374530178708355],
                [0.0007608018428344705, 0.0006844576992744829, 0.0002059670617245624],
                [0.00038851230966021114, 0.00031530704091149026, 0.00011351804054311],
            ]
        ),
        rtol=1e-10,
        atol=1e-12,
    )


def test_local_polynomial_estimates_match_illustration_baseline():
    data = illustration_data()
    chol1 = data["chol1"]
    ev = np.linspace(np.quantile(chol1, 0.10), np.quantile(chol1, 0.90), 5)
    result = lprobust(data["cholf"][data["t"] == 0], chol1[data["t"] == 0], eval=ev)

    np.testing.assert_allclose(
        result.Estimate[["tau.us", "tau.bc", "se.rb"]].to_numpy(dtype=float),
        np.array(
            [
                [255.95286891586267, 256.76417242564287, 2.040838516518859],
                [270.64313126572574, 269.6318188272758, 1.515342825504498],
                [286.1054829343728, 284.9600508390116, 1.7321051671367078],
                [301.4668753229664, 300.8733226018188, 2.3485973922394],
                [316.15087913019784, 316.38090877070096, 3.323344234485067],
            ]
        ),
        rtol=1e-10,
        atol=1e-9,
    )


def test_bandwidths_match_illustration_baseline():
    data = illustration_data()
    chol1 = data["chol1"]
    ev = np.linspace(np.quantile(chol1, 0.10), np.quantile(chol1, 0.90), 5)

    kd_bw = kdbwselect(chol1[data["t"] == 0], bwselect="imse-dpi")
    lp_bw = lpbwselect(
        data["cholf"][data["t"] == 0],
        chol1[data["t"] == 0],
        eval=ev,
        bwselect="mse-dpi",
    )

    np.testing.assert_allclose(
        kd_bw.bws[["h", "b"]].iloc[0].to_numpy(dtype=float),
        np.array([29.291065515881648, 76.06998841086319]),
        rtol=1e-10,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        lp_bw.bws[["h", "b"]].to_numpy(dtype=float),
        np.array(
            [
                [34.27342329281759, 70.69875419601624],
                [42.808199958599815, 157.79999999999995],
                [28.898171976015586, 77.16454879573325],
                [34.646693757172145, 83.67941185241095],
                [66.72734012471214, 100.64110885754822],
            ]
        ),
        rtol=1e-10,
        atol=1e-10,
    )
