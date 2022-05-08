"""
===================================================
Case 15: Circular F-16 flight around North Pole
===================================================

==============  ===============
Verifies        Propagation, geodetic transforms
Gravitation     J2
Geodesy         WGS-84 rotating
Atmosphere      US 1976 STD
Winds           still air
Vehicle         F-16 with circumnavigating auto-pilot
Notes           Initially straight & level and engage auto-pilot
==============  ===============

"""

from simupy import systems
from simupy.block_diagram import BlockDiagram
import os
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

from nesc_testcase_helper import (
    get_baselines,
    plot_nesc_comparisons,
    plot_F16_controls,
    nesc_options,
    nesc_colors,
    benchmark,
)
from nesc_testcase_helper import ft_per_m

from nesc_case11 import (
    int_opts,
    BD,
    earth,
    rho_0,
    eval_trim,
    run_trimmer,
    knots_per_mps,
)

import F16_model
from F16_gnc import F16_gnc, trimmedKEAS

F16_vehicle = F16_model.F16()

spec_ic_args = dict(
    phi_D=89.95 * np.pi / 180,  # latitude
    lamda_D=-45 * np.pi / 180,  # longitude
    h=10_000 / ft_per_m,
    V_N=0.0 / ft_per_m,
    V_E=563.643 / ft_per_m,
    V_D=0.0 / ft_per_m,
    psi=90.0 * np.pi / 180,
    theta=0.0 * np.pi / 180,
    phi=0.0 * np.pi / 180,
    p_B=0.0 * np.pi / 180,
    q_B=0.0 * np.pi / 180,
    r_B=0.0 * np.pi / 180,
)

earth_output_for_gnc_select = np.array(
    [
        earth.phi_D_idx,
        earth.lamda_D_idx,
        earth.h_D_idx,
        earth.V_T_idx,
        earth.alpha_idx,
        earth.beta_idx,
        earth.psi_idx,
        earth.theta_idx,
        earth.phi_idx,
        earth.p_B_idx,
        earth.q_B_idx,
        earth.r_B_idx,
        earth.rho_idx,
    ]
)


dim_feedback = len(earth_output_for_gnc_select)


def get_gnc_function(
    throttleTrim,
    longStkTrim,
    keasCmd=trimmedKEAS,
    altCmd=spec_ic_args["h"] * ft_per_m,
    sasOn=True,
    apOn=True,
    circlePoleSW=True,
):
    def gnc_function(t, u):
        throttle, longStk, latStk, pedal = 0.0, 0.0, 0.0, 0.0  # pilot command

        (
            latitude,
            longitude,
            alt,
            V_T,
            alpha,
            beta,
            psi,
            theta,
            phi,
            pb,
            qb,
            rb,  # feedback
            rho,
        ) = u  # rho to calculate equivalent airspeed

        Vequiv = V_T * np.sqrt(rho / rho_0)
        angles = np.array([latitude, longitude, alpha, beta, phi, theta, psi])
        latitude, longitude, alpha, beta, phi, theta, psi = angles * 180 / np.pi

        control_eart = F16_gnc(
            throttle,
            longStk,
            latStk,
            pedal,
            sasOn,
            apOn,
            circlePoleSW,
            latitude,
            longitude,
            keasCmd,
            altCmd,
            alt * ft_per_m,
            Vequiv * knots_per_mps,
            alpha,
            beta,
            phi,
            theta,
            psi,
            pb,
            qb,
            rb,
            throttleTrim,
            longStkTrim,
        )
        return control_eart

    return gnc_function


int_opts["nsteps"] = 100_000

# %%

if __name__ == "__main__":

    opt_args, opt_ctrl = run_trimmer(
        spec_ic_args, throttle_ic=0.0, longStk_ic=0.0, allow_roll=False
    )

    trimmed_flight_condition = earth.ic_from_planetodetic(**opt_args)

    trimmed_KEAS = (
        earth.output_equation_function(0, trimmed_flight_condition)[earth.V_T_idx]
        * np.sqrt(
            earth.output_equation_function(0, trimmed_flight_condition)[earth.rho_idx]
            / rho_0
        )
        * knots_per_mps
    )

    earth.initial_condition = trimmed_flight_condition

    gnc_block = systems.SystemFromCallable(
        get_gnc_function(
            *opt_ctrl,
            keasCmd=trimmed_KEAS,
        ),
        dim_feedback,
        4,
    )

    BD = BlockDiagram(earth, F16_vehicle, gnc_block)
    BD.connect(earth, F16_vehicle, inputs=np.arange(earth.dim_output))
    BD.connect(F16_vehicle, earth, inputs=np.arange(F16_vehicle.dim_output))
    BD.connect(
        gnc_block,
        F16_vehicle,
        inputs=np.arange(earth.dim_output, earth.dim_output + 4),
    )
    BD.connect(earth, gnc_block, outputs=earth_output_for_gnc_select)

    with benchmark() as b:
        res = BD.simulate(180, integrator_options=int_opts)

    plot_nesc_comparisons(res, "15")
    plot_F16_controls(res, "15", y_idx_offset=0)

    def xy_for_north_pole_ground_track(lat, long):
        xx = (90 - lat * 180 / np.pi) * np.cos(long)
        yy = (90 - lat * 180 / np.pi) * np.sin(long)
        return xx, yy

    plt.subplots(constrained_layout=True)
    plt.axis("equal")

    sim_lat = res.y[:, earth.phi_D_idx]
    sim_long = res.y[:, earth.lamda_D_idx]

    # refence_lats = np.array([90.0, 90-0.02, 90-0.04, 90-0.06])*np.pi/180
    refence_lats = (90 - np.arange(5) / 60) * np.pi / 180
    reference_longs = np.arange(0, 360, 30) * np.pi / 180
    ref_steps = 100

    ref_lag_long_list = []
    for i in range(len(refence_lats)):
        ref_lag_long_list.append(
            (refence_lats[i], np.arange(ref_steps + 1) * 2 * np.pi / ref_steps)
        )
    for j in range(len(reference_longs)):
        ref_lag_long_list.append(
            ((90 - np.arange(5) / 60) * np.pi / 180, reference_longs[j])
        )

    for lat, long in ref_lag_long_list:
        plt.plot(*xy_for_north_pole_ground_track(lat, long), "k--", alpha=0.25)

    baseline_pds, baseline_pd_labels = get_baselines("15")

    plt.plot(
        *xy_for_north_pole_ground_track(sim_lat, sim_long),
        "k-",
        alpha=1.0,
        label="SimuPy",
    )
    plt.plot(
        *xy_for_north_pole_ground_track(sim_lat[0], sim_long[0]),
        "o",
        alpha=0.5,
        markerfacecolor="None",
        markeredgecolor="k",
    )
    plt.plot(
        *xy_for_north_pole_ground_track(sim_lat[-1], sim_long[-1]),
        "x",
        alpha=0.5,
        markerfacecolor="None",
        markeredgecolor="k",
    )
    for baseline_idx, baseline_pd in enumerate(baseline_pds):
        plt.plot(
            *xy_for_north_pole_ground_track(
                *(baseline_pd[["latitude_deg", "longitude_deg"]] * np.pi / 180).values.T
            ),
            nesc_colors[baseline_pd_labels[baseline_idx]],
            alpha=0.5,
            label="NESC %s" % (baseline_pd_labels[baseline_idx]),
        )
        plt.plot(
            *xy_for_north_pole_ground_track(
                *(
                    baseline_pd[["latitude_deg", "longitude_deg"]].iloc[0] * np.pi / 180
                ).values.T
            ),
            "o",
            alpha=0.5,
            markerfacecolor="None",
            markeredgecolor=nesc_colors[baseline_pd_labels[baseline_idx]],
        )
        plt.plot(
            *xy_for_north_pole_ground_track(
                *(
                    baseline_pd[["latitude_deg", "longitude_deg"]].iloc[-1]
                    * np.pi
                    / 180
                ).values.T
            ),
            "x",
            alpha=0.5,
            markerfacecolor="None",
            markeredgecolor=nesc_colors[baseline_pd_labels[baseline_idx]],
        )

    plt.legend()

    if nesc_options["interactive_mode"]:
        plt.show()
    else:
        plt.savefig(
            os.path.join(nesc_options["save_relative_path"], "15_groundtrack.pdf")
        )
