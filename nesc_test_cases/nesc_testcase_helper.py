import numpy as np
import matplotlib.pyplot as plt
from simupy.block_diagram import DEFAULT_INTEGRATOR_OPTIONS, SimulationResult
from simupy.utils import isclose
from scipy import interpolate
from simupy_flight import Planet, Vehicle
import glob
import pandas as pd
import argparse
import sys
import os
import time
from contextlib import contextmanager
from dataclasses import dataclass

ft_per_m = 3.28084
kg_per_slug = 14.5939
N_per_lbf = 4.44822

int_opts = DEFAULT_INTEGRATOR_OPTIONS.copy()
int_opts["max_step"] = 2**-4

frame_rate_for_differencing = 10

_here = os.path.abspath(os.path.dirname(__file__))
regression_data_path = os.path.join(_here, "..", "regression_data")

nesc_options = dict(
    interactive_mode=True,
    include_simupy_in_autoscale=True,
    only_baseline05=False,
    data_relative_path=os.path.join(_here, "..", "NESC_data"),
    save_relative_path="plots/",
    regression_test=True,
    write_regression_data=False,
)

nesc_colors = {"%02d" % (sim_idx + 1): "C%d" % (sim_idx) for sim_idx in range(10)}

# if not in an interactive interpreter session, handle command line args
if not hasattr(sys, "ps1"):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        help="show plots in interactive mode rather than saving to disk",
    )
    parser.add_argument(
        "--simupy-scale",
        action="store_true",
        help="include simupy in plot autoscale or not",
    )
    parser.add_argument(
        "--baseline05",
        action="store_true",
        help="Use SIM 05 as the baseline rather than the ensemble (total) average",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        help="Path to save plot outputs",
    )
    parser.add_argument(
        "--nesc-data-path",
        help="Path to parent directory of Atmospheric_checkcases data folder",
    )
    parser.add_argument(
        "-w",
        "--write-regression-data",
        action="store_true",
        help="Write simulation results to regression data path",
    )
    parser.add_argument(
        "--no-test",
        action="store_true",
        help="skip regression testing",
    )

    args = parser.parse_args()

    nesc_options["interactive_mode"] = args.interactive
    nesc_options["include_simupy_in_autoscale"] = args.simupy_scale
    nesc_options["only_baseline05"] = args.baseline05

    if args.nesc_data_path:
        nesc_options["data_relative_path"] = args.nesc_data_path
    if args.output_path:
        nesc_options["save_relative_path"] = args.output_path

    if not os.path.exists(nesc_options["save_relative_path"]):
        os.makedirs(nesc_options["save_relative_path"])

    if args.write_regression_data:
        nesc_options["write_regression_data"] = True
        nesc_options["regression_test"] = False

        if not os.path.exists(regression_data_path):
            os.mkdir(regression_data_path)
    else:
        nesc_options["regression_test"] = not args.no_test


def deg_diff(sim, baseline):
    """
    Compute the NESC angle difference, sim and baseline is a numpy array of angles,
    note that the difference is not symmetric
    """
    beta = baseline * np.pi / 180
    alpha = sim * np.pi / 180

    diff_rad = np.arctan2(
        np.sin(alpha) * np.cos(beta) - np.cos(alpha) * np.sin(beta),
        np.cos(alpha) * np.cos(beta) + np.sin(alpha) * np.sin(beta),
    )
    return diff_rad * 180 / np.pi


def plot_cols(
    simupy_res, baseline_pds, baseline_pd_labels, sfvt_idxs, baseline_cols, col_labels
):
    """
    simupy_res is the simupy simulation result object
    baseline_pds is an iterable of pandas dataframes of the baseline simulation results
    baseline_pd_labels is an iterable of labels used in the figure legends 'SIM %s'
    sfvt_idxs is an iterable of integers of the column index to compare
    baseline_cols is the list of column names in the baseline pds to compare
    col_labels are the y-axis labels for each subplot
    """
    include_simupy_in_autoscale = nesc_options["include_simupy_in_autoscale"]

    if len(sfvt_idxs) != len(baseline_cols) or len(sfvt_idxs) != len(col_labels):
        raise ValueError("Mismatched column lengths")

    num_plots = len(sfvt_idxs)

    abs_fig, abs_axes = plt.subplots(num_plots, sharex=True, constrained_layout=True)
    if len(baseline_pds) > 0:
        delta_fig, delta_axes = plt.subplots(
            num_plots, sharex=True, constrained_layout=True
        )
    else:
        delta_axes = abs_axes

    tf = simupy_res.t[-1]
    num_for_average = int(tf * frame_rate_for_differencing) + 1
    times_for_average = np.arange(num_for_average) / frame_rate_for_differencing

    for sfvt_idx, baseline_col, col_label, abs_ax, delta_ax in zip(
        sfvt_idxs, baseline_cols, col_labels, abs_axes, delta_axes
    ):
        # collect and plot direct value of simupy result
        simupy_y = simupy_res.y[:, sfvt_idx]
        if "deg" in baseline_col:
            simupy_y = simupy_y * 180 / np.pi
        elif "ft" in baseline_col:
            simupy_y = simupy_y * ft_per_m

        # prepare to compute average
        num_of_averages = 0
        average_value = np.zeros(num_for_average)

        # iterate over NESC results
        for baseline_idx, baseline_pd in enumerate(baseline_pds):
            try:  # try to  collect the columnr esult
                baseline_y = baseline_pd[baseline_col]
            except KeyError:
                pass
            else:  # if available, plot and add to ensemble average
                plot_baseline_sel = baseline_pd.index <= tf
                baseline_t = baseline_pd.index[plot_baseline_sel]
                abs_ax.plot(
                    baseline_t,
                    baseline_y[plot_baseline_sel],
                    nesc_colors[baseline_pd_labels[baseline_idx]],
                    alpha=0.5,
                    label="NESC %s" % (baseline_pd_labels[baseline_idx]),
                )
                to_interpolate = interpolate.make_interp_spline(
                    baseline_pd.index, baseline_y
                )
                if (not nesc_options["only_baseline05"]) or (
                    baseline_pd_labels[baseline_idx] == "05"
                ):
                    num_of_averages = num_of_averages + 1
                    average_value[:] = average_value[:] + to_interpolate(
                        times_for_average
                    )

        abs_ax_ylim = abs_ax.get_ylim()
        abs_ax.plot(simupy_res.t, simupy_y, "k--", label="SimuPy")
        if not include_simupy_in_autoscale:
            abs_ax.set_ylim(*abs_ax_ylim)

        # compute average
        if num_of_averages > 0:
            average_value[:] = average_value[:] / num_of_averages

        # compute difference of simupy result from average
        to_interpolate = interpolate.make_interp_spline(simupy_res.t, simupy_y)
        simupy_y = to_interpolate(times_for_average)
        if "deg" in baseline_col:
            simupy_y = deg_diff(simupy_y, average_value)
        else:
            simupy_y = simupy_y - average_value

        for baseline_idx, baseline_pd in enumerate(baseline_pds):
            try:
                baseline_y = baseline_pd[baseline_col]
            except KeyError:
                pass
            else:
                to_interpolate = interpolate.make_interp_spline(
                    baseline_pd.index, baseline_y
                )
                baseline_y = to_interpolate(times_for_average)
                if "deg" in baseline_col:
                    baseline_y = deg_diff(baseline_y, average_value)
                else:
                    baseline_y = baseline_y - average_value
                delta_ax.plot(
                    times_for_average,
                    baseline_y,
                    nesc_colors[baseline_pd_labels[baseline_idx]],
                    alpha=0.5,
                    label="NESC %d" % (baseline_idx + 1),
                )

        if len(baseline_pds) > 0:
            delta_ax_ylim = delta_ax.get_ylim()
            delta_ax.plot(times_for_average, simupy_y, "k--", label="SimuPy")
            if not include_simupy_in_autoscale:
                delta_ax.set_ylim(*delta_ax_ylim)
            delta_ax.set_ylabel("$\\Delta$ " + col_label)
            delta_ax.grid(True)

        abs_ax.set_ylabel(col_label)

        abs_ax.grid(True)

    abs_axes[0].legend(ncol=2)
    abs_axes[-1].set_xlabel("time, s")

    if len(baseline_pds) > 0:
        delta_axes[-1].set_xlabel("time, s")
        return abs_fig, delta_fig

    return abs_fig, None


def get_baselines(case):
    interactive_mode = nesc_options["interactive_mode"]

    baseline_pds = []
    baseline_pd_labels = []

    glob_path = os.path.join(
        nesc_options["data_relative_path"],
        "Atmospheric_checkcases",
        f"Atmos_{case}_*",
        f"Atmos_{case}_sim_*.csv",
    )
    for fname in sorted(glob.glob(glob_path)):
        sim = fname.rsplit("_", maxsplit=1)[-1].replace(".csv", "")
        if interactive_mode:
            print(fname)
        baseline_pds.append(pd.read_csv(fname, index_col=0))
        baseline_pd_labels.append(sim)

    return baseline_pds, baseline_pd_labels


def regression_test(res, case):
    fp = os.path.join(regression_data_path, f"case_{case}.npz")

    if nesc_options["write_regression_data"]:
        res.to_file(fp)
        print(f"Outputs saved to {os.path.basename(fp)}")
        return

    if not os.path.exists(fp):
        print("No regression file found, skipping test")
        return

    test_res = SimulationResult.from_file(fp)
    cols = slice(7) # translational position and quaternion state
    if case == "11":
        atol = 1.E-8
        rtol = 1.E-9
        p = np.inf
    elif case in ["13p1", "13p2", "13p3", "13p4"]:
        atol = 2.5E-7
        rtol = 1.E-9
        p = 2
    elif case in ["15", "16"]:
        atol = 1.E-7
        rtol = 1.E-7
        p = 2
    else:
        atol = 1.E-12
        rtol = 1.E-9
        p = np.inf

    matches = isclose(test_res, res, cols=cols, p=p, atol=atol, rtol=rtol, mode="pep485")
    passed = np.all(matches)
    result = "passed" if passed else "failed"
    print(f"Regression test {result.upper()}")

    if not passed:
        for i, match in enumerate(matches):
            if not match:
                print(f"col {i:02d} failed")
        print("Writing output data for comparison")
        fp = os.path.join(regression_data_path, f"case_{case}_fail.npz")
        res.to_file(fp)


def plot_nesc_comparisons(simupy_res, case, plot_name=""):
    if nesc_options["regression_test"]:
        regression_test(simupy_res, case)

    if plot_name == "":
        plot_name = case

    save_relative_path = nesc_options["save_relative_path"]
    interactive_mode = nesc_options["interactive_mode"]

    baseline_pds, baseline_pd_labels = get_baselines(case)

    abs_fig, delta_fig = plot_cols(
        simupy_res,
        baseline_pds,
        baseline_pd_labels,
        [Planet.lamda_D_idx, Planet.phi_D_idx, Planet.h_D_idx],
        ["longitude_deg", "latitude_deg", "altitudeMsl_ft"],
        ["longitude, deg", "latitude, deg", "altitude, ft"],
    )
    if not interactive_mode:
        abs_fig.set_size_inches(4, 6)
        abs_fig.savefig(
            os.path.join(save_relative_path, f"{plot_name}_geodetic_pos.pdf")
        )

        if delta_fig is not None:
            delta_fig.set_size_inches(4, 6)
            delta_fig.savefig(
                os.path.join(save_relative_path, f"{plot_name}_geodetic_pos_delta.pdf")
            )

    abs_fig, delta_fig = plot_cols(
        simupy_res,
        baseline_pds,
        baseline_pd_labels,
        [Planet.psi_idx, Planet.theta_idx, Planet.phi_idx],
        ["eulerAngle_deg_Yaw", "eulerAngle_deg_Pitch", "eulerAngle_deg_Roll"],
        # ['yaw, deg', 'pitch, deg', 'roll, deg'],
        ["$\\psi$, deg", "$\\theta$, deg", "$\\phi$, deg"],
    )

    if not interactive_mode:
        abs_fig.set_size_inches(4, 6)
        abs_fig.savefig(os.path.join(save_relative_path, plot_name + "_eulerangle.pdf"))

        if delta_fig is not None:
            delta_fig.set_size_inches(4, 6)
            delta_fig.savefig(
                os.path.join(save_relative_path, plot_name + "_eulerangle_delta.pdf")
            )

    abs_fig, delta_fig = plot_cols(
        simupy_res,
        baseline_pds,
        baseline_pd_labels,
        [Planet.omega_X_idx, Planet.omega_Y_idx, Planet.omega_Z_idx],
        [
            "bodyAngularRateWrtEi_deg_s_Roll",
            "bodyAngularRateWrtEi_deg_s_Pitch",
            "bodyAngularRateWrtEi_deg_s_Yaw",
        ],
        ["$p$, deg/s", "$q$, deg/s", "$r$, deg/s"],
    )

    if not interactive_mode:
        abs_fig.set_size_inches(4, 6)
        abs_fig.savefig(os.path.join(save_relative_path, plot_name + "_body_rates.pdf"))

        if delta_fig is not None:
            delta_fig.set_size_inches(4, 6)
            delta_fig.savefig(
                os.path.join(save_relative_path, plot_name + "_body_rates_delta.pdf")
            )

    abs_fig, delta_fig = plot_cols(
        simupy_res,
        baseline_pds,
        baseline_pd_labels,
        [Planet.p_x_idx, Planet.p_y_idx, Planet.p_z_idx],
        ["eiPosition_ft_X", "eiPosition_ft_Y", "eiPosition_ft_Z"],
        ["inertial $x$, ft", "inertial $y$, ft", "inertial $z$, ft"],
    )

    if not interactive_mode:
        abs_fig.set_size_inches(4, 6)
        abs_fig.savefig(
            os.path.join(save_relative_path, plot_name + "_inertial_pos.pdf")
        )

        if delta_fig is not None:
            delta_fig.set_size_inches(4, 6)
            delta_fig.savefig(
                os.path.join(save_relative_path, plot_name + "_inertial_pos_delta.pdf")
            )

    abs_fig, delta_fig = plot_cols(
        simupy_res,
        baseline_pds,
        baseline_pd_labels,
        [Planet.V_N_idx, Planet.V_E_idx, Planet.V_D_idx],
        ["feVelocity_ft_s_X", "feVelocity_ft_s_Y", "feVelocity_ft_s_Z"],
        [
            "relative velocity\nN, ft/s",
            "relative velocity\nE, ft/s",
            "relative velocity\nD, ft/s",
        ],
    )

    if not interactive_mode:
        abs_fig.set_size_inches(4, 6)
        abs_fig.savefig(
            os.path.join(save_relative_path, plot_name + "_velocity_NED.pdf")
        )

        if delta_fig is not None:
            delta_fig.set_size_inches(4, 6)
            delta_fig.savefig(
                os.path.join(save_relative_path, plot_name + "_velocity_NED_delta.pdf")
            )

    if interactive_mode:
        plt.show()


def plot_F16_controls(simupy_res, plot_name="", y_idx_offset=-4):
    """ """
    save_relative_path = nesc_options["save_relative_path"]
    interactive_mode = nesc_options["interactive_mode"]
    old_scale_opt = nesc_options["include_simupy_in_autoscale"]
    nesc_options["include_simupy_in_autoscale"] = True

    abs_fig, delta_fig = plot_cols(
        simupy_res,
        [],
        [],
        np.array([-4, -3, -2, -1]) + y_idx_offset,
        [
            "",
        ]
        * 4,  # no NESC test cases provide the control output!
        ["elevator, deg", "aileron, deg", "rudder, deg", "throttle, %"],
    )

    if not old_scale_opt:

        for ax, lims in zip(
            abs_fig.axes, ((-30, 30), (-26.5, 26.5), (-35, 35), (-5, 105))
        ):
            ax.set_ylim(*lims)

    nesc_options["include_simupy_in_autoscale"] = old_scale_opt

    if not interactive_mode:
        abs_fig.set_size_inches(4, 6)
        abs_fig.savefig(os.path.join(save_relative_path, plot_name + "_controls.pdf"))

    if interactive_mode:
        plt.show()


@dataclass
class BenchmarkInfo:
    """Class data for passing timing data to an enclosing benchmark context."""

    tfinal: float = None


@contextmanager
def benchmark():
    """Context manager for timing and printing runtime of code within the context.

    A ``BenchmarkInfo`` object is yielded so the enclosed code block can pass
    information back to the context manager for printing.
    """
    b = BenchmarkInfo()
    ts = time.time()

    yield b

    dt = time.time() - ts

    eval_msg = ""
    if b.tfinal is not None:
        r = b.tfinal / dt
        eval_msg = f"    eval time to run time: {r:.3f}"

    print(f"time to simulate: {dt:.3f} s{eval_msg}")
