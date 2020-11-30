import numpy as np
import matplotlib.pyplot as plt
from simupy.block_diagram import DEFAULT_INTEGRATOR_OPTIONS
from scipy import interpolate
from simupy_flight import Planet
import glob
import pandas as pd
import argparse
import sys

int_opts = DEFAULT_INTEGRATOR_OPTIONS.copy()
int_opts['max_step'] = 2**-4

data_relative_path = '..'

ft_per_m = 3.28084
kg_per_slug = 14.5939

nesc_options = dict(
    interactive_mode = True,
    include_simupy_in_autoscale = True,
    only_baseline05 = False,
)

nesc_colors = {'%02d' % (sim_idx+1) : 'C%d' % (sim_idx) for sim_idx in range(10) }

parser = argparse.ArgumentParser()
parser.add_argument("--interactive", help="show plots in interactive mode rather than saving to disk",
                    action="store_true")
parser.add_argument("--simupy_scale", help="include simupy in plot autoscale or not",
                    action="store_true")
parser.add_argument("--baseline05", help="Use SIM 05 as the baseline rather than the ensemble (total) average",
                    action="store_true")

def do_argparse():
    if not hasattr(sys, 'ps1'):
        args = parser.parse_args()

        if args.interactive:
            nesc_options['interactive_mode'] = True
        else:
            nesc_options['interactive_mode'] = False

        if args.simupy_scale:
            nesc_options['include_simupy_in_autoscale'] = True
        else:
            nesc_options['include_simupy_in_autoscale'] = False
        
        if args.baseline05:
            nesc_options['only_baseline05'] = True
        else:
            nesc_options['only_baseline05'] = False


frame_rate_for_differencing = 10

def deg_diff(sim, baseline):
    """
    Compute the NESC angle difference, sim and baseline is a numpy array of angles, 
    note that the difference is not symmetric
    """
    beta = baseline*np.pi/180
    alpha = sim*np.pi/180
    
    diff_rad = np.arctan2(np.sin(alpha)*np.cos(beta) - np.cos(alpha)*np.sin(beta), np.cos(alpha)*np.cos(beta)+np.sin(alpha)*np.sin(beta))
    return diff_rad*180/np.pi
    

def plot_cols(simupy_res, baseline_pds, baseline_pd_labels, sfvt_idxs, baseline_cols, col_labels):
    """
    simupy_res is the simupy simulation result object
    baseline_pds is an iterable of pandas dataframes of the baseline simulation results
    baseline_pd_labels is an iterable of labels used in the figure legends 'SIM %s'
    sfvt_idxs is an iterable of integers of the column index to compare
    baseline_cols is the list of column names in the baseline pds to compare
    col_labels are the y-axis labels for each subplot
    """
    include_simupy_in_autoscale = nesc_options['include_simupy_in_autoscale']

    if len(sfvt_idxs) != len(baseline_cols) or len(sfvt_idxs) != len(col_labels):
        raise ValueError("Mismatched column lengths")

    num_plots = len(sfvt_idxs)

    abs_fig, abs_axes = plt.subplots(num_plots, sharex=True, constrained_layout=True)
    delta_fig, delta_axes = plt.subplots(num_plots, sharex=True, constrained_layout=True)

    tf = simupy_res.t[-1]
    num_for_average = int(tf*frame_rate_for_differencing)
    times_for_average = np.arange(0, tf, 1/frame_rate_for_differencing)

    for sfvt_idx, baseline_col, col_label, abs_ax, delta_ax in zip(sfvt_idxs, baseline_cols, col_labels, abs_axes, delta_axes):
        # collect and plot direct value of simupy result
        simupy_y = simupy_res.y[:, sfvt_idx]
        if 'deg' in baseline_col:
            simupy_y = simupy_y*180/np.pi
        elif 'ft' in baseline_col:
            simupy_y = simupy_y*ft_per_m

        # prepare to compute average
        num_of_averages = 0
        average_value = np.zeros(num_for_average)

        # iterate over NESC results
        for baseline_idx, baseline_pd in enumerate(baseline_pds):
            try: # try to  collect the columnr esult
                baseline_y = baseline_pd[baseline_col]
            except KeyError:
                print("missing %s for SIM %s" % (baseline_col, baseline_pd_labels[baseline_idx]))
            else: # if available, plot and add to ensemble average
                abs_ax.plot(baseline_pd.index, baseline_y, nesc_colors[baseline_pd_labels[baseline_idx]], alpha=0.5, label='NESC %s' % (baseline_pd_labels[baseline_idx]))
                to_interpolate = interpolate.make_interp_spline(baseline_pd.index, baseline_y)
                if (not nesc_options['only_baseline05']) or (baseline_pd_labels[baseline_idx] == '05'):
                    num_of_averages = num_of_averages + 1
                    average_value[:] = average_value[:] + to_interpolate(times_for_average)

        abs_ax_ylim = abs_ax.get_ylim()
        abs_ax.plot(simupy_res.t, simupy_y, 'k--', label='SimuPy')
        if not include_simupy_in_autoscale:
            abs_ax.set_ylim(*abs_ax_ylim)

        # compute average
        average_value[:] = average_value[:]/num_of_averages

        # compute difference of simupy result from average
        to_interpolate = interpolate.make_interp_spline(simupy_res.t, simupy_y)
        simupy_y = to_interpolate(times_for_average)
        if 'deg' in baseline_col:
            simupy_y = deg_diff(simupy_y, average_value)
        else:
            simupy_y = simupy_y - average_value

        for baseline_idx, baseline_pd in enumerate(baseline_pds): 
            try:
                baseline_y = baseline_pd[baseline_col]
            except KeyError:
                pass
            else:
                to_interpolate = interpolate.make_interp_spline(baseline_pd.index, baseline_y)
                baseline_y = to_interpolate(times_for_average)
                if 'deg' in baseline_col:
                    baseline_y = deg_diff(baseline_y, average_value)
                else:
                    baseline_y = baseline_y - average_value
                delta_ax.plot(times_for_average, baseline_y, nesc_colors[baseline_pd_labels[baseline_idx]], alpha=0.5, label='NESC %d' % (baseline_idx+1))

        delta_ax_ylim = delta_ax.get_ylim()
        delta_ax.plot(times_for_average, simupy_y, 'k--', label='SimuPy')
        if not include_simupy_in_autoscale:
            delta_ax.set_ylim(*delta_ax_ylim)

        abs_ax.set_ylabel(col_label)
        delta_ax.set_ylabel('$\\Delta$ ' + col_label)
        abs_ax.grid(True)
        delta_ax.grid(True)

    abs_axes[0].legend(ncol=3)

    abs_axes[-1].set_xlabel('time, s')
    delta_axes[-1].set_xlabel('time, s')

    return abs_fig, delta_fig

    
    

def plot_nesc_comparisons(simupy_res, baseline_pd_glob, plot_name=''):
    """
    """
    do_argparse()

    interactive_mode = nesc_options['interactive_mode']

    baseline_pds = []
    baseline_pd_labels = []

    prefix, suffix = baseline_pd_glob.split('*')
    for fname in sorted(glob.glob(baseline_pd_glob)):
        if interactive_mode:
            print(fname)
        baseline_pds.append(pd.read_csv(fname, index_col=0))
        baseline_pd_labels.append(fname.replace(prefix, '').replace(suffix, ''))


    abs_fig, delta_fig = plot_cols(simupy_res, baseline_pds, baseline_pd_labels,
    [Planet.lamda_D_idx, Planet.phi_D_idx, Planet.h_D_idx],
    ['longitude_deg', 'latitude_deg', 'altitudeMsl_ft'],
    ['longitude, deg', 'latitude, deg', 'altitude, ft'],
    )
    if not interactive_mode:
        abs_fig.set_size_inches(4, 6)
        delta_fig.set_size_inches(4, 6)

        abs_fig.savefig(plot_name + '_geodetic_pos.pdf')
        delta_fig.savefig(plot_name + '_geodetic_pos_delta.pdf')

    abs_fig, delta_fig = plot_cols(simupy_res, baseline_pds, baseline_pd_labels,
    [Planet.psi_idx, Planet.theta_idx, Planet.phi_idx],
    ['eulerAngle_deg_Yaw', 'eulerAngle_deg_Pitch', 'eulerAngle_deg_Roll'],
    #['yaw, deg', 'pitch, deg', 'roll, deg'],
    ['$\\psi$, deg', '$\\theta$, deg', '$\\phi$, deg']
    )

    if not interactive_mode:
        abs_fig.set_size_inches(4, 6)
        delta_fig.set_size_inches(4, 6)

        abs_fig.savefig(plot_name + '_eulerangle.pdf')
        delta_fig.savefig(plot_name + '_eulerangle_delta.pdf')

    abs_fig, delta_fig = plot_cols(simupy_res, baseline_pds, baseline_pd_labels,
    [Planet.omega_X_idx, Planet.omega_Y_idx, Planet.omega_Z_idx],
    ['bodyAngularRateWrtEi_deg_s_Roll', 'bodyAngularRateWrtEi_deg_s_Pitch', 'bodyAngularRateWrtEi_deg_s_Yaw'],
    ['$p$, deg/s', '$q$, deg/s', '$r$, deg/s']
    )

    if not interactive_mode:
        abs_fig.set_size_inches(4, 6)
        delta_fig.set_size_inches(4, 6)

        abs_fig.savefig(plot_name + '_body_rates.pdf')
        delta_fig.savefig(plot_name + '_body_rates_delta.pdf')

    abs_fig, delta_fig = plot_cols(simupy_res, baseline_pds, baseline_pd_labels,
    [Planet.p_x_idx, Planet.p_y_idx, Planet.p_z_idx],
    ['eiPosition_ft_X', 'eiPosition_ft_Y', 'eiPosition_ft_Z'],
    ['inertial $x$, ft', 'inertial $y$, ft', 'inertial $z$, ft']
    )

    if not interactive_mode:
        abs_fig.set_size_inches(4, 6)
        delta_fig.set_size_inches(4, 6)

        abs_fig.savefig(plot_name + '_inertial_pos.pdf')
        delta_fig.savefig(plot_name + '_inertial_pos_delta.pdf')

    abs_fig, delta_fig = plot_cols(simupy_res, baseline_pds, baseline_pd_labels,
    [Planet.V_N_idx, Planet.V_E_idx, Planet.V_D_idx],
    ['feVelocity_ft_s_X', 'feVelocity_ft_s_Y', 'feVelocity_ft_s_Z'],
    ['relative velocity\nN, ft/s', 'relative velocity\nE, ft/s', 'relative velocity\nD, ft/s']
    )

    if not interactive_mode:
        abs_fig.set_size_inches(4, 6)
        delta_fig.set_size_inches(4, 6)

        abs_fig.savefig(plot_name + '_velocity_NED.pdf')
        delta_fig.savefig(plot_name + '_velocity_NED_delta.pdf')

    if interactive_mode:
        plt.show()
