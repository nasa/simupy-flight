import numpy as np
import matplotlib.pyplot as plt
from simupy.block_diagram import DEFAULT_INTEGRATOR_OPTIONS
from scipy import interpolate
from simupy_flight import Planet

int_opts = DEFAULT_INTEGRATOR_OPTIONS.copy()
int_opts['max_step'] = 2**-4

data_relative_path = '..'

ft_per_m = 3.28084
kg_per_slug = 14.5939

frame_rate_for_differencing = 10

def deg_diff(sim, baseline):
    beta = baseline*np.pi/180
    alpha = sim*np.pi/180
    
    diff_rad = np.arctan2(np.sin(alpha)*np.cos(beta) - np.cos(alpha)*np.sin(beta), np.cos(alpha)*np.cos(beta)+np.sin(alpha)*np.sin(beta))
    return diff_rad*180/np.pi
    

def plot_cols(simupy_res, baseline_pds, baseline_cols, sfvt_idxs, labels):
    
    abs_fig, abs_axes = plt.subplots(3, sharex=True, constrained_layout=True)
    delta_fig, delta_axes = plt.subplots(3, sharex=True, constrained_layout=True)
    
    
    tf = simupy_res.t[-1]
    num_for_average = int(tf*frame_rate_for_differencing)
    times_for_average = np.arange(0, tf, 1/frame_rate_for_differencing)
    
    
    for sfvt_idx, baseline_col, label, abs_ax, delta_ax in zip(baseline_cols, sfvt_idxs, labels, abs_axes, delta_axes):
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
                print("missing %s for sim %d" % (baseline_col, baseline_idx))
            else: # if available, plot and add to ensemble average
                abs_ax.plot(baseline_pd.index, baseline_y, 'C%d'%baseline_idx, alpha=0.5, label='NESC %d' % (baseline_idx+1))
                num_of_averages = num_of_averages + 1
                to_interpolate = interpolate.make_interp_spline(baseline_pd.index, baseline_y)
                average_value[:] = average_value[:] + to_interpolate(times_for_average)

        abs_ax_ylim = abs_ax.get_ylim()
        abs_ax.plot(simupy_res.t, simupy_y, 'k--', label='SimuPy')
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
                delta_ax.plot(times_for_average, baseline_y, 'C%d'%baseline_idx, alpha=0.5, label='NESC %d' % (baseline_idx+1))
        
        delta_ax_ylim = delta_ax.get_ylim()
        delta_ax.plot(times_for_average, simupy_y, 'k--', label='SimuPy')
        delta_ax.set_ylim(*delta_ax_ylim)
        
        abs_ax.set_ylabel(label)
        delta_ax.set_ylabel('$\\Delta$ ' + label)
        abs_ax.grid(True)
        delta_ax.grid(True)
        
    abs_axes[0].legend(ncol=3)
    
    abs_axes[-1].set_xlabel('time, s')
    delta_axes[-1].set_xlabel('time, s')
    
    return abs_fig, delta_fig

    
    

def plot_nesc_comparisons(simupy_res, baseline_pds, name=''):
    abs_fig, delta_fig = plot_cols(simupy_res, baseline_pds,
    [Planet.lamda_D_idx, Planet.phi_D_idx, Planet.h_D_idx],
    ['longitude_deg', 'latitude_deg', 'altitudeMsl_ft'],
    ['longitude, deg', 'latitude, deg', 'altitude, ft'],
    )

    abs_fig.set_size_inches(4, 6)
    delta_fig.set_size_inches(4, 6)

    abs_fig.savefig(name + 'geodetic_pos.pdf')
    delta_fig.savefig(name + 'geodetic_pos_delta.pdf')
    
    abs_fig, delta_fig = plot_cols(simupy_res, baseline_pds,
    [Planet.psi_idx, Planet.theta_idx, Planet.phi_idx],
    ['eulerAngle_deg_Yaw', 'eulerAngle_deg_Pitch', 'eulerAngle_deg_Roll'],
    #['yaw, deg', 'pitch, deg', 'roll, deg'],
    ['$\\psi$, deg', '$\\theta$, deg', '$\\phi$, deg']
    )

    abs_fig.set_size_inches(4, 6)
    delta_fig.set_size_inches(4, 6)

    abs_fig.savefig(name + 'eulerangle.pdf')
    delta_fig.savefig(name + 'eulerangle_delta.pdf')
    
    abs_fig, delta_fig = plot_cols(simupy_res, baseline_pds,
    # [Planet.p_B_idx, Planet.q_B_idx, Planet.r_B_idx],
    [Planet.omega_X_idx, Planet.omega_Y_idx, Planet.omega_Z_idx],
    ['bodyAngularRateWrtEi_deg_s_Roll', 'bodyAngularRateWrtEi_deg_s_Pitch', 'bodyAngularRateWrtEi_deg_s_Yaw'],
    ['$p$, deg/s', '$q$, deg/s', '$r$, deg/s']
    )

    abs_fig.set_size_inches(4, 6)
    delta_fig.set_size_inches(4, 6)

    abs_fig.savefig(name + 'body_rates.pdf')
    delta_fig.savefig(name + 'body_rates_delta.pdf')
    
    # abs_fig, delta_fig = plot_cols(simupy_res, baseline_pds,
    # [Planet.p_x_idx, Planet.p_y_idx, Planet.p_z_idx],
    # ['eiPosition_ft_X', 'eiPosition_ft_Y', 'eiPosition_ft_X'],
    # ['inertial $x$, ft', 'inertial $y$, ft', 'inertial $z$, ft']
    # )

    # abs_fig.set_size_inches(4, 6)
    # delta_fig.set_size_inches(4, 6)

    # abs_fig.savefig(name + 'inertial_pos.pdf')
    # delta_fig.savefig(name + 'inertial_pos_delta.pdf')
    
    abs_fig, delta_fig = plot_cols(simupy_res, baseline_pds,
    [Planet.V_N_idx, Planet.V_E_idx, Planet.V_D_idx],
    ['feVelocity_ft_s_X', 'feVelocity_ft_s_Y', 'feVelocity_ft_s_Z'],
    ['relative velocity\nN, ft/s', 'relative velocity\nE, ft/s', 'relative velocity\nD, ft/s']
    )

    abs_fig.set_size_inches(4, 6)
    delta_fig.set_size_inches(4, 6)

    abs_fig.savefig(name + 'velocity_NED.pdf')
    delta_fig.savefig(name + 'velocity_NED_delta.pdf')
    
    #plt.show()