from simupy.block_diagram import BlockDiagram
from simupy import systems
import simupy_flight
import pandas as pd
import numpy as np
import os, sys
import glob
from scipy import optimize
from nesc_testcase_helper import plot_nesc_comparisons, nesc_options, int_opts, ft_per_m, kg_per_slug

# from F16_model import F16_vehicle, N_per_lbf
from F16_model_00 import F16_vehicle, N_per_lbf
from F16_control import F16_control, trimmedKEAS

planet = simupy_flight.Planet(
    gravity=simupy_flight.earth_J2_gravity,
    winds=simupy_flight.get_constant_winds(),
    atmosphere=simupy_flight.atmosphere_1976,
    planetodetics=simupy_flight.Planetodetic(a=simupy_flight.earth_equitorial_radius, omega_p=simupy_flight.earth_rotation_rate, f=simupy_flight.earth_f)
)


def get_ic_args_from_spec():
    return dict(
        phi_E = 36.01916667*np.pi/180,  # latitude
        lamda_E = -75.67444444*np.pi/180, # longitude
        h = 10_013/ft_per_m,
        
        V_N = 400./ft_per_m,
        V_E = 400./ft_per_m,
        V_D = 0./ft_per_m,
        
        psi = 45.0*np.pi/180,
        
        theta = 2.653814*np.pi/180,
        phi = 0.0*np.pi/180,
        
        p_B = 0.*np.pi/180,
        q_B = 0.*np.pi/180,
        r_B = 0.*np.pi/180
    )

spec_args = get_ic_args_from_spec()
ang_vel_1 = planet.ic_from_planetodetic(**spec_args)[-3:]
spec_args['V_N'] = 0.
spec_args['V_E'] = 0.
ang_vel_2 = planet.ic_from_planetodetic(**spec_args)[-3:]
spec_args['h'] = 0.
ang_vel_2 = planet.ic_from_planetodetic(**spec_args)[-3:]

data_relative_path = nesc_options['data_relative_path']
glob_path = os.path.join(data_relative_path, 'Atmospheric_checkcases', 'Atmos_11_TrimCheckSubsonicF16', 'Atmos_11_sim_*.csv')



def get_ic_args_from_baseline(idx):
    baseline_df = pd.read_csv(glob_path.replace('*','%02d' % idx ), index_col=0).iloc[0]
    
    try:
        long, lat, h = planet.planetodetics.pcf2pd(*baseline_df[['eiPosition_ft_X', 'eiPosition_ft_Y', 'eiPosition_ft_Z']]/ft_per_m)
    except KeyError:
        long, lat, h = planet.planetodetics.pcf2pd(*get_ic_from_spec()[:3])
        print("Missing inertial position...")
    psi, theta, phi = baseline_df[['eulerAngle_deg_Yaw', 'eulerAngle_deg_Pitch', 'eulerAngle_deg_Roll']]*np.pi/180
    flight_condition = planet.ic_from_planetodetic(
        lamda_E=long, phi_E=lat, h=h, 
        psi=psi, theta=theta, phi=phi,
    )
    try:
        flight_condition[:3] = baseline_df[['eiPosition_ft_X', 'eiPosition_ft_Y', 'eiPosition_ft_Z']]/ft_per_m
    except KeyError:
        pass
    try:
        flight_condition[7:10] = baseline_df[['eiVelocity_ft_s_X', 'eiVelocity_ft_s_Y', 'eiVelocity_ft_s_Z']]/ft_per_m
    except KeyError:
        print("Missing inertial velocity...")
    flight_condition[10:] = baseline_df[['bodyAngularRateWrtEi_deg_s_Roll', 'bodyAngularRateWrtEi_deg_s_Pitch', 'bodyAngularRateWrtEi_deg_s_Yaw']]*np.pi/180
    
    baseline_output = planet.output_equation_function(0, flight_condition)
    
    return dict(
        phi_E = lat,  # latitude
        lamda_E = long, # longitude
        h = h,
        
        V_N = baseline_output[planet.V_N_idx],
        V_E = baseline_output[planet.V_E_idx],
        V_D = baseline_output[planet.V_D_idx],
        
        psi = psi,
        
        theta = theta,
        phi = phi,
        
        p_B = baseline_output[planet.p_B_idx],
        q_B = baseline_output[planet.q_B_idx],
        r_B = baseline_output[planet.r_B_idx]
    )
    
rho_0 = planet.atmosphere(0, 0, 0, 0)[0]
knots_per_mps = 1.94384

def get_controller_function(throttleTrim, longStkTrim, 
throttle=0., longStk=0., latStk=0., pedal=0., 
keasCmd=trimmedKEAS, altCmd=10_013, latOffset=0.0, baseChiCmd=45.0, 
sasOn=False, apOn=False):
    
    def controller_function(t, u):
        # throttle, longStk, latStk, pedal = 0., 0., 0., 0. # pilot command
        (alt, V_T, alpha, beta, psi, theta, phi, pb, qb, rb, # feedback
            rho) = u # rho to calculate equivalent airspeed
    
        Vequiv = V_T * np.sqrt(rho/rho_0)
        angles = np.array([alpha, beta, phi, theta, psi, pb, qb, rb])
        alpha, beta, phi, theta, psi, pb, qb, rb = angles*180/np.pi
        return F16_control(throttle, longStk, latStk, pedal,
            sasOn, apOn, 
            keasCmd, altCmd, latOffset, baseChiCmd, 
            alt*ft_per_m, Vequiv*knots_per_mps, 
            alpha, beta, phi, theta, psi, pb, qb, rb, throttleTrim, longStkTrim)
    return controller_function


def eval_trim(flight_condition, elevator, aileron, rudder, throttle):
    kin_out = planet.output_equation_function(0, flight_condition)

    controller_func = get_controller_function(throttleTrim=throttle/100, longStkTrim=elevator/-25., throttle=0., longStk=0., latStk=aileron, pedal=rudder)
    aero_plus_prop_acceleration = simupy_flight.dynamics.dynamics_output_function(F16_vehicle, 0, *kin_out, *controller_func(0, np.zeros(11)))
    
    gen_accel = aero_plus_prop_acceleration
    
    gen_accel[:3] = simupy_flight.kinematics.local_translational_trim_residual(planet, *flight_condition[:-3], *aero_plus_prop_acceleration[:-3]).squeeze()

    return gen_accel


def run_trimmer(flight_ic_args, throttle_ic=13.9, elevator_ic=-3.241, allow_roll=False, rudder_ic=None, aileron_ic=None):
    len_vars = 3 + allow_roll + (rudder_ic is not None) + (aileron_ic is not None)
    
    psi, theta_ic, phi_ic = flight_ic_args['psi'], flight_ic_args['theta'], flight_ic_args['phi']
    
    initial_guess = np.zeros(len_vars)
    initial_guess[0] = theta_ic
    initial_guess[1] = throttle_ic
    initial_guess[2] = elevator_ic
    
    extra_index = 3
    
    if allow_roll:
        initial_guess[extra_index] = phi_ic
        extra_index += 1        

    if rudder_ic is not None:
        initial_guess[extra_index] = rudder_ic
        extra_index += 1
    
    if aileron_ic is not None:
        initial_guess[extra_index] = aileron_ic
        extra_index += 1
    
    def parse_x(x):
        theta, throttle, elevator = x[:3]
        extra_index = 3
        
        if allow_roll:
            phi = x[extra_index]
            extra_index += 1
        else:
            phi = phi_ic
    
        if rudder_ic is not None:
            rudder = x[extra_index]
            extra_index += 1
        else:
            rudder = 0.0
        
        if aileron_ic is not None:
            aileron = x[extra_index]
            extra_index += 1
        else:
            aileron = 0.0
        
        return theta, phi, elevator, aileron, rudder, throttle
    
    weighting_matrix = np.eye(6)
    # weighting_matrix = np.diag([10, 10, 10, 0, 0, 0])
    # weighting_matrix[3,3] = 10
    # weighting_matrix[4,4] = 10
    # weighting_matrix[5,5] = 10
    
    
    def trim_opt_func(x):
        eval_args = flight_ic_args.copy()
        theta, phi, elevator, aileron, rudder, throttle = parse_x(x)
        eval_args['theta'] = theta
        eval_args['phi'] = phi
        flight_condition = planet.ic_from_planetodetic(**eval_args)
        return np.linalg.norm(weighting_matrix@eval_trim(flight_condition, elevator, aileron, rudder, throttle), ord=2)

    opt_res = optimize.minimize(trim_opt_func, initial_guess, tol=1E-12, options={'disp': True, 'adaptive': True, 'fatol': 1E-12, 'maxiter': 20_000, 'xatol': 1E-12}, method='Nelder-Mead')
    # opt_res = optimize.minimize(trim_opt_func, initial_guess, tol=1E-12, options={'disp': True, 'ftol': 1E-12,}, method='SLSQP')
    
    opt_theta, opt_phi, opt_elevator, opt_aileron, opt_rudder, opt_throttle = opt_result = parse_x(opt_res.x)
    opt_args = flight_ic_args.copy()
    opt_args['theta'] = opt_theta
    opt_args['phi'] = opt_phi
    opt_flight_condition = planet.ic_from_planetodetic(**opt_args)
    print("pitch: %.4e  roll: %.4e  elevator: %.4f  aileron: %.4f  rudder: %.4f  throttle: %.4f" % (opt_theta*180/np.pi, opt_phi*180/np.pi, opt_elevator, opt_aileron, opt_rudder, opt_throttle))
    print("accelerations:\n", eval_trim(opt_flight_condition, opt_elevator, opt_aileron, opt_rudder, opt_throttle).reshape((2,3)) )
    return opt_args, np.array([opt_elevator, opt_aileron, opt_rudder, opt_throttle])
    
# opt_args, opt_ctrl = run_trimmer(get_ic_args_from_spec(), throttle_ic=13.9, elevator_ic=-3.241, allow_roll=True, rudder_ic=0., aileron_ic=0.0)
opt_args, opt_ctrl = run_trimmer(get_ic_args_from_spec(), throttle_ic=13.9, elevator_ic=-3.241, allow_roll=False, rudder_ic=None, aileron_ic=None)
# opt_args, opt_ctrl = run_trimmer(get_ic_args_from_spec(), throttle_ic=13.9, elevator_ic=-3.241, allow_roll=True, rudder_ic=None, aileron_ic=None)

# run_trimmer(get_ic_from_spec(), throttle_ic=13.9, elevator_ic=-3.241, allow_roll=False, rudder_ic=None, aileron_ic=None)
# run_trimmer(get_ic_from_baseline(5), throttle_ic=13.9, elevator_ic=-3.241, allow_roll=True, rudder_ic=0., aileron_ic=0.0)

# opt_args, opt_ctrl = run_trimmer(get_ic_args_from_spec(), throttle_ic=13.9, elevator_ic=-3.241, allow_roll=False, rudder_ic=None, aileron_ic=None)

# run_trimmer(get_ic_args_from_spec(), throttle_ic=13.9, elevator_ic=-3.241, allow_roll=True, rudder_ic=None, aileron_ic=None)
# run_trimmer(get_ic_from_baseline(5), throttle_ic=13.9, elevator_ic=-3.241, allow_roll=True, rudder_ic=0., aileron_ic=0.0)
# run_trimmer(get_ic_from_baseline(5), throttle_ic=13.9, elevator_ic=-3.241, allow_roll=True, rudder_ic=None, aileron_ic=None)
# ESD_cond = update_flight_condition(get_ic_from_spec(), psi=0., theta=0., phi=0.,)
##


int_opts['nsteps'] = 5_000
# int_opts['max_step'] = 2**-5

flight_condition = planet.ic_from_planetodetic(**opt_args)

planet.initial_condition = flight_condition

trim_ctrl_func = lambda t: opt_ctrl
trim_ctrl_block = systems.SystemFromCallable(trim_ctrl_func, 0, 4)

BD = BlockDiagram(planet, F16_vehicle, trim_ctrl_block)
BD.connect(planet, F16_vehicle, inputs=np.arange(planet.dim_output))
BD.connect(F16_vehicle, planet, inputs=np.arange(F16_vehicle.dim_output))
BD.connect(trim_ctrl_block, F16_vehicle, inputs=np.arange(planet.dim_output, planet.dim_output+4))

import cProfile
cProfile.run('res = BD.simulate(180, integrator_options=int_opts)')

import time
tstart = time.time()
res = BD.simulate(180, integrator_options=int_opts)
tend = time.time()
tdelta = tend - tstart
print("time to simulate: %f    eval time to run time: %f" % (tdelta, res.t[-1]/tdelta))

plot_nesc_comparisons(res, glob_path, '11')

