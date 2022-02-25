"""
===================================================
Case 11: Subsonic F-16 trimmed flight across planet
===================================================

==============  ===============
Verifies        Atmosphere, air-data calculations
Gravitation     J2
Geodesy         WGS-84 rotating
Atmosphere      US 1976 STD
Winds           still air
Vehicle         F-16 (unaugmented)
Notes           Initial position is 10,000 ft above KFFA airport on a 45 degree true
                course. 335.15 KTAS. Stability augmentation off. Test of trim solution.
==============  ===============

"""

from simupy.block_diagram import BlockDiagram
from simupy import systems
import simupy_flight
import numpy as np
from scipy import optimize

from nesc_testcase_helper import plot_nesc_comparisons, int_opts, benchmark
from nesc_testcase_helper import ft_per_m, kg_per_slug

import F16_model
import F16_model_01
from F16_control import F16_control

F16_vehicle = F16_model.F16_vehicle
F16_vehicle = F16_model_01.F16()

spec_ic_args = dict(
    phi_D = 36.01916667*np.pi/180,  # latitude
    lamda_D = -75.67444444*np.pi/180, # longitude
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

knots_per_mps = 1.94384

planet = simupy_flight.Planet(
    gravity=simupy_flight.earth_J2_gravity,
    winds=simupy_flight.get_constant_winds(),
    atmosphere=simupy_flight.atmosphere_1976,
    planetodetics=simupy_flight.Planetodetic(
        a=simupy_flight.earth_equitorial_radius,
        omega_p=simupy_flight.earth_rotation_rate,
        f=simupy_flight.earth_f
    )
)

rho_0 = planet.atmosphere(0, 0, 0, 0)[0]

controller_feedback_indices = np.array([
    planet.h_D_idx, planet.V_T_idx, planet.alpha_idx, planet.beta_idx,
    planet.psi_idx, planet.theta_idx, planet.phi_idx,
    planet.p_B_idx, planet.q_B_idx, planet.r_B_idx,
    planet.rho_idx])

dim_feedback = len(controller_feedback_indices)


def get_controller_function(throttleTrim, longStkTrim, sasOn=False, apOn=False):
    def controller_function(t, u):
        throttle, longStk, latStk, pedal = 0., 0., 0., 0. # pilot command

        (alt, V_T, alpha, beta, psi, theta, phi, pb, qb, rb, # feedback
            rho, # rho to calculate equivalent airspeed
            keasCmd, altCmd, latOffset, baseChiCmd) = u # commands

        Vequiv = V_T * np.sqrt(rho/rho_0)
        angles = np.array([alpha, beta, phi, theta, psi])
        alpha, beta, phi, theta, psi = angles*180/np.pi

        control_eart = F16_control(throttle, longStk, latStk, pedal,
            sasOn, apOn,
            keasCmd, altCmd, latOffset, baseChiCmd,
            alt*ft_per_m, Vequiv*knots_per_mps,
            alpha, beta, phi, theta, psi, pb, qb, rb, throttleTrim, longStkTrim)
        return control_eart
    return controller_function


def eval_trim(flight_condition, longStk, throttle):
    kin_out = planet.output_equation_function(0, flight_condition)
    controller_func = get_controller_function(throttleTrim=throttle, longStkTrim=longStk)
    #print("Eval...", F16_vehicle, F16_vehicle.tot_aero_forces_moments, sep="\n")
    aero_plus_prop_acceleration = simupy_flight.dynamics.dynamics_output_function(F16_vehicle, 0, *kin_out, *controller_func(0, np.zeros(dim_feedback+4)))

    gen_accel = aero_plus_prop_acceleration

    gen_accel[:3] = simupy_flight.kinematics.local_translational_trim_residual(planet, *flight_condition[:-3], *aero_plus_prop_acceleration[:-3]).squeeze()

    return gen_accel


def run_trimmer(flight_ic_args, throttle_ic=0., longStk_ic=0., allow_roll=False):
    len_vars = 3 + allow_roll

    psi, theta_ic, phi_ic = flight_ic_args['psi'], flight_ic_args['theta'], flight_ic_args['phi']

    initial_guess = np.zeros(len_vars)
    initial_guess[0] = theta_ic
    initial_guess[1] = throttle_ic
    initial_guess[2] = longStk_ic

    extra_index = 3

    if allow_roll:
        initial_guess[extra_index] = phi_ic
        extra_index += 1

    def parse_x(x):
        theta, throttle, longStk = x[:3]
        extra_index = 3

        if allow_roll:
            phi = x[extra_index]
            extra_index += 1
        else:
            phi = phi_ic

        return theta, phi, longStk, throttle

    weighting_matrix = np.eye(6)

    aileron, rudder = 0.0, 0.0

    def trim_opt_func(x):
        eval_args = flight_ic_args.copy()
        theta, phi, longStk, throttle = parse_x(x)
        eval_args['theta'] = theta
        eval_args['phi'] = phi
        flight_condition = planet.ic_from_planetodetic(**eval_args)
        return np.linalg.norm(weighting_matrix@eval_trim(flight_condition, longStk, throttle), ord=2)

    opt_res = optimize.minimize(trim_opt_func, initial_guess, tol=1E-12, options={'disp': True, 'adaptive': True, 'fatol': 1E-12, 'maxiter': 20_000, 'xatol': 1E-12}, method='Nelder-Mead')

    opt_theta, opt_phi, opt_longStk, opt_throttle = opt_result = parse_x(opt_res.x)
    opt_args = flight_ic_args.copy()
    opt_args['theta'] = opt_theta
    opt_args['phi'] = opt_phi
    opt_flight_condition = planet.ic_from_planetodetic(**opt_args)
    print("pitch: %.4e  roll: %.4e  longStk: %.4f  throttle: %.4f" % (opt_theta*180/np.pi, opt_phi*180/np.pi, opt_longStk*100, opt_throttle*100))
    print("accelerations:\n", eval_trim(opt_flight_condition, opt_longStk, opt_throttle).reshape((2,3)) )
    return opt_args, np.array([opt_throttle, opt_longStk])


opt_args, opt_ctrl = run_trimmer(spec_ic_args, throttle_ic=0.0, longStk_ic=0.0, allow_roll=False)

trimmed_flight_condition = planet.ic_from_planetodetic(**opt_args)

trimmed_KEAS = planet.output_equation_function(0, trimmed_flight_condition)[planet.V_T_idx]*np.sqrt(planet.output_equation_function(0, trimmed_flight_condition)[planet.rho_idx]/rho_0) * knots_per_mps
int_opts['nsteps'] = 5_000

planet.initial_condition = trimmed_flight_condition

controller_block = systems.SystemFromCallable(get_controller_function(*opt_ctrl), dim_feedback + 4, 4)

keasCmdOutput = np.array([trimmed_KEAS])
keasCmdBlock = systems.SystemFromCallable(lambda *args: keasCmdOutput, 0, 1)

altCmdOutput = np.array([spec_ic_args['h']*ft_per_m])
altCmdBlock = systems.SystemFromCallable(lambda *args: altCmdOutput, 0, 1)

def latOffsetStateEquation(t, x, u):
    chi_cmd, V_N, V_E = u
    V_ground_magnitude = np.sqrt(V_N**2 + V_E**2)
    V_ground_heading = np.arctan2(V_E, V_N)
    return V_ground_magnitude*np.sin( V_ground_heading - chi_cmd*np.pi/180)

def latOffsetOutputEquation(t, x):
    return x*ft_per_m

latOffsetBlock = systems.DynamicalSystem(state_equation_function=latOffsetStateEquation, output_equation_function=latOffsetOutputEquation, dim_state=1, dim_input=3, dim_output=1)

baseChiCmdOutput = np.array([spec_ic_args['psi']*180/np.pi])
baseChiCmdBlock = systems.SystemFromCallable(lambda *args: baseChiCmdOutput, 0, 1)

BD = BlockDiagram(planet, F16_vehicle, controller_block, keasCmdBlock, altCmdBlock, latOffsetBlock, baseChiCmdBlock)
BD.connect(planet, F16_vehicle, inputs=np.arange(planet.dim_output))
BD.connect(F16_vehicle, planet, inputs=np.arange(F16_vehicle.dim_output))
BD.connect(controller_block, F16_vehicle, inputs=np.arange(planet.dim_output, planet.dim_output+4))
BD.connect(planet, controller_block, outputs=controller_feedback_indices, inputs=np.arange(dim_feedback))

BD.connect(keasCmdBlock, controller_block, inputs=[dim_feedback+0])
BD.connect(altCmdBlock, controller_block, inputs=[dim_feedback+1])
BD.connect(latOffsetBlock, controller_block, inputs=[dim_feedback+2])
BD.connect(baseChiCmdBlock, controller_block, inputs=[dim_feedback+3])

if __name__ == "__main__":
    with benchmark() as b:
        res = BD.simulate(180, integrator_options=int_opts)

    plot_nesc_comparisons(res, '11')
