from simupy.block_diagram import BlockDiagram
import simupy_flight
import numpy as np
import os
import glob
from nesc_testcase_helper import plot_nesc_comparisons, data_relative_path, int_opts, ft_per_m

planet = simupy_flight.Planet(
    gravity=simupy_flight.earth_J2_gravity,
    winds=simupy_flight.get_constant_winds(),
    atmosphere=simupy_flight.get_constant_atmosphere(),
    planetodetics = simupy_flight.Planetodetic(a=simupy_flight.earth_equitorial_radius, omega_p=simupy_flight.earth_rotation_rate, f=simupy_flight.earth_f),
)

BD = BlockDiagram(planet)

lat_ic = 0.*np.pi/180
long_ic = 0.*np.pi/180
h_ic = 30_000/ft_per_m
V_N_ic = 0.
V_E_ic = 0.
V_D_ic = 0.
psi_ic = 0.*np.pi/180
theta_ic = 0.*np.pi/180
phi_ic = 0.*np.pi/180
omega_X_ic = 0.
omega_Y_ic = 0.
omega_Z_ic = 0.

planet.initial_condition = planet.ic_from_planetodetic(long_ic, lat_ic, h_ic, V_N_ic, V_E_ic, V_D_ic, psi_ic, theta_ic, phi_ic)
planet.initial_condition[-3:] = omega_X_ic, omega_Y_ic, omega_Z_ic


import time
tstart = time.time()
res = BD.simulate(30, integrator_options=int_opts)
tend = time.time()
tdelta = tend - tstart
print("time to simulate: %f    eval time to run time: %f" % (tdelta, res.t[-1]/tdelta))

glob_path = os.path.join(data_relative_path, 'Atmospheric_checkcases', 'Atmos_01_DroppedSphere', 'Atmos_01_sim_*.csv')
plot_nesc_comparisons(res, glob_path, '01')
