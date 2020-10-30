from simupy.block_diagram import BlockDiagram
import simupy_flight
import pandas as pd, matplotlib.pyplot as plt
import numpy as np
import os
import glob
from nesc_testcase_helper import plot_nesc_comparisons, data_relative_path, int_opts, ft_per_m, kg_per_slug

planet = simupy_flight.Planet(
    gravity=simupy_flight.earth_J2_gravity,
    winds=simupy_flight.get_constant_winds(),
    density=simupy_flight.density_1976_atmosphere,
    speed_of_sound=simupy_flight.get_constant_speed_of_sound(),
    viscosity=simupy_flight.get_constant_viscosity(),
    planetodetics=simupy_flight.Planetodetic(a=simupy_flight.earth_equitorial_radius, omega_p=simupy_flight.earth_rotation_rate, f=simupy_flight.earth_f)
)

Ixx = 0.001894220*kg_per_slug/(ft_per_m**2) #slug-ft2
Iyy = 0.006211019*kg_per_slug/(ft_per_m**2) #slug-ft2
Izz = 0.007194665*kg_per_slug/(ft_per_m**2) #slug-ft2
Ixy = 0.0*kg_per_slug/(ft_per_m**2) #slug-ft2
Iyz = 0.0*kg_per_slug/(ft_per_m**2) #slug-ft2
Izx = 0.0*kg_per_slug/(ft_per_m**2) #slug-ft2
m = 0.155404754*kg_per_slug #slug

x = 0.
y = 0.
z = 0.

S_A = 0.22222/(ft_per_m**2)
b_l = 1/(3*ft_per_m)
c_l = 2/(3*ft_per_m)
a_l = b_l
aero_model = simupy_flight.get_constant_aero(Cp_b=-1.0, Cq_b=-1.0, Cr_b=-1.0)
vehicle =  simupy_flight.Vehicle(base_aero_coeffs=aero_model, m=m, I_xx=Ixx, I_yy=Iyy, I_zz=Izz, I_xy=Ixy, I_yz=Iyz, I_xz=Izx, x_com=x, y_com=y, z_com=z, x_mrc=x, y_mrc=y, z_mrc=z, S_A=S_A, a_l=a_l, b_l=b_l, c_l=c_l, d_l=0.,)

BD = BlockDiagram(planet, vehicle)
BD.connect(planet, vehicle, inputs=np.arange(planet.dim_output))
BD.connect(vehicle, planet, inputs=np.arange(vehicle.dim_output))

lat_ic = 0.*np.pi/180
long_ic = 0.*np.pi/180
h_ic = 30_000/ft_per_m
V_N_ic = 0.
V_E_ic = 0.
V_D_ic = 0.
psi_ic = 0.*np.pi/180
theta_ic = 0.*np.pi/180
phi_ic = 0.*np.pi/180
omega_X_ic = 10.*np.pi/180
omega_Y_ic = 20.*np.pi/180
omega_Z_ic = 30.*np.pi/180

planet.initial_condition = planet.ic_from_planetodetic(long_ic, lat_ic, h_ic, V_N_ic, V_E_ic, V_D_ic, psi_ic, theta_ic, phi_ic)
planet.initial_condition[-3:] = omega_X_ic, omega_Y_ic, omega_Z_ic

res = BD.simulate(30, integrator_options=int_opts)

baseline_pds = []
for fname in glob.glob(os.path.join(data_relative_path, 'Atmospheric_checkcases', 'Atmos_03_TumblingBrickDamping', 'Atmos_03_sim_*.csv'),):
    baseline_pds.append(pd.read_csv(fname, index_col=0))

plot_nesc_comparisons(res, baseline_pds)