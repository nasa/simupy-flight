"""
=====================================================================
Case 10: Sphere launched ballistically northward along Prime Meridian
=====================================================================

==============  ===============
Verifies        Coriolis
Gravitation     J2
Geodesy         WGS-84 rotating
Atmosphere      US 1976 STD
Winds           still air
Vehicle         Sphere with constant :math:`C_D`
Notes           Initial velocity is :math:`\\sqrt{2000}` ft/s aligned 45 degrees from
                vertical, heading north; zero angular rate relative to launch platform.
==============  ===============
"""

from simupy.block_diagram import BlockDiagram
import simupy_flight
import numpy as np

from nesc_testcase_helper import plot_nesc_comparisons, int_opts, benchmark
from nesc_testcase_helper import ft_per_m, kg_per_slug

Ixx = 3.6 * kg_per_slug / (ft_per_m**2)  # slug-ft2
Iyy = 3.6 * kg_per_slug / (ft_per_m**2)  # slug-ft2
Izz = 3.6 * kg_per_slug / (ft_per_m**2)  # slug-ft2
Ixy = 0.0 * kg_per_slug / (ft_per_m**2)  # slug-ft2
Iyz = 0.0 * kg_per_slug / (ft_per_m**2)  # slug-ft2
Izx = 0.0 * kg_per_slug / (ft_per_m**2)  # slug-ft2
m = 1.0 * kg_per_slug  # slug

x = 0.0
y = 0.0
z = 0.0

S_A = 0.1963495 / (ft_per_m**2)
b_l = 1.0
c_l = 1.0
a_l = b_l

lat_ic = 0.0 * np.pi / 180
long_ic = 0.0 * np.pi / 180
h_ic = 0.0 / ft_per_m
V_N_ic = 1000.0 / ft_per_m
V_E_ic = 000.0 / ft_per_m
V_D_ic = -1000.0 / ft_per_m
psi_ic = 0.0 * np.pi / 180
theta_ic = 0.0 * np.pi / 180
phi_ic = 0.0 * np.pi / 180
p_b_ic = 0.0 * np.pi / 180
q_b_ic = 0.0 * np.pi / 180
r_b_ic = 0.0 * np.pi / 180
# omega_X_ic = 0.004178073*np.pi/180
# omega_Y_ic = 0.*np.pi/180
# omega_Z_ic = 0.*np.pi/180

planet = simupy_flight.Planet(
    gravity=simupy_flight.earth_J2_gravity,
    winds=simupy_flight.get_constant_winds(),
    atmosphere=simupy_flight.atmosphere_1976,
    planetodetics=simupy_flight.Planetodetic(
        a=simupy_flight.earth_equitorial_radius,
        omega_p=simupy_flight.earth_rotation_rate,
        f=simupy_flight.earth_f,
    ),
)

vehicle = simupy_flight.Vehicle(
    base_aero_coeffs=simupy_flight.get_constant_aero(CD_b=0.1),
    m=m,
    I_xx=Ixx,
    I_yy=Iyy,
    I_zz=Izz,
    I_xy=Ixy,
    I_yz=Iyz,
    I_xz=Izx,
    x_com=x,
    y_com=y,
    z_com=z,
    x_mrc=x,
    y_mrc=y,
    z_mrc=z,
    S_A=S_A,
    a_l=a_l,
    b_l=b_l,
    c_l=c_l,
    d_l=0.0,
)

BD = BlockDiagram(planet, vehicle)
BD.connect(planet, vehicle, inputs=np.arange(planet.dim_output))
BD.connect(vehicle, planet, inputs=np.arange(vehicle.dim_output))

planet.initial_condition = planet.ic_from_planetodetic(
    lamda_D=long_ic,
    phi_D=lat_ic,
    h=h_ic,
    V_N=V_N_ic,
    V_E=V_E_ic,
    V_D=V_D_ic,
    psi=psi_ic,
    theta=theta_ic,
    phi=phi_ic,
    p_B=p_b_ic,
    q_B=q_b_ic,
    r_B=r_b_ic,
)
# planet.initial_condition[-3:] = omega_X_ic, omega_Y_ic, omega_Z_ic
planet.initial_condition[-2] = 0.0

with benchmark() as b:
    res = BD.simulate(30, integrator_options=int_opts)

# %%

plot_nesc_comparisons(res, "10")
