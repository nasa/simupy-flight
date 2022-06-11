"""
===================================================
Case 1: Dropped sphere with no drag
===================================================

==============  ===============
Verifies        Gravitation, translational EOM
Gravitation     J2
Geodesy         WGS-84 rotating
Atmosphere      US 1976 STD
Winds           still air
Vehicle         Dragless sphere
Notes           Drag coefficient set to zero
==============  ===============
"""

from simupy.block_diagram import BlockDiagram
import simupy_flight
import numpy as np

from nesc_testcase_helper import plot_nesc_comparisons, int_opts, benchmark
from nesc_testcase_helper import ft_per_m

# %%
# Construct the Planet model according to the
# configuration listed above

planet = simupy_flight.Planet(
    gravity=simupy_flight.earth_J2_gravity,
    winds=simupy_flight.get_constant_winds(),
    atmosphere=simupy_flight.get_constant_atmosphere(),
    planetodetics=simupy_flight.Planetodetic(
        a=simupy_flight.earth_equitorial_radius,
        omega_p=simupy_flight.earth_rotation_rate,
        f=simupy_flight.earth_f,
    ),
)

# %%
# Since there is no dynamics model (aerodynamic drag is set to zero)
# there is no need to define a Vehicle model. So we construct the
# BlockDiagram with just the planet object

BD = BlockDiagram(planet)

# %%
# List the planetodetic initial conditions and use the
# `ic_from_planetodetic` method to transform it to the inertial
# coordinate system needed for the initial condition

lat_ic = 0.0 * np.pi / 180
long_ic = 0.0 * np.pi / 180
h_ic = 30_000 / ft_per_m
V_N_ic = 0.0
V_E_ic = 0.0
V_D_ic = 0.0
psi_ic = 0.0 * np.pi / 180
theta_ic = 0.0 * np.pi / 180
phi_ic = 0.0 * np.pi / 180
omega_X_ic = 0.0
omega_Y_ic = 0.0
omega_Z_ic = 0.0

planet.initial_condition = planet.ic_from_planetodetic(
    long_ic, lat_ic, h_ic, V_N_ic, V_E_ic, V_D_ic, psi_ic, theta_ic, phi_ic
)
planet.initial_condition[-3:] = omega_X_ic, omega_Y_ic, omega_Z_ic

# %%
# Simulate the dropped sphere
#
# .. note::
#
#    The generated code performs a divide by zero if the velocity is zero, generating
#    ``RuntimeWarning``\s. However, this condition is checked and handled correctly

with benchmark() as b:
    res = BD.simulate(30, integrator_options=int_opts)

# %%
# Plot the results of the simulation

plot_nesc_comparisons(res, "01")
