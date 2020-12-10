# run after process_NESC_DaveML
import F16_aero
import F16_prop
import F16_inertia
import numpy as np
import simupy_flight
from nesc_testcase_helper import ft_per_m, kg_per_slug


def transformed_aero_coefficients(vt, alpha, beta, p, q, r, el, ail, rdr):
    cbar, bspan, sref, cx, cy, cz, cl, cm, cn = F16_aero.F16_aero(vt, alpha, beta, p, q, r, el, ail, rdr)
    cwx, cwy, cwz = simupy_flight.dynamics.body_to_wind_dcm(alpha, beta)@np.array([cx, cy, cz])
    return np.array(np.broadcast_array(-cwx, cwy, -cwz, cl, cm, cn, 0, 0, 0))

S_A = 300.0/(ft_per_m**2)
b_l = 30./ft_per_m
c_l = 11.32/ft_per_m
a_l = b_l

Ixx = 9_496.0*kg_per_slug/(ft_per_m**2) #slug-ft2
Iyy = 55_814.0*kg_per_slug/(ft_per_m**2) #slug-ft2
Izz = 63_100.0*kg_per_slug/(ft_per_m**2) #slug-ft2
Ixy = 0.0*kg_per_slug/(ft_per_m**2) #slug-ft2
Iyz = 0.0*kg_per_slug/(ft_per_m**2) #slug-ft2
Izx = 982.0*kg_per_slug/(ft_per_m**2) #slug-ft2
m = 637.26*kg_per_slug #slug

longitudinal_com_pos_percent = 25.
x_com = 0.01*c_l*(35.-longitudinal_com_pos_percent)
y_com = 0.
z_com = 0.

x_mrc = 0.
y_mrc = 0.
z_mrc = 0.

# Ixx, Iyy, Izz, Ixy, Iyz, Izx, m, x_com, y_com, z_com = F16_inertia.F16_inertia(0.25)



def F16_aero_model(self, qbar, Ma, Re, V_T, alpha, beta, p_B, q_B, r_B, el, ail, rdr, *args):
    # el, ail, rdr = args[0][:-1]
    cbar, bspan, sref, cx, cy, cz, cl, cm, cn = F16_aero.F16_aero(V_T*ft_per_m, alpha*180/np.pi, beta*180/np.pi, p_B, q_B, r_B, el, ail, rdr)
    forces_moments = np.empty(6)
    body_force_coeffs = np.array([cx, cy, cz])
    body_moment_coeffs = np.array([cl, cm, cn])
    dynamic_force = qbar * self.S_A
    ref_lengths = np.array([self.a_l, self.c_l, self.b_l])
    forces_moments[:3] = dynamic_force * body_force_coeffs
    forces_moments[3:] = (dynamic_force * ref_lengths * body_moment_coeffs) - (simupy_flight.dynamics.mrc_to_com_cpm(self) @ forces_moments[:3])
    return forces_moments

N_per_lbf = 4.44822
def F16_prop_model(*args):
    throttle = args[-1]
    alt = args[simupy_flight.Vehicle.h_D_arg_idx]*ft_per_m
    Ma = args[simupy_flight.Vehicle.Ma_arg_idx]
    return F16_prop.F16_prop(throttle, alt, Ma)*N_per_lbf
    

aero_model = simupy_flight.get_constant_aero()
F16_vehicle = simupy_flight.Vehicle(
    
    m=m, 
    I_xx=Ixx, I_yy=Iyy, I_zz=Izz, 
    I_xy=Ixy, I_yz=Iyz, I_xz=Izx, 
    x_com=x_com, y_com=y_com, z_com=z_com, 
    
    base_aero_coeffs=aero_model, 
    x_mrc=x_mrc, y_mrc=y_mrc, z_mrc=z_mrc, 
    S_A=S_A, a_l=a_l, b_l=b_l, c_l=c_l, d_l=0.,

    input_force_moment = F16_prop_model,
    dim_additional_input=4,
)

setattr(F16_vehicle, 
    F16_vehicle.tot_aero_forces_moments.__name__, 
    F16_aero_model.__get__(F16_vehicle, F16_vehicle.__class__))


