# run after process_NESC_DaveML
import F16_aero
import F16_prop
import F16_inertia
import numpy as np
import simupy_flight
from nesc_testcase_helper import ft_per_m, kg_per_slug, N_per_lbf

x_mrc = 0.
y_mrc = 0.
z_mrc = 0.

inertia_output = F16_inertia.F16_inertia(25)
Ixx, Iyy, Izz, Ixy, Iyz, Izx = inertia_output[:6]*kg_per_slug/(ft_per_m**2)
m = inertia_output[6]*kg_per_slug
x_com, y_com, z_com = inertia_output[[-1, -3, -2]]/ft_per_m

aero_output = F16_aero.F16_aero(*np.zeros(9))
S_A = aero_output[2]/(ft_per_m**2)
a_l, b_l, c_l = aero_output[[1, 1, 0]]/ft_per_m

def F16_prop_model(*args):
    """
    Wrapper to convert units appropriately
    """
    throttle = args[-1]
    alt = args[simupy_flight.Vehicle.h_D_arg_idx]*ft_per_m
    Ma = args[simupy_flight.Vehicle.Ma_arg_idx]
    return F16_prop.F16_prop(throttle, alt, Ma)*N_per_lbf

class VehicleWithAlternateAeroSturcture(simupy_flight.Vehicle):
    def tot_aero_forces_moments(self, qbar, Ma, Re, V_T, alpha, beta, p_B, q_B, r_B, el, ail, rdr, *args):
        """
        Alternate aero structure that depends only on base_aero_coeffs with three inputs
        (elevator, aileron, and rudder), forces in body-fixed coordinates instead of wind-fixed,
        uses impertial units, degrees for angles, and depends on body rates.
        """
        aero_out = self.base_aero_coeffs(
            V_T*ft_per_m, alpha*180/np.pi, beta*180/np.pi,
            p_B, q_B, r_B,
            el, ail, rdr
        )
        cx, cy, cz, cl, cm, cn = aero_out[-6:]
        forces_moments = np.empty(6)
        body_force_coeffs = np.array([cx, cy, cz])
        body_moment_coeffs = np.array([cl, cm, cn])
        dynamic_force = qbar * self.S_A
        ref_lengths = np.array([self.a_l, self.c_l, self.b_l])
        forces_moments[:3] = dynamic_force * body_force_coeffs
        forces_moments[3:] = (dynamic_force * ref_lengths * body_moment_coeffs) - (simupy_flight.dynamics.mrc_to_com_cpm(self) @ forces_moments[:3])
        return forces_moments

F16_vehicle = VehicleWithAlternateAeroSturcture(
    m=m,
    I_xx=Ixx, I_yy=Iyy, I_zz=Izz,
    I_xy=Ixy, I_yz=Iyz, I_xz=Izx,
    x_com=x_com, y_com=y_com, z_com=z_com,

    base_aero_coeffs=F16_aero.F16_aero,
    x_mrc=x_mrc, y_mrc=y_mrc, z_mrc=z_mrc,
    S_A=S_A, a_l=a_l, b_l=b_l, c_l=c_l, d_l=0.,

    input_force_moment = F16_prop_model,
    dim_additional_input=4,
)
