"""
===========
F16 Model
===========



"""
import F16_aero
import F16_prop
import F16_inertia
import numpy as np
import simupy_flight
from nesc_testcase_helper import ft_per_m, kg_per_slug


N_per_lbf = 4.44822

S_A = F16_aero.sref/(ft_per_m**2)
b_l = F16_aero.bspan/ft_per_m
c_l = F16_aero.cbar/ft_per_m
a_l = b_l

class F16(simupy_flight.Vehicle):
    def __init__(self, CG_PCT_MAC=25):
        inertia_output = F16_inertia.F16_inertia(CG_PCT_MAC)
        Ixx, Iyy, Izz, Izx, Ixy, Iyz = inertia_output[:6]*kg_per_slug/(ft_per_m**2)
        #m = inertia_output[6]*kg_per_slug
        m = 637.26*kg_per_slug #slug
        y_com, z_com, x_com = inertia_output[7:]/ft_per_m
        super().__init__(
            m=m,
            I_xx=Ixx, I_yy=Iyy, I_zz=Izz,
            I_xy=Ixy, I_yz=Iyz, I_xz=Izx,
            x_com=x_com, y_com=y_com, z_com=z_com,

            base_aero_coeffs=0.,
            x_mrc=0., y_mrc=0., z_mrc=0.,
            S_A=S_A, a_l=a_l, b_l=b_l, c_l=c_l, d_l=0.,

            input_force_moment = self.prop_model,
            dim_additional_input=4,
        )

    def tot_aero_forces_moments(self, qbar, Ma, Re, V_T, alpha, beta, p_B, q_B, r_B, el, ail, rdr, *args):
        #print("hello")
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

    def prop_model(self, *args):
        throttle = args[-1]
        alt = args[simupy_flight.Vehicle.h_D_arg_idx]*ft_per_m
        Ma = args[simupy_flight.Vehicle.Ma_arg_idx]
        return F16_prop.F16_prop(throttle, alt, Ma)*N_per_lbf


