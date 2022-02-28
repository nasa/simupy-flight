import numpy


def tot_aero_forces_moments(self, qbar, Ma, Re, V_T, alpha, beta, p_B, q_B, r_B, *args):
    [CD_b, CS_b, CL_b, CLcal_b, CMcal_b, CNcal_b, Cp_b, Cq_b, Cr_b] = self.base_aero_coeffs(alpha,beta,Ma,Re).ravel()
    [CD_e, CS_e, CL_e, CLcal_e, CMcal_e, CNcal_e, Cp_e, Cq_e, Cr_e] = self.input_aero_coeffs(alpha,beta,Ma,Re,*args).ravel()
    x0 = -CL_b - CL_e
    x1 = numpy.sin(alpha)
    x2 = numpy.cos(alpha)
    x3 = numpy.cos(beta)
    x4 = -CD_b - CD_e
    x5 = x3*x4
    x6 = CS_b + CS_e
    x7 = numpy.sin(beta)
    x8 = x6*x7
    x9 = -x0*x1 + x2*x5 - x2*x8
    x10 = self.S_A*qbar
    x11 = x3*x6 + x4*x7
    x12 = x0*x2 + x1*x5 - x1*x8
    x13 = self.z_com - self.z_mrc
    x14 = self.y_com - self.y_mrc
    x15 = numpy.select([numpy.greater(V_T, 0.0)], [(1/2)/V_T], default=0.0)
    x16 = self.x_com - self.x_mrc
    return (numpy.array([x10*x9, x10*x11, x10*x12, self.a_l*x10*(CLcal_b + CLcal_e + self.a_l*p_B*x15*(Cp_b + Cp_e) - (x11*x13 - x12*x14)/self.a_l), self.c_l*x10*(CMcal_b + CMcal_e + self.c_l*q_B*x15*(Cq_b + Cq_e) - (x12*x16 - x13*x9)/self.c_l), self.b_l*x10*(CNcal_b + CNcal_e + self.b_l*r_B*x15*(Cr_b + Cr_e) - (-x11*x16 + x14*x9)/self.b_l)]))

def dynamics_output_function(self, t, p_x, p_y, p_z, q_0, q_1, q_2, q_3, v_x, v_y, v_z, omega_X, omega_Y, omega_Z, lamda_D, phi_D, h_D, psi, theta, phi, rho, c_s, mu, V_T, alpha, beta, p_B, q_B, r_B, V_N, V_E, V_D, W_N, W_E, W_D, *args):
    qbar = (1/2)*V_T**2*rho
    Ma = V_T/c_s
    Re = V_T*self.d_l*rho/mu
    [F_ax, F_ay, F_az, Lcal, Mcal, Ncal] = self._tot_aero_forces_moments(qbar,Ma,Re,V_T,alpha,beta,p_B,q_B,r_B,*args).ravel()
    [Phi_x, Phi_y, Phi_z, tau_x, tau_y, tau_z] = self._input_force_moment(t,p_x,p_y,p_z,q_0,q_1,q_2,q_3,v_x,v_y,v_z,omega_X,omega_Y,omega_Z,lamda_D,phi_D,h_D,psi,theta,phi,rho,c_s,mu,V_T,alpha,beta,p_B,q_B,r_B,V_N,V_E,V_D,W_N,W_E,W_D,qbar,Ma,Re,*args).ravel()
    [F_x, F_y, F_z, M_x, M_y, M_z] = numpy.array([F_ax + Phi_x, F_ay + Phi_y, F_az + Phi_z, Lcal + tau_x, Mcal + tau_y, Ncal + tau_z]).ravel()
    x0 = 1/self.m
    x1 = self.I_yz**2
    x2 = self.I_xz**2
    x3 = self.I_xy**2
    x4 = self.I_yy*self.I_zz
    x5 = self.I_xy*self.I_yz
    x6 = 1/(self.I_xx*x1 - self.I_xx*x4 + 2*self.I_xz*x5 + self.I_yy*x2 + self.I_zz*x3)
    x7 = omega_Y**2
    x8 = omega_Z**2
    x9 = omega_X*omega_Y
    x10 = omega_Y*omega_Z
    x11 = omega_X*omega_Z
    x12 = -self.I_xy*x11 + self.I_xz*x9 + self.I_yy*x10 + self.I_yz*x7 - self.I_yz*x8 - self.I_zz*x10 + M_x
    x13 = self.I_xz*self.I_yy + x5
    x14 = omega_X**2
    x15 = self.I_xx*x9 + self.I_xy*x14 - self.I_xy*x7 - self.I_xz*x10 - self.I_yy*x9 + self.I_yz*x11 + M_z
    x16 = self.I_xy*self.I_zz + self.I_xz*self.I_yz
    x17 = -self.I_xx*x11 + self.I_xy*x10 - self.I_xz*x14 + self.I_xz*x8 - self.I_yz*x9 + self.I_zz*x11 + M_y
    x18 = self.I_xx*self.I_yz + self.I_xy*self.I_xz
    return (numpy.array([F_x*x0, F_y*x0, F_z*x0, -x6*(x12*(-x1 + x4) + x13*x15 + x16*x17), -x6*(x12*x16 + x15*x18 + x17*(self.I_xx*self.I_zz - x2)), -x6*(x12*x13 + x15*(self.I_xx*self.I_yy - x3) + x17*x18)]))

def mrc_to_com_cpm(self):
    return (numpy.array([[0, -self.z_com + self.z_mrc, self.y_com - self.y_mrc], [self.z_com - self.z_mrc, 0, -self.x_com + self.x_mrc], [-self.y_com + self.y_mrc, self.x_com - self.x_mrc, 0]]))

def body_to_wind_dcm(alpha, beta):
    x0 = numpy.cos(alpha)
    x1 = numpy.cos(beta)
    x2 = numpy.sin(beta)
    x3 = numpy.sin(alpha)
    return (numpy.array([[x0*x1, x2, x1*x3], [-x0*x2, x1, -x2*x3], [-x3, 0, x0]]))
