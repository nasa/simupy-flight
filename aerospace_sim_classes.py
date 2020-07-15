import numpy
import numpy as np
import ndsplines
import astropy._erfa
import fluids.atmosphere

earth_spherical_gravity_constant = 14076443110000000/(3.28084**3)
earth_f = 1/298.257223563
earth_J2_constant = 0.00108262982
earth_rotation_rate_deg = 0.004178073
earth_rotation_rate = earth_rotation_rate_deg*np.pi/180
earth_equitorial_radius = 6378137.0

# kinematic functions
def get_nonflat_gc2gd(a, f):
    def gc2gde(px, py, pz):
        return np.array(astropy._erfa.gc2gde(a, f, np.array([px, py, pz])))
    return gc2gde

def get_flat_gc2gd():
    def gc2gde(px, py, pz):
        return np.array([px, py, pz])
    return gc2gde

def get_nonflat_gd2gc(a, f):
    def gd2gce(longitude, latitude, altitude):
        return np.array(astropy._erfa.gd2gce(a, f, longitude, latitude, altitude))
    return gd2gce

def get_flat_gd2gc():
    def gd2gce(longitude, latitude, altitude):
        return np.array([longitude, latitude, altitude])
    return gd2gce

# constant atmosphere functions
def get_constant_viscosity(viscosity_val=0.):
    def viscosity_function(t,ap_x,ap_y,ap_z):
        return viscosity_val
    return viscosity_function

def get_constant_density(density_val=0.):
    def density_function(t,ap_x,ap_y,ap_z):
        return density_val
    return density_function

def get_constant_speed_of_sound(speed_of_sound_val=1.):
    def speed_of_sound_function(t,ap_x,ap_y,ap_z):
        return speed_of_sound_val
    return speed_of_sound_function

def density_1976_atmosphere(t,ap_x,ap_y,ap_z):
    atmo = fluids.atmosphere.ATMOSPHERE_1976(ap_z)
    return atmo.rho

# gravity functions
def get_spherical_gravity(gravitational_constant):
    def gravity_function(px, py, pz):
        pos_vec = np.array([px, py, pz])
        return -gravitational_constant*pos_vec/np.linalg.norm(pos_vec)**3
    return gravity_function

def earth_J2_gravity(px, py, pz):
    pos_vec = np.array([px, py, pz])
    r = np.linalg.norm(pos_vec)
    return -earth_spherical_gravity_constant*(pos_vec/r**3)*(1-3*earth_J2_constant*earth_equitorial_radius**2*(5*pz**2-r**2)/(2*r**4))

# wind functions
def get_constant_winds(wx=0., wy=0., wz=0.):
    wind_val = np.array([wx, wy, wz])
    def winds_function(t,ap_x,ap_y,ap_z):
        return wind_val
    return winds_function

# aero functions
def get_constant_aero(CD_b=0., CS_b=0., CL_b=0., CLcal_b=0., CMcal_b=0., CNcal_b=0., Cp_b=0., Cq_b=0., Cr_b=0.):
    aero_vals = np.array([CD_b, CS_b, CL_b, CLcal_b, CMcal_b, CNcal_b, Cp_b, Cq_b, Cr_b])
    def aero_function(alpha,beta,Ma,Re):
        return aero_vals
    return aero_function


class KinematicsBlock(object):
    dim_state = 13
    dim_output = 34
    dim_input = 7

    def __init__(self, gravity, winds, density, speed_of_sound, viscosity, omega_p=0., a=0., f=0.,):
        self.omega_p, self.a, self.f= omega_p, a, f
        if a<=0:
            self.gc2gd = get_flat_gc2gd()
            self.gd2gc = get_flat_gd2gc(a, f)
        else:
            self.gc2gd = get_nonflat_gc2gd(a, f)
            self.gd2gc = get_nonflat_gd2gc(a, f)

        self.gravity, self.winds, self.density, self.speed_of_sound, self.viscosity = gravity, winds, density, speed_of_sound, viscosity
    
    def prepare_to_integrate(self, *args, **kwargs):
        return

    def state_equation_function(self, t, x, u):
        return self.kinematics_state_function(t, *x, *u)
    
    def output_equation_function(self, t, x):
        return self.kinematics_output_function(t, *x)

    # TODO: need to be able to construct inertia to body angular velocities given local to body

    # below is from direct_eom output
    def kinematics_state_function(self, t, p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z, A_X, A_Y, A_Z, alpha_X, alpha_Y, alpha_Z, c_q):
        [g_x, g_y, g_z] = self.gravity(p_x,p_y,p_z).ravel()
        x0 = q_0**2
        x1 = q_1**2
        x2 = q_2**2
        x3 = q_3**2
        x4 = (x0 + x1 + x2 + x3)**(-1.0)
        x5 = 2*A_Y
        x6 = q_0*q_3
        x7 = q_1*q_2
        x8 = 2*A_Z
        x9 = q_0*x8
        x10 = q_3*x8
        x11 = 2*A_X
        x12 = (1/2)*omega_X
        x13 = (1/2)*q_2
        x14 = (1/2)*q_3
        x15 = (1/2)*omega_Y
        x16 = (1/2)*omega_Z
        return (numpy.array([v_x, v_y, v_z, x4*(A_X*x0 + A_X*x1 - A_X*x2 - A_X*x3 + g_x*x0 + g_x*x1 + g_x*x2 + g_x*x3 + q_1*x10 + q_2*x9 - x5*x6 + x5*x7), x4*(A_Y*x0 - A_Y*x1 + A_Y*x2 - A_Y*x3 + g_y*x0 + g_y*x1 + g_y*x2 + g_y*x3 - q_1*x9 + q_2*x10 + x11*x6 + x11*x7), x4*(A_Z*x0 - A_Z*x1 - A_Z*x2 + A_Z*x3 + g_z*x0 + g_z*x1 + g_z*x2 + g_z*x3 + q_0*q_1*x5 - q_0*q_2*x11 + q_1*q_3*x11 + q_2*q_3*x5), c_q*q_0 - omega_Y*x13 - omega_Z*x14 - q_1*x12, c_q*q_1 - omega_Y*x14 + omega_Z*x13 + q_0*x12, c_q*q_2 + q_0*x15 - q_1*x16 + q_3*x12, c_q*q_3 + q_0*x16 + q_1*x15 - q_2*x12, alpha_X, alpha_Y, alpha_Z]))

    def ic_from_geodetic(self, t, lamda_E, phi_E, h, V_N, V_E, V_D, psi, theta, phi, omega_X, omega_Y, omega_Z):
        [p_x, p_y, p_z] = self.gd2gc(lamda_E,phi_E,h).ravel()
        x0 = lamda_E + self.omega_p*t
        x1 = numpy.sin(x0)
        x2 = numpy.cos(x0)
        x3 = numpy.cos(phi_E)
        x4 = V_D*x3
        x5 = numpy.sin(phi_E)
        x6 = V_N*x5
        x7 = numpy.cos(phi)
        x8 = numpy.cos(lamda_E)
        x9 = numpy.cos(psi)
        x10 = x8*x9
        x11 = x10*x7
        x12 = numpy.sin(theta)
        x13 = x12*x3
        x14 = numpy.sin(phi)
        x15 = numpy.sin(psi)
        x16 = x15*x3
        x17 = numpy.cos(theta)
        x18 = x17*x5
        x19 = numpy.sin(lamda_E)
        x20 = x15*x17
        x21 = x14*x15
        x22 = x12*x8
        x23 = x13*x9
        x24 = x19*x7
        x25 = x15*x24
        x26 = x17*x3
        x27 = x14*x26
        x28 = x12*x5
        x29 = x19*x9
        x30 = x14*x29
        x31 = numpy.sqrt(-x10*x18 + x11 + x13*x8 + x14*x16 - x18*x7 - x19*x20 - x19*x27 + x21*x22 + x23*x7 + x25*x5 - x28*x30 + 1)
        x32 = x10*x14
        x33 = x21*x5
        x34 = x15*x7
        x35 = x29*x7
        x36 = (1/2)/x31
        return (numpy.array([p_x, p_y, p_z, -V_E*x1 - self.omega_p*p_y - x2*x4 - x2*x6, V_E*x2 + self.omega_p*p_x - x1*x4 - x1*x6, -V_D*x5 + V_N*x3, (1/2)*x31, x36*(-x14*x18 + x14*x23 - x16*x7 + x19*x33 - x22*x34 + x24*x26 + x28*x35 + x32), -x36*(x11*x28 + x12*x25 + x26*x7*x8 + x26*x9 + x28 - x30 + x33*x8), x36*(x12*x19*x21 + x13*x19 - x18*x29 + x20*x8 + x27*x8 + x28*x32 - x34*x5*x8 + x35), omega_X, omega_Y, omega_Z]))

    def kinematics_output_function(self, t, p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z):
        [lamda_I, phi_E, h] = self.gc2gd(p_x,p_y,p_z).ravel()
        lamda_E = lamda_I - self.omega_p*t
        ap_x = numpy.select([numpy.less_equal(self.a, 0.0)], [p_x], default=phi_E)
        ap_y = numpy.select([numpy.less_equal(self.a, 0.0)], [p_y], default=lamda_E)
        ap_z = h
        [g_x, g_y, g_z] = self.gravity(p_x,p_y,p_z).ravel()
        [W_N, W_E, W_D] = self.winds(t,ap_x,ap_y,ap_z).ravel()
        rho = self.density(t,ap_x,ap_y,ap_z)
        c_s = self.speed_of_sound(t,ap_x,ap_y,ap_z)
        mu = self.viscosity(t,ap_x,ap_y,ap_z)
        x0 = numpy.cos(lamda_E)
        x1 = self.omega_p*t
        x2 = numpy.cos(x1)
        x3 = x0*x2
        x4 = numpy.sin(lamda_E)
        x5 = numpy.sin(x1)
        x6 = x4*x5
        x7 = q_0*q_3
        x8 = q_1**2
        x9 = q_2**2
        x10 = q_0**2
        x11 = q_3**2
        x12 = x10 + x11
        x13 = (x12 + x8 + x9)**(-1.0)
        x14 = 2*x13
        x15 = x14*x7
        x16 = q_1*q_2
        x17 = x14*x16
        x18 = x15 + x17
        x19 = x0*x5
        x20 = x2*x4
        x21 = x13*x8
        x22 = x13*x9
        x23 = -x22
        x24 = x10*x13
        x25 = x11*x13
        x26 = x24 - x25
        x27 = x21 + x23 + x26
        x28 = numpy.sin(phi_E)
        x29 = q_0*q_2
        x30 = x14*x29
        x31 = q_1*q_3
        x32 = x14*x31
        x33 = -x30 + x32
        x34 = numpy.cos(phi_E)
        x35 = -x3*x34 + x34*x6
        x36 = -x19*x34 - x20*x34
        x37 = -x21
        x38 = q_0*q_1
        x39 = x14*x38
        x40 = q_2*q_3
        x41 = x14*x40
        x42 = self.omega_p**2
        x43 = lamda_E + x1
        x44 = numpy.cos(x43)
        x45 = W_E*x44
        x46 = 2*v_y
        x47 = W_N*x34
        x48 = 2*v_z
        x49 = self.omega_p*p_x
        x50 = W_D*x28
        x51 = numpy.sin(x43)
        x52 = W_E*x51
        x53 = 2*v_x
        x54 = self.omega_p*p_y
        x55 = W_D*x34
        x56 = x44*x55
        x57 = x51*x55
        x58 = 2*x49
        x59 = 2*x54
        x60 = W_N*x28
        x61 = x44*x60
        x62 = x51*x60
        x63 = W_D**2 + W_E**2 + W_N**2 + p_x**2*x42 + p_y**2*x42 + v_x**2 + v_y**2 + v_z**2 - x45*x46 + x45*x58 - x46*x49 + x46*x57 + x46*x62 - x47*x48 + x48*x50 + x52*x53 + x52*x59 + x53*x54 + x53*x56 + x53*x61 + x56*x59 - x57*x58 - x58*x62 + x59*x61
        x64 = numpy.sqrt(x63*((x63>0.0)+0.5*(x63==0.0)))
        x65 = v_z - x47 + x50
        x66 = -x8
        x67 = -x9
        x68 = v_x + x52 + x54 + x56 + x61
        x69 = 2*x68
        x70 = v_y - x45 - x49 + x57 + x62
        x71 = 2*x70
        x72 = x10 - x11
        x73 = 2*x65
        x74 = x13*(x69*(x16 - x7) + x70*(x66 + x72 + x9) + x73*(x38 + x40))/x64
        x75 = 2*self.omega_p
        x76 = v_x*x44
        x77 = v_y*x51
        x78 = x49*x51
        x79 = x44*x54
        return (numpy.array([p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z, lamda_E, phi_E, h, numpy.arctan2(x18*(x3 - x6) + x27*(-x19 - x20), x18*(-x19*x28 - x20*x28) + x27*(-x28*x3 + x28*x6) + x33*x34), -numpy.arcsin(x18*x36 + x27*x35 - x28*x33), numpy.arctan2(-x28*(x39 + x41) + x35*(-x15 + x17) + x36*(x22 + x26 + x37), -x28*(x23 + x24 + x25 + x37) + x35*(x30 + x32) + x36*(-x39 + x41)), rho, c_s, mu, x64, numpy.arctan2(x13*(x65*(x12 + x66 + x67) + x69*(x29 + x31) + x71*(-x38 + x40)), x13*(x68*(x67 + x72 + x8) + x71*(x16 + x7) + x73*(-x29 + x31))), numpy.arcsin(numpy.select([numpy.greater(x74, 1.0),numpy.less(x74, -1.0),numpy.greater(x64, 0.0)], [1.0,-1.0,x74], default=0.0)), x13*(omega_X*x10 + omega_X*x11 + omega_X*x8 + omega_X*x9 + x29*x75 - x31*x75), x13*(omega_Y*x10 + omega_Y*x11 + omega_Y*x8 + omega_Y*x9 - x38*x75 - x40*x75), x13*(omega_Z*x10 + omega_Z*x11 + omega_Z*x8 + omega_Z*x9 - self.omega_p*x10 - self.omega_p*x11 + self.omega_p*x8 + self.omega_p*x9), v_z*x34 - x28*x76 - x28*x77 + x28*x78 - x28*x79, -v_x*x51 + v_y*x44 - x44*x49 - x51*x54, -v_z*x28 - x34*x76 - x34*x77 + x34*x78 - x34*x79, W_N, W_E, W_D]))



class DynamicsBlock(object):
    dim_state = 0
    dim_output = 6
    dim_input = 49
    initial_condition = np.empty(0)

    def __init__(self, base_aero_coeffs, m, I_xx, I_yy, I_zz, I_xy, I_yz, I_xz, x_com, y_com, z_com, x_mrc, y_mrc, z_mrc, S_A, a_l, b_l, c_l, d_l):
        # TODO: make mrc, S_A, ref lengths, etc a property of aero model, mass stuff property of mass model?
        self.base_aero_coeffs = base_aero_coeffs
        self.m, self.I_xx, self.I_yy, self.I_zz, self.I_xy, self.I_yz, self.I_xz, self.x_com, self.y_com, self.z_com, self.x_mrc, self.y_mrc, self.z_mrc, self.S_A, self.a_l, self.b_l, self.c_l, self.d_l = m, I_xx, I_yy, I_zz, I_xy, I_yz, I_xz, x_com, y_com, z_com, x_mrc, y_mrc, z_mrc, S_A, a_l, b_l, c_l, d_l

    def prepare_to_integrate(self, *args, **kwargs):
        return

    def output_equation_function(self, t, x):
        return self.dynamics_output_function(t, *x)
    
    # below is from direct_eom output
    def tot_aero_forces_moments(self, qbar, V_T, alpha, beta, p_B, q_B, r_B, CD_b, CS_b, CL_b, CLcal_b, CMcal_b, CNcal_b, Cp_b, Cq_b, Cr_b, CD_e, CS_e, CL_e, CLcal_e, CMcal_e, CNcal_e, Cp_e, Cq_e, Cr_e):
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

    def dynamics_output_function(self, t, p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z, lamda_E, phi_E, h, psi, theta, phi, rho, c_s, mu, V_T, alpha, beta, p_B, q_B, r_B, V_N, V_E, V_D, W_N, W_E, W_D, Phi_x, Phi_y, Phi_z, tau_x, tau_y, tau_z, CD_e, CS_e, CL_e, CLcal_e, CMcal_e, CNcal_e, Cp_e, Cq_e, Cr_e):
        qbar = (1/2)*V_T**2*rho
        Ma = V_T/c_s
        Re = V_T*self.d_l*rho/mu
        [CD_b, CS_b, CL_b, CLcal_b, CMcal_b, CNcal_b, Cp_b, Cq_b, Cr_b] = self.base_aero_coeffs(alpha,beta,Ma,Re).ravel()
        [F_ax, F_ay, F_az, Lcal, Mcal, Ncal] = self.tot_aero_forces_moments(qbar,V_T,alpha,beta,p_B,q_B,r_B,CD_b,CS_b,CL_b,CLcal_b,CMcal_b,CNcal_b,Cp_b,Cq_b,Cr_b,CD_e,CS_e,CL_e,CLcal_e,CMcal_e,CNcal_e,Cp_e,Cq_e,Cr_e).ravel()
        [F_x, F_y, F_z, M_x, M_y, M_z] = numpy.array([F_ax + Phi_x, F_ay + Phi_y, F_az + Phi_z, Lcal + tau_x, Mcal + tau_y, Ncal + tau_z]).ravel()
        x0 = self.m**(-1.0)
        x1 = self.I_yz**2
        x2 = self.I_xz**2
        x3 = self.I_xy**2
        x4 = self.I_yy*self.I_zz
        x5 = self.I_xy*self.I_yz
        x6 = (self.I_xx*x1 - self.I_xx*x4 + 2*self.I_xz*x5 + self.I_yy*x2 + self.I_zz*x3)**(-1.0)
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

