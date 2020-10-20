from .version import __version__
import numpy
import numpy as np
import erfa
import fluids.atmosphere
from derivation import kinematics, dynamics

earth_spherical_gravity_constant = 14076443110000000/(3.28084**3)
earth_f = 1/298.257223563
earth_J2_constant = 0.00108262982
earth_rotation_rate_deg = 0.004178073
earth_rotation_rate = earth_rotation_rate_deg*np.pi/180
earth_equitorial_radius = 6378137.0

# kinematic functions
def get_nonflat_gc2gd(a, f):
    def gc2gde(px, py, pz):
        return np.array(erfa.gc2gde(a, f, np.array([px, py, pz])))
    return gc2gde

def get_flat_gc2gd():
    def gc2gde(px, py, pz):
        return np.array([px, py, pz])
    return gc2gde

def get_nonflat_gd2gc(a, f):
    def gd2gce(longitude, latitude, altitude):
        return np.array(erfa.gd2gce(a, f, longitude, latitude, altitude))
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
    """
    The KinematicsBlock is the state dynamics block that provides the kinematic equations of motion for integration and an
    output equation for commonly used variables for flight vehicles according to the planet model. The KinematicsBlock 
    planet model is parameterized based on the following components:

    gravity model: translational acceleration due to gravity as a function of planet-fixed position in rectangular coordinates

    atmospheric models: density, speed of sound, and viscocity outputs of the atmosphere model as a function of time (i.e.,
    for stochasticity) and position in planet-fixed frame

    planetodetic model parameters: semi-major axis length, flattening ratio
    TODO: replace this with to_detic and to_centric transformation functions

    planet rotation rate: used to adjust the local wind and to calculate a simple rotation transformation between planet-fixed
    and planet-centered inertial
    TODO: add generic rotation matrix call (i.e., to allow models including precession and nutation)

    TODO: document EOM's?


    The state components are:
        [0:3] p_x, p_y, p_z
        translational position (of vehicle center of mass) relative to inertial origin expressed in inertial coordinate system

        [3:6] v_x, v_y, v_z
        translational velocity (of vehicle center of mass) in the inertial coordinate system expressed in inertial coordinates

        [6:10] q_0, q_1, q_2, q_3
        quaternion components representing the rotation from the inertial to the body-fixed coordinate systems

        [10:13] omega_X, omega_Y, omega_Z
        angular velocity from inertial coordinate system to body-fixed coordinate system expressed in body-fixed coordinates

    The input components are:
        [0:3] A_X, A_Y, A_Z
        translational acceleration of vehicle center of mass due to non-gravitaitonal forces in the inertial coordinate system 
        expressed in body-fixed coordinates
              
        [3:6] alpha_X, alpha_Y, alpha_Z
        angular acceleration of vehicle coordinate system in the inertial coordinate system expressed in body-fixed coordinates

        [6] c_q
        scaling parameter for non-unit quaternion kinematics integration; this is a control for numerical behavior of the integration
        and can be ignored in most circumstances

    The output components are:
        [0:13] p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z
        The state components listed above

        [13:16] lamda_E, phi_E, h
        translational position of the vehicle center of mass in the planet-fixed frame expressed in
        planetodetic spherical position coordinates: lattitude, longitude, and altitude (TODO: verify order)

        [16:19] psi, theta, phi
        euler-angles relating the NED frame to the body-fixed frame: yaw, pitch, roll

        [19:22] rho, c_s, mu
        density, speed of sound, and viscocity outputs of the atmosphere model as a function of time and position in planet-fixed frame

        [22:25] V_T, alpha, beta
        true air speed, angle of attack, and angle of sideslip used for aerodynamics (includes wind)

        [25:28] p_B, q_B, r_B
        angular velocity components of vehicle-fixed coordinate system relative to NED, used for aerodynamic damping derivatives

        [28:31] V_N, V_E, V_D
        velocity of the vehicle center of mass in the planet-fixed frame expressed in NED frame

        [31:34] W_N, W_E, W_D
        Output of wind model as a function of time and position in planet-fixed frame

    """
    dim_state = 13
    dim_output = 34
    dim_input = 7

    def __init__(self, gravity, winds, density, speed_of_sound, viscosity, omega_p=0., a=0., f=0.,):
        self.omega_p, self.a, self.f = omega_p, a, f
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

    def kinematics_state_function(self, t, p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z, A_X, A_Y, A_Z, alpha_X, alpha_Y, alpha_Z, c_q=0.,):
        return kinematics.kinematics_state_function(self, t, p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z, A_X, A_Y, A_Z, alpha_X, alpha_Y, alpha_Z, c_q)

    def ic_from_geodetic(self, t, lamda_E, phi_E, h, V_N, V_E, V_D, psi, theta, phi, omega_X, omega_Y, omega_Z):
        return kinematics.ic_from_geodetic(self, t, lamda_E, phi_E, h, V_N, V_E, V_D, psi, theta, phi, omega_X, omega_Y, omega_Z)

    def kinematics_output_function(self, t, p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z):
        return kinematics.kinematics_output_function(self, t, p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z)


class DynamicsBlock(object):
    """
    The DynamicsBlock is a state-less dynamics block that provides a common calculation for the dynamics of a flight vehicle. The DynamicsBlock 
    vehicle model is parameterized based on the following components:

    Inertial properties: m, I_xx, I_yy, I_zz, I_xy, I_yz, I_xz, x_com, y_com, z_com
    Inertia is expressed in body-fixed coordinates about the center of mass
    center of mass (x_com, y_com, z_com) position is relative to an arbitrary reference such as the CAD origin
    TODO: verify sign of products of inertia

    Aerodynamics properties: base_aero_coeffs, S_A, a_l, b_l, c_l, d_l, x_mrc, y_mrc, z_mrc,
    base_aero_coeffs is assumed to be a function of angles of attack and sideslip, and Mach and Reynolds number
    Should return coefficients of:
        drag, sideforce, and lift (force coefficients expressed in wind coordinate system), 
        static roll, pitch, and yaw moments (moment coefficients expressed in body-fixed coordinate system), and
        dynamic or damping roll, pitch, and yaw moments (also moment coefficients expressed in body-fixed coordinate system)
    
    S_A is the reference area
    a_l, b_l, and c_l are reference lengths for roll, pitch, and yaw
    d_l is the reference length for Reynolds number calculation

    Moment coefficients are assumed to be about the moment reference center (x_com, y_com, z_com) the position of 
    which is defined relative to the same arbitrary reference as the center of mass. Damping coefficients are transformed
    with the static moment coefficients which does not usually hold; to model damping coefficients it is recommended
    to set the center of mass and moment reference center to the same location at the arbitrary reference (i.e., [0,0,0])



    The input components are:
        [0:13] p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z
        The vehicle state components

        [13:16] lamda_E, phi_E, h
        translational position of the vehicle center of mass in the planet-fixed frame expressed in
        planetodetic spherical position coordinates: latitude, longitude, and altitude (TODO: verify order)

        [16:19] psi, theta, phi
        euler-angles relating the NED frame to the body-fixed frame: yaw, pitch, roll

        [19:22] rho, c_s, mu
        density, speed of sound, and viscocity outputs of the atmosphere model as a function of time and position in planet-fixed frame

        [22:25] V_T, alpha, beta
        true air speed, angle of attack, and angle of sideslip used for aerodynamics (includes wind)

        [25:28] p_B, q_B, r_B
        angular velocity components of body relative to NED used for aerodynamic damping derivatives

        [28:31] V_N, V_E, V_D
        velocity of the vehicle center of mass in the planet-fixed frame expressed in NED frame

        [31:34] W_N, W_E, W_D
        winds acting on the vehicle expressed in NED frame

        [34:37] Phi_x, Phi_y, Phi_z
        Extra force acting on vehicle (center of mass) expressed in body-fixed coordinates

        [37:40] tau_x, tau_y, tau_z
        Extra moment acting on vehicle expressed in body-fixed coordinates

        [40:49] CD_e, CS_e, CL_e, CLcal_e, CMcal_e, CNcal_e, Cp_e, Cq_e, Cr_e
        Extra aerodynamic coefficient components (additive)

    The output components are:
        [0:3] A_X, A_Y, A_Z
        translational acceleration of vehicle center of mass due to non-gravitaitonal forces in the inertial coordinate system 
        expressed in body-fixed coordinates, computed assuming constant mass and total non-gravitational forces due to aerodynamics
        (base model and extra aerodynamic coefficients input components) and the extra force input components.
              
        [3:6] alpha_X, alpha_Y, alpha_Z
        angular acceleration of vehicle coordinate system in the inertial coordinate system expressed in body-fixed coordinates, 
        computed assuming constant inertia and total moments due to aerodynamics
        (base model and extra aerodynamic coefficients input components) and the extra moment input components.
    """
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
        return dynamics.tot_aero_forces_moments(self, qbar, V_T, alpha, beta, p_B, q_B, r_B, CD_b, CS_b, CL_b, CLcal_b, CMcal_b, CNcal_b, Cp_b, Cq_b, Cr_b, CD_e, CS_e, CL_e, CLcal_e, CMcal_e, CNcal_e, Cp_e, Cq_e, Cr_e)

    def dynamics_output_function(self, t, p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z, lamda_E, phi_E, h, psi, theta, phi, rho, c_s, mu, V_T, alpha, beta, p_B, q_B, r_B, V_N, V_E, V_D, W_N, W_E, W_D, Phi_x, Phi_y, Phi_z, tau_x, tau_y, tau_z, CD_e, CS_e, CL_e, CLcal_e, CMcal_e, CNcal_e, Cp_e, Cq_e, Cr_e):
        return dynamics.dynamics_output_function(self, t, p_x, p_y, p_z, v_x, v_y, v_z, q_0, q_1, q_2, q_3, omega_X, omega_Y, omega_Z, lamda_E, phi_E, h, psi, theta, phi, rho, c_s, mu, V_T, alpha, beta, p_B, q_B, r_B, V_N, V_E, V_D, W_N, W_E, W_D, Phi_x, Phi_y, Phi_z, tau_x, tau_y, tau_z, CD_e, CS_e, CL_e, CLcal_e, CMcal_e, CNcal_e, Cp_e, Cq_e, Cr_e)

