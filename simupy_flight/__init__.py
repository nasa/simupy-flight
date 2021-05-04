from .version import __version__
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
earth_mean_radius = 6371007.1809

# kinematic functions
def get_nonflat_pc2pd(a, f):
    def pcf2pd(px, py, pz):
        return np.array(erfa.gc2gde(a, f, np.array([px, py, pz])))
    return pcf2pd

def get_flat_pc2pd():
    def pcf2pd(px, py, pz):
        return np.array([0.0, np.pi/2, pz])
    return pcf2pd

def get_nonflat_pd2pc(a, f):
    def pd2pcf(longitude, latitude, altitude):
        return np.array(erfa.gd2gce(a, f, longitude, latitude, altitude))
    return pd2pcf

def get_flat_pd2pc():
    def pd2pcf(longitude, latitude, altitude):
        return np.array([longitude, latitude, altitude])
    return pd2pcf

class Planetodetic:
    def __init__(self, a=0., f=0., omega_p=0.):
        self.omega_p, self.a, self.f = omega_p, a, f
        if a<=0:
            self.pcf2pd = get_flat_pc2pd()
            self.pd2pcf = get_flat_pd2pc()
        else:
            self.pcf2pd = get_nonflat_pc2pd(a, f)
            self.pd2pcf = get_nonflat_pd2pc(a, f)


# constant atmosphere functions
def get_constant_atmosphere(density_val=0., speed_of_sound_val=1., viscocity_val=1.,):
    atmosphere_val = np.array([density_val, speed_of_sound_val, viscocity_val])
    def atmosphere_constant(t,ap_x,ap_y,ap_z):
        return atmosphere_val
    return atmosphere_constant

def atmosphere_1976(t,ap_x,ap_y,ap_z):
    atmo = fluids.atmosphere.ATMOSPHERE_1976(ap_z)
    return np.array([atmo.rho, atmo.v_sonic, atmo.mu])

# gravity functions
# these are "gravitation" functions according to WGS-84 notation
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

def get_constant_force_moments(FX=0., FY=0., FZ=0., MX=0., MY=0., MZ=0.,):
    force_moment_vals = np.array([FX, FY, FZ, MX, MY, MZ])
    def force_moment_function(*args):
        return force_moment_vals
    return force_moment_function

class Planet(object):
    """
    The Planet is the state dynamics block that provides the kinematic equations of motion for integration and an
    output equation for commonly used variables for flight vehicles according to the planet model. The Planet 
    planet model is parameterized based on the following components:

    ``gravity`` model: translational acceleration due to gravity as a function of planet-fixed position in rectangular coordinates.
    For 

    ``atmosphere`` models: density, speed of sound, and viscocity outputs of the atmosphere model as a function of time (i.e.,
    for stochasticity) and position in planet-fixed frame

    ``winds`` model: wind in local NED frame as a function of time and planetodetic position. Positive wind indicates wind in specified
    direction so wind "from west" is a positive W_E component.

    ``planetodetic`` model, must provide rotation rate in rad/s as ``omega_p`` and functions ``pd2pcf`` and ``pcf2pd`` to
    convert between planetodetic rectangular coordinates and planetocentric spherical coordinates. Future versions may support
    nutation and precession. The planetodetic model should assume pcf excludes sidereel rotation
    which is accounted for by the kinematics model.

    The state components are:
        [0:3] p_x, p_y, p_z
        translational position (of vehicle center of mass) relative to inertial origin expressed in inertial coordinate system

        [3:7] q_0, q_1, q_2, q_3
        quaternion components representing the rotation from the inertial to the body-fixed coordinate systems

        [7:10] v_x, v_y, v_z
        translational velocity (of vehicle center of mass) in the inertial coordinate system expressed in inertial coordinates

        [10:13] omega_X, omega_Y, omega_Z
        angular velocity from inertial coordinate system to body-fixed coordinate system expressed in body-fixed coordinates

    The input components are:
        [0:3] A_X, A_Y, A_Z
        translational acceleration of vehicle center of mass due to non-gravitaitonal forces in the inertial coordinate system 
        expressed in body-fixed Forward-Right-Down (FRD) coordinate system
              
        [3:6] alpha_X, alpha_Y, alpha_Z
        angular acceleration of vehicle coordinate system in the inertial coordinate system expressed in body-fixed FRD coordinates

        [6] c_q
        scaling parameter for non-unit quaternion kinematics integration; this is a control for numerical behavior of the integration
        and can be ignored in most circumstances

    The output components are:
        [0:13] p_x, p_y, p_z, q_0, q_1, q_2, q_3, v_x, v_y, v_z, omega_X, omega_Y, omega_Z
        The state components listed above

        [13:16] lamda_E, phi_E, h
        translational position of the vehicle center of mass in the planet-fixed frame expressed in
        planetodetic spherical position coordinates: longitude, latitude, and altitude

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

    p_x_idx = 0
    p_y_idx = 1
    p_z_idx = 2
    q_0_idx = 3
    q_1_idx = 4
    q_2_idx = 5
    q_3_idx = 6
    v_x_idx = 7
    v_y_idx = 8
    v_z_idx = 9
    omega_X_idx = 10
    omega_Y_idx = 11
    omega_Z_idx = 12
    lamda_D_idx = 13
    phi_D_idx = 14
    h_D_idx = 15
    psi_idx = 16
    theta_idx = 17
    phi_idx = 18
    rho_idx = 19
    c_s_idx = 20
    mu_idx = 21
    V_T_idx = 22
    alpha_idx = 23
    beta_idx = 24
    p_B_idx = 25
    q_B_idx = 26
    r_B_idx = 27
    V_N_idx = 28
    V_E_idx = 29
    V_D_idx = 30
    W_N_idx = 31
    W_E_idx = 32
    W_D_idx = 33

    output_column_names = ['p_x', 'p_y', 'p_z',
        'q_0', 'q_1', 'q_2', 'q_3', 
        'v_x', 'v_y', 'v_z', 
        'omega_X', 'omega_Y', 'omega_Z', 
        'lamda_D', 'phi_D', 'h_D',
        'psi', 'theta', 'phi', 
        'rho', 'c_s', 'mu',
        'V_T', 'alpha', 'beta',
        'p_B', 'q_B', 'r_B',
        'V_N', 'V_E', 'V_D',
        'W_N', 'W_E', 'W_D',]
    
    output_column_names_latex = ['$p_{x}$', '$p_{y}$', '$p_{z}$', 
        '$q_{0}$', '$q_{1}$', '$q_{2}$', '$q_{3}$',
        '$v_{x}$', '$v_{y}$', '$v_{z}$', 
        '$\\omega_{X}$', '$\\omega_{Y}$', '$\\omega_{Z}$', 
        '$\\lambda_{D}$', '$\\phi_{D}$', '$h_{D}$', 
        '$\\psi$', '$\\theta$', '$\\phi$', 
        '$\\rho$', '$c_{s}$', '$\\mu$', 
        '$V_{T}$', '$\\alpha$', '$\\beta$', 
        '$p_{B}$', '$q_{B}$', '$r_{B}$',
        '$V_{N}$', '$V_{E}$', '$V_{D}$',
        '$W_{N}$', '$W_{E}$', '$W_{D}$',]

    def __init__(self, gravity, winds, atmosphere, planetodetics):
        (self.gravity, self.winds, self.atmosphere, self.planetodetics) =\
            gravity, winds, atmosphere, planetodetics
    
    def prepare_to_integrate(self, *args, **kwargs):
        return

    def state_equation_function(self, t, x, u):
        return kinematics.kinematics_state_function(self, t, *x, *u)
    
    def output_equation_function(self, t, x):
        return kinematics.kinematics_output_function(self, t, *x)

    def ic_from_planetodetic(self, lamda_E=0., phi_E=0., h=0., V_N=0., V_E=0., V_D=0., psi=0., theta=0., phi=0., p_B=0., q_B=0., r_B=0.):
        return kinematics.ic_from_planetodetic(self, lamda_E, phi_E, h, V_N, V_E, V_D, psi, theta, phi, p_B, q_B, r_B)
    
    def local_translational_trim_residual(self, p_x, p_y, p_z, q_0, q_1, q_2, q_3, v_x, v_y, v_z, A_X, A_Y, A_Z):
        return kinematics.local_translational_trim_residual(p_x, p_y, p_z, q_0, q_1, q_2, q_3, v_x, v_y, v_z, A_X, A_Y, A_Z)


class Vehicle(object):
    """
    The Vehicle is a state-less dynamics block that provides a common calculation for the dynamics of a flight vehicle. The Vehicle 
    vehicle model is parameterized based on the following components:

    Inertial properties: m, I_xx, I_yy, I_zz, I_xy, I_yz, I_xz, x_com, y_com, z_com
    Inertia is expressed in body-fixed coordinates about the center of mass
    center of mass (x_com, y_com, z_com) position is relative to an arbitrary reference such as the CAD origin
    TODO: verify sign of products of inertia

    Aerodynamics properties: base_aero_coeffs, S_A, a_l, b_l, c_l, d_l, x_mrc, y_mrc, z_mrc,
    base_aero_coeffs is assumed to be a function of angles of attack and sideslip, and Mach and Reynolds number
    Should return coefficients of:
        drag, sideforce, and lift (force coefficients expressed in wind coordinate system), 
        static roll, pitch, and yaw moments (moment coefficients expressed in body-fixed Forward-Right-Down coordinate system), and
        dynamic or damping roll, pitch, and yaw moments (also moment coefficients expressed in body-fixed coordinate system)
    
    S_A is the reference area
    a_l, b_l, and c_l are reference lengths for roll, pitch, and yaw
    d_l is the reference length for Reynolds number calculation

    Moment coefficients are assumed to be about the moment reference center (x_com, y_com, z_com) the position of 
    which is defined relative to the same arbitrary reference as the center of mass. Damping coefficients are transformed
    with the static moment coefficients which does not usually hold; to model damping coefficients it is recommended
    to set the center of mass and moment reference center to the same location at the arbitrary reference (i.e., [0,0,0])

    Processing of ``input_aero_coeffs`` and ``input_force_moment`` parameterized models follow the same logic:
        if ``None``: assume that a child class will over-write
        if callable: use directly (so it should have the correct signature)
        else: try to build a function that returns constant value(s)

    Attributes ``input_aero_coeffs_idx`` and ``input_force_moment_idx`` are used to route extra (control) inputs to the input aero and foce/moment
    functions respectively. Use ``None`` to route all inputs, use a list of integers to route particular inputs including an empty list
    for no routing.


    The input components are:
        [0:13] p_x, p_y, p_z, q_0, q_1, q_2, q_3, v_x, v_y, v_z, omega_X, omega_Y, omega_Z
        The vehicle state components

        [13:16] lamda_E, phi_E, h
        translational position of the vehicle center of mass in the planet-fixed frame expressed in
        planetodetic spherical position coordinates: longitude, latitude, and altitude

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

    The output components are:
        [0:3] A_X, A_Y, A_Z
        translational acceleration of vehicle center of mass due to non-gravitaitonal forces in the inertial coordinate system 
        expressed in body-fixed Forward-Right-Down (FRD) coordinates, computed assuming constant mass and total non-gravitational forces due to aerodynamics
        (base model and extra aerodynamic coefficients input components) and the extra force input components.
              
        [3:6] alpha_X, alpha_Y, alpha_Z
        angular acceleration of vehicle coordinate system in the inertial coordinate system expressed in body-fixed FRD coordinates, 
        about the center of mass computed assuming constant inertia and total moments due to aerodynamics
        (base model and extra aerodynamic coefficients input components) and the extra moment input components.

    """
    dim_state = 0
    dim_output = 6
    dim_input = 34

    initial_condition = np.empty(0)

    output_column_names = ['A_X', 'A_Y', 'A_Z', 
        'alpha_X', 'alpha_Y', 'alpha_Z']

    output_column_names_latex = ['$A_{X}$', '$A_{Y}$', '$A_{Z}$', 
        '$\\alpha_{X}$', '$\\alpha_{Y}$', '$\\alpha_{Z}$']

    t_arg_idx = 0
    p_x_arg_idx = 1
    p_y_arg_idx = 2
    p_z_arg_idx = 3
    q_0_arg_idx = 4
    q_1_arg_idx = 5
    q_2_arg_idx = 6
    q_3_arg_idx = 7
    v_x_arg_idx = 8
    v_y_arg_idx = 9
    v_z_arg_idx = 10
    omega_X_arg_idx = 11
    omega_Y_arg_idx = 12
    omega_Z_arg_idx = 13
    lamda_D_arg_idx = 14
    phi_D_arg_idx = 15
    h_D_arg_idx = 16
    psi_arg_idx = 17
    theta_arg_idx = 18
    phi_arg_idx = 19
    rho_arg_idx = 20
    c_s_arg_idx = 21
    mu_arg_idx = 22
    V_T_arg_idx = 23
    alpha_arg_idx = 24
    beta_arg_idx = 25
    p_B_arg_idx = 26
    q_B_arg_idx = 27
    r_B_arg_idx = 28
    V_N_arg_idx = 29
    V_E_arg_idx = 30
    V_D_arg_idx = 31
    W_N_arg_idx = 32
    W_E_arg_idx = 33
    W_D_arg_idx = 34
    qbar_arg_idx = 35
    Ma_arg_idx = 36
    Re_arg_idx = 37
    starargs_idx = 38

    def __init__(self,  
            m, I_xx, I_yy, I_zz, I_xy, I_yz, I_xz, x_com, y_com, z_com, 
            base_aero_coeffs, x_mrc, y_mrc, z_mrc, S_A, a_l, b_l, c_l, d_l, 
            input_aero_coeffs=0.0, input_force_moment=0.0,
            input_aero_coeffs_idx=None, input_force_moment_idx=None,
            dim_additional_input=0,
            ):
        
        self.dim_input = Vehicle.dim_input + dim_additional_input

        if callable(base_aero_coeffs):
            pass
        else:
            base_aero_coeff_vals = np.array(base_aero_coeffs, dtype=np.float_).reshape(-1)
            if base_aero_coeff_vals.size == 1:
                base_aero_coeff_vals = np.tile(base_aero_coeff_vals, 9)
            base_aero_coeffs = get_constant_aero(*base_aero_coeff_vals)
        self.base_aero_coeffs = base_aero_coeffs

        if callable(input_aero_coeffs):
            pass
        else:
            input_aero_coeff_vals = np.array(input_aero_coeffs, dtype=np.float_).reshape(-1)
            if input_aero_coeff_vals.size == 1:
                input_aero_coeff_vals = np.tile(input_aero_coeff_vals, 9)
            input_aero_coeffs = get_constant_aero(*input_aero_coeff_vals)
        self.input_aero_coeffs = input_aero_coeffs

        if callable(input_force_moment):
            pass
        else:
            input_force_moment_vals = np.array(input_force_moment, dtype=np.float_).reshape(-1)
            if input_force_moment_vals.size == 1:
                input_force_moment_vals = np.tile(input_force_moment_vals, 6)
            input_force_moment = get_constant_force_moments(*input_force_moment_vals)
        self.input_force_moment = input_force_moment
        
        if input_aero_coeffs_idx is None:
            input_aero_coeffs_idx = np.arange(self.dim_input - Vehicle.dim_input)
        # TODO: defensive checking for routing indices
        self.input_aero_coeffs_idx = input_aero_coeffs_idx

        if input_force_moment_idx is None:
            input_force_moment_idx = np.arange(self.dim_input - Vehicle.dim_input)
        self.input_force_moment_idx = input_force_moment_idx

        self.m, self.I_xx, self.I_yy, self.I_zz, self.I_xy, self.I_yz, self.I_xz, self.x_com, self.y_com, self.z_com, self.x_mrc, self.y_mrc, self.z_mrc, self.S_A, self.a_l, self.b_l, self.c_l, self.d_l = m, I_xx, I_yy, I_zz, I_xy, I_yz, I_xz, x_com, y_com, z_com, x_mrc, y_mrc, z_mrc, S_A, a_l, b_l, c_l, d_l

    def prepare_to_integrate(self, *args, **kwargs):
        return
    
    def _input_aero_coeffs(self, alpha, beta, Ma, Re, *args):
        filtered_args = np.array(args)[self.input_aero_coeffs_idx]
        return self.input_aero_coeffs(alpha, beta, Ma, Re, *filtered_args)
    
    def _input_force_moment(self,t,p_x,p_y,p_z,q_0,q_1,q_2,q_3,v_x,v_y,v_z,omega_X,omega_Y,omega_Z,lamda_D,phi_D,h_D,psi,theta,phi,rho,c_s,mu,V_T,alpha,beta,p_B,q_B,r_B,V_N,V_E,V_D,W_N,W_E,W_D,qbar,Ma,Re,*args):
        filtered_args = np.array(args)[self.input_force_moment_idx]
        return self.input_force_moment(t,p_x,p_y,p_z,v_x,v_y,v_z,q_0,q_1,q_2,q_3,omega_X,omega_Y,omega_Z,lamda_D,phi_D,h_D,psi,theta,phi,rho,c_s,mu,V_T,alpha,beta,p_B,q_B,r_B,V_N,V_E,V_D,W_N,W_E,W_D,qbar,Ma,Re,*filtered_args)

    def output_equation_function(self, t, u):
        # TODO: test that an inherited class that overwrites dim_input still has access to Vehicle.dim_input for this to work.
        uu = u[..., :Vehicle.dim_input]
        u_extra = u[..., Vehicle.dim_input:]
        return dynamics.dynamics_output_function(self, t, *uu, *u_extra)

    def tot_aero_forces_moments(self, qbar, Ma, Re, V_T, alpha, beta, p_B, q_B, r_B, *args):
        return dynamics.tot_aero_forces_moments(self, qbar, Ma, Re, V_T, alpha, beta, p_B, q_B, r_B, *args)
