from .version import __version__
import numpy as np
import erfa
import fluids.atmosphere
from . import kinematics, dynamics

earth_spherical_gravity_constant = 14076443110000000 / (3.28084**3)
earth_f = 1 / 298.257223563
earth_J2_constant = 0.00108262982
earth_rotation_rate_deg = 0.004178073
earth_rotation_rate = earth_rotation_rate_deg * np.pi / 180
earth_equitorial_radius = 6378137.0
earth_mean_radius = 6371007.1809

# kinematic functions
def get_nonflat_pc2pd(a, f):
    """
    Generate a planetocentric-to-planetodetic coordinate transformation function for an
    ellipsoidal planetary model, parameterized by equatorial radius `a` and flattening
    `f`.

    The transformation function takes in the planetocentric rectangular coordinates
    `px`, `py`, `pz` and returns the planetodetic `longitude`, `latitude`,
    and `altitude`. The units of `a`, `px`, `py`, `pz`, and `altitude` are all the same
    and assumed to be meters to be consistent with other kinematics.

    Uses the ERFA implementation of Fukushima's method.

    Parameters
    ----------
    a : float
        equatorial radius
    f : float
        flattening parameter

    Returns
    -------
    pcf2pd : callable
        The planetocentric-to-planetodetic coordinate transformation function.
    """

    def pcf2pd(px, py, pz):
        return np.array(erfa.gc2gde(a, f, np.array([px, py, pz])))

    return pcf2pd


def get_flat_pc2pd():
    """
    Generate a planetocentric-to-planetodetic coordinate transformation function for a
    flat planetary model.

    The transformation function takes in the planetocentric rectangular coordinates
    `px`, `py`, `pz` and returns a constant 0 for `longitude`, -pi/2 for `latitude`, and
    uses `pz` for `altitude`.

    Returns
    -------
    pcf2pd : callable
        The planetocentric-to-planetodetic coordinate transformation function.
    """

    def pcf2pd(px, py, pz):
        return np.array([0.0, -np.pi / 2, pz])

    return pcf2pd


def get_nonflat_pd2pc(a, f):
    """
    Generate a planetodetic-to-planetocentric coordinate transformation function for an
    ellipsoidal planetary model, parameterized by equatorial radius `a` and flattening
    `f`.

    The transformation function takes in the position defined by the planetodetic
    `longitude`, `latitude`, and `altitude` and returns the planetocentric rectangular
    coordinates `px`, `py`, `pz`.  The units of `a`, `px`, `py`, `pz`, and `altitude`
    are all the same and assumed to be meters to be consistent with other kinematics.

    Uses the ERFA implementation of Fukushima's method.

    Parameters
    ----------
    a : float
        equatorial radius
    f : float
        flattening parameter

    Returns
    -------
    pd2pcf : callable
        The planetodetic-to-planetocentric coordinate transformation function.
    """

    def pd2pcf(longitude, latitude, altitude):
        return np.array(erfa.gd2gce(a, f, longitude, latitude, altitude))

    return pd2pcf


def get_flat_pd2pc():
    """
    Generate a planetodetic-to-planetocentric coordinate transformation function for a
    flat planetary model.

    The transformation function takes in the planetocentric rectangular coordinates
    `px`, `py`, `pz` and returns a constant 0 longitude, -90 degree latitude, and uses
    `pz` for altitude.

    Returns
    -------
    pd2pcf : callable
        The planetocentric to planetodetic coordinate transformation function.
    """

    def pd2pcf(longitude, latitude, altitude):
        return np.array([longitude, latitude, altitude])

    return pd2pcf


class Planetodetic:
    """
    Namespace class used to model the planetodetics used by the Planet class.

    Parameters
    ----------
    a : float
        equatorial radius
    f : float
        flattening parameter
    omega_p : float
        Angular rate of planetary rotation


    .. note::
        Use `a<=0` and `omega_p=0` for a flat planet model; `f` will be ignored
    """

    def __init__(self, a=0.0, f=0.0, omega_p=0.0):
        self.omega_p, self.a, self.f = omega_p, a, f
        if a <= 0:
            self.pcf2pd = get_flat_pc2pd()
            self.pd2pcf = get_flat_pd2pc()
        else:
            self.pcf2pd = get_nonflat_pc2pd(a, f)
            self.pd2pcf = get_nonflat_pd2pc(a, f)


# constant atmosphere functions
def get_constant_atmosphere(
    density_val=0.0,
    speed_of_sound_val=1.0,
    viscocity_val=1.0,
):
    """
    Generate an atmosphere callback that returns a constant and specified density, speed
    of sound, and viscosity regardless of the input time or planetodetic position
    """
    atmosphere_val = np.array([density_val, speed_of_sound_val, viscocity_val])

    def atmosphere_constant(t, longitude, latitude, altitude):
        return atmosphere_val

    return atmosphere_constant


def atmosphere_1976(t, longitude, latitude, altitude):
    """Atmospheric model callback for the US Standard 1976 Atmosphere."""
    atmo = fluids.atmosphere.ATMOSPHERE_1976(altitude)
    return np.array([atmo.rho, atmo.v_sonic, atmo.mu])


# gravity functions
# these are "gravitation" functions according to WGS-84 notation
def get_spherical_gravity(gravitational_constant):
    """
    Generate a spherical gravitational model with the specified gravitational constant

    units?
    """

    def gravity_function(px, py, pz):
        pos_vec = np.array([px, py, pz])
        return -gravitational_constant * pos_vec / np.linalg.norm(pos_vec) ** 3

    return gravity_function


def earth_J2_gravity(px, py, pz):
    """
    Gravitational model callback for the J2 gravity model
    """
    pos_vec = np.array([px, py, pz])
    r = np.linalg.norm(pos_vec)
    return (
        -earth_spherical_gravity_constant
        * (pos_vec / r**3)
        * (
            1
            - 3
            * earth_J2_constant
            * earth_equitorial_radius**2
            * (5 * pz**2 - r**2)
            / (2 * r**4)
        )
    )


# wind functions
def get_constant_winds(wx=0.0, wy=0.0, wz=0.0):
    """
    Generate a wind model callback that returns a constant, specified wind vector
    regardless of the input time or geodetic position
    """
    wind_val = np.array([wx, wy, wz])

    def winds_function(t, ap_x, ap_y, ap_z):
        return wind_val

    return winds_function


# aero functions
def get_constant_aero(
    CD_b=0.0,
    CS_b=0.0,
    CL_b=0.0,
    CLcal_b=0.0,
    CMcal_b=0.0,
    CNcal_b=0.0,
    Cp_b=0.0,
    Cq_b=0.0,
    Cr_b=0.0,
):
    """
    Generate an aerodynamics model callback that returns a constant and specified set of
    force, moment, and damping coefficients regardless of flight condition.
    """
    aero_vals = np.array(
        [CD_b, CS_b, CL_b, CLcal_b, CMcal_b, CNcal_b, Cp_b, Cq_b, Cr_b]
    )

    def aero_function(alpha, beta, Ma, Re):
        return aero_vals

    return aero_function


def get_constant_force_moments(
    FX=0.0,
    FY=0.0,
    FZ=0.0,
    MX=0.0,
    MY=0.0,
    MZ=0.0,
):
    """
    Generate a force and moment callback that returns a constant and specified set of
    forces and moments regardless of flight condition.
    """
    force_moment_vals = np.array([FX, FY, FZ, MX, MY, MZ])

    def force_moment_function(*args):
        return force_moment_vals

    return force_moment_function


class Planet(object):
    """
    The Planet is the state dynamics block that provides the kinematic equations of
    motion for integration and an output equation for commonly used variables for flight
    vehicles according to the planet model.

    The planet model is parameterized based on the following components:

    `gravity` model:
        translational acceleration due to gravity as a function of planet-fixed position
        in rectangular coordinates. See :func:`earth_J2_gravity` for example.
    `winds` model:
        wind in local NED frame as a function of time and planetodetic position.
        Positive wind indicates wind in specified direction so wind "from west" is a
        positive W_E component. See :func:`get_constant_winds` for generating a simple
        wind model.
    `atmosphere` model:
        density, speed of sound, and viscosity outputs of the atmosphere model as a
        function of time (i.e. for stochasticity) and position in planet-fixed frame
        expressed in geodetic coordinates (longitude, latitude, altitude). See
        :func:`get_constant_atmosphere` for generating a simple atmosphere model and
        :func:`atmosphere_1976` for the US Standard 1976 Atmosphere.
    `planetodetic` model:
        must provide rotation rate in rad/s as ``omega_p`` and functions ``pd2pcf`` and
        ``pcf2pd`` to convert between planetodetic rectangular coordinates and
        planetocentric spherical coordinates. Future versions may support nutation and
        precession. The planetodetic model should assume pc2pd excludes sidereal
        rotation which is accounted for by the kinematics model. See
        :class:`Planetodetic`.

    The **state** components are:

    ``[0:3] p_x, p_y, p_z``
        translational position (of vehicle center of mass) relative to inertial origin
        expressed in inertial coordinate system
    ``[3:7] q_0, q_1, q_2, q_3``
        quaternion components representing the rotation from the inertial to the
        body-fixed coordinate systems
    ``[7:10] v_x, v_y, v_z``
        translational velocity (of vehicle center of mass) in the inertial coordinate
        system expressed in inertial coordinates
    ``[10:13] omega_X, omega_Y, omega_Z``
        angular velocity from inertial coordinate system to body-fixed coordinate system
        expressed in body-fixed coordinates

    The **input** components are:

    ``[0:3] A_X, A_Y, A_Z``
        translational acceleration of vehicle center of mass due to non-gravitational
        forces in the inertial coordinate system expressed in body-fixed
        Forward-Right-Down (FRD) coordinate system
    ``[3:6] alpha_X, alpha_Y, alpha_Z``
        angular acceleration of vehicle coordinate system in the inertial coordinate
        system expressed in body-fixed FRD coordinates
    ``[6] c_q``
        scaling parameter for non-unit quaternion kinematics integration; this is a
        control for numerical behavior of the integration and can be ignored in most
        circumstances

    The **output** components are:

    ``[0:13] p_x, p_y, p_z, q_0, q_1, q_2, q_3, v_x, v_y, v_z, omega_X, omega_Y, omega_Z``
        The state components listed above
    ``[13:16] lamda_D, phi_D, h``
        translational position of the vehicle center of mass in the planet-fixed frame
        expressed in planetodetic spherical position coordinates: longitude, latitude,
        and altitude
    ``[16:19] psi, theta, phi``
        Euler-angles relating the NED frame to the body-fixed frame: yaw, pitch, roll
    ``[19:22] rho, c_s, mu``
        density, speed of sound, and viscosity outputs of the atmosphere model as a
        function of time and position in planet-fixed frame
    ``[22:25] V_T, alpha, beta``
        true air speed, angle of attack, and angle of sideslip used for aerodynamics
        (includes wind)
    ``[25:28] p_B, q_B, r_B``
        angular velocity components of vehicle-fixed coordinate system relative to NED,
        used for aerodynamic damping derivatives
    ``[28:31] V_N, V_E, V_D``
        velocity of the vehicle center of mass in the planet-fixed frame expressed in
        NED frame
    ``[31:34] W_N, W_E, W_D``
        Output of wind model as a function of time and position in planet-fixed frame

    Additional class attributes:

    ``dim_state, dim_output, dim_input, num_events``
        Used for the SimuPy interface; must be updated for any sub-classes
    ``<variable_name>_idx``
        Provides a named index value for indexing the state or output data columns
    ``output_column_names``
        A list of strings of the output names, useful for constructing data structures
    ``output_column_names_latex``
        A list of strings of the LaTeX output names, useful for labeling plots

    """

    dim_state = 13
    dim_output = 34
    dim_input = 7
    num_events = 0

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

    output_column_names = [
        "p_x",
        "p_y",
        "p_z",
        "q_0",
        "q_1",
        "q_2",
        "q_3",
        "v_x",
        "v_y",
        "v_z",
        "omega_X",
        "omega_Y",
        "omega_Z",
        "lamda_D",
        "phi_D",
        "h_D",
        "psi",
        "theta",
        "phi",
        "rho",
        "c_s",
        "mu",
        "V_T",
        "alpha",
        "beta",
        "p_B",
        "q_B",
        "r_B",
        "V_N",
        "V_E",
        "V_D",
        "W_N",
        "W_E",
        "W_D",
    ]

    output_column_names_latex = [
        "$p_{x}$",
        "$p_{y}$",
        "$p_{z}$",
        "$q_{0}$",
        "$q_{1}$",
        "$q_{2}$",
        "$q_{3}$",
        "$v_{x}$",
        "$v_{y}$",
        "$v_{z}$",
        "$\\omega_{X}$",
        "$\\omega_{Y}$",
        "$\\omega_{Z}$",
        "$\\lambda_{D}$",
        "$\\phi_{D}$",
        "$h_{D}$",
        "$\\psi$",
        "$\\theta$",
        "$\\phi$",
        "$\\rho$",
        "$c_{s}$",
        "$\\mu$",
        "$V_{T}$",
        "$\\alpha$",
        "$\\beta$",
        "$p_{B}$",
        "$q_{B}$",
        "$r_{B}$",
        "$V_{N}$",
        "$V_{E}$",
        "$V_{D}$",
        "$W_{N}$",
        "$W_{E}$",
        "$W_{D}$",
    ]

    def __init__(self, gravity, winds, atmosphere, planetodetics):
        (self.gravity, self.winds, self.atmosphere, self.planetodetics) = (
            gravity,
            winds,
            atmosphere,
            planetodetics,
        )

    def prepare_to_integrate(self, *args, **kwargs):
        """
        SimuPy calls each Block's ``prepare_to_integrate`` function prior to simulation
        Returns the first output vector.
        """
        return self.output_equation_function(*args, **kwargs)

    def state_equation_function(self, t, x, u):
        """
        Compute the kinematic rates used for simulation.
        """
        return kinematics.kinematics_state_function(self, t, *x, *u)

    def output_equation_function(self, t, x):
        """
        Compute the kinematic outputs used for simulation.
        """
        return kinematics.kinematics_output_function(self, t, *x)

    def ic_from_planetodetic(
        self,
        lamda_D=0.0,
        phi_D=0.0,
        h=0.0,
        V_N=0.0,
        V_E=0.0,
        V_D=0.0,
        psi=0.0,
        theta=0.0,
        phi=0.0,
        p_B=0.0,
        q_B=0.0,
        r_B=0.0,
    ):
        """
        A helper function for defining the inertial initial conditions from planetodetic
        definitions

        Parameters
        ----------
        lamda_D, phi_D, h
            planetodetic longitude, latitude, and altitude
        V_N, V_E, V_D
            relative velocity expressed in the North, East, Down (NED) frame
        psi, theta, phi
            Euler-angles relating the NED frame to the body-fixed frame: yaw, pitch,
            roll
        p_B, q_B, r_B
            angular velocity components of vehicle-fixed coordinate system relative to
            NED frame
        """
        return kinematics.ic_from_planetodetic(
            self, lamda_D, phi_D, h, V_N, V_E, V_D, psi, theta, phi, p_B, q_B, r_B
        )

    def local_translational_trim_residual(
        self, p_x, p_y, p_z, q_0, q_1, q_2, q_3, v_x, v_y, v_z, A_X, A_Y, A_Z
    ):
        """A helper function that computes the trim residual in the planet-fixed frame

        Use an optimizer to zero the output of this function to achieve trim. The first
        nine inputs, which are components of the state, might come from an initial
        condition calculation. It is assumed that the trim condition has no relative
        angular velocity. The last three inputs might come from the vehicle dynamics
        model. The control inputs of the vehicle dynamics model is a good candidate for
        the trim-optimizer's free variables.

        Parameters
        ----------
        p_x, p_y, p_z
            translational position (of vehicle center of mass) relative to inertial
            origin expressed in inertial coordinate system
        q_0, q_1, q_2, q_3
            quaternion components representing the rotation from the inertial to the
            body-fixed coordinate systems
        v_x, v_y, v_z
            translational velocity (of vehicle center of mass) in the inertial
            coordinate system expressed in inertial coordinates
        A_X, A_Y, A_Z
            translational acceleration of vehicle center of mass due to
            non-gravitational forces in the inertial coordinate system expressed in
            body-fixed Forward-Right-Down (FRD) coordinate system
        """
        return kinematics.local_translational_trim_residual(
            self, p_x, p_y, p_z, q_0, q_1, q_2, q_3, v_x, v_y, v_z, A_X, A_Y, A_Z
        )

    def inertial_to_NED_dcm(self, t, lamda_D, phi_D):
        """
        Construct a direction cosine matrix (DCM) relating the inertial coordinate
        system to the North-East-Down (NED) coordinate system

        Parameters
        ----------
        t
            time, used to account for sidereal rotation
        lamda_D, phi_D
            planetodetic longitude and latitude
        """
        return kinematics.inertial_to_NED_dcm(self, t, lamda_D, phi_D)


def inertial_to_body_dcm(q_0, q_1, q_2, q_3):
    """
    Construct a direction cosine matrix (DCM) relating the inertial coordinate system
    to the body-fixed Forward-Right-Down (FRD) coordinate system

    Parameters
    ----------
    q_0, q_1, q_2, q_3
        quaternion components defining the DCM
    """
    return kinematics.inertial_to_body_dcm(q_0, q_1, q_2, q_3)


def body_to_NED_dcm(phi, theta, psi):
    """
    Constructs a direction cosine matrix (DCM) relating the body-fixed
    Forward-Right-Down (FRD) coordinate system to the North-East-Down (NED) coordinate
    system

    Parameters
    ----------
    phi, theta, phi
        Euler-angles relating the NED frame to the body-fixed frame: roll, pitch, yaw
    """
    return kinematics.body_to_NED_dcm(phi, theta, psi)


# TODO: verify sign of products of inertia
# TODO: provide pattern to support time varying inertial properties -- should be stateful;
class Vehicle(object):
    """
    The Vehicle is a state-less dynamics block that provides a common calculation for
    the dynamics of a flight vehicle.

    The Vehicle vehicle model is parameterized based on the following components:

    Inertial properties: ``m, I_xx, I_yy, I_zz, I_xy, I_yz, I_xz, x_com, y_com, z_com``
        Inertia is expressed in body-fixed coordinates about the center of mass
        position (``x_com, y_com, z_com``). Center of mass position is relative to an
        arbitrary reference, such as the CAD origin.

    Aerodynamics model:

    ``base_aero_coeffs``:
        assumed to be a function of angles of attack and sideslip, and Mach and Reynolds
        number with signature ``base_aero_coeffs(alpha, beta, Ma, Re)`` and should
        return an array of coefficients in the following order: drag, sideforce, and
        lift (force coefficients expressed in wind coordinate system), static roll,
        pitch, and yaw moments (moment coefficients expressed in body-fixed
        Forward-Right-Down coordinate system), and dynamic or damping roll, pitch, and
        yaw moments (also moment coefficients expressed in body-fixed coordinate
        system). To translate an aerodynamic database that provides forces expressed in
        body-fixed coordinates, the transpose of the ``dynamics.body_to_wind_dcm(alpha,
        beta)`` direction cosine matrix can be used to transform into the wind
        coordinate system. See :func:`get_constant_aero`.
    ``x_mrc, y_mrc, z_mrc``:
        Moment reference center position, defined relative to the same arbitrary
        reference as the center of mass. Moment coefficients are assumed to be about
        this position.
    ``S_A``:
        the reference area
    ``a_l, b_l, c_l``:
        reference lengths for roll, pitch, and yaw
    ``d_l``
        reference length for Reynolds number calculation

    Damping coefficients are transformed with the static moment coefficients, which does
    not usually hold; to model damping coefficients, it is recommended to set the center
    of mass and moment reference center to the same location at the arbitrary reference
    (i.e. (0, 0, 0)).

    Additional, controlled input can be modeled:

    ``input_force_moment, input_force_moment_idx``:
        Non-aerodynamic controlled inputs (e.g. propulsion)
    ``input_aero_coeffs, input_aero_coeffs_idx``:
        Controlled aerodynamics

    Processing of ``input_aero_coeffs`` and ``input_force_moment`` parameterized models
    follow the same logic:

    - if ``None``: assume that a child class will over-write
    - if callable: use directly (so it should have the correct signature)
        - The ``input_aero_coeffs`` callback is assumed to take the flight condition
          arguments of the aerodynamics model, before any additional routed control
          inputs, and return aerodynamic coefficient increments.

        - The ``input_force_moment`` callback takes only routed control inputs and
          returns the translational forces and moments (about the center of mass) being
          modeled, such as for a propulsion model.
    - else: try to build a function that returns constant value(s)

    Attributes ``input_aero_coeffs_idx`` and ``input_force_moment_idx`` are used to
    route extra (control) inputs to the input aero and foce/moment functions
    respectively. Use ``None`` to route all inputs or use a list of integers to route
    particular inputs including an empty list for no routing. ``dim_additional_input``
    as an input argument to ``__init__`` is used to allocate extra control input
    channels in the block diagram.

    The **input** components are:

    ``[0:13] p_x, p_y, p_z, q_0, q_1, q_2, q_3, v_x, v_y, v_z, omega_X, omega_Y, omega_Z``
        The vehicle state components
    ``[13:16] lamda_D, phi_D, h``
        translational position of the vehicle center of mass in the planet-fixed frame
        expressed in planetodetic spherical position coordinates: longitude, latitude,
        and altitude
    ``[16:19] psi, theta, phi``
        Euler-angles relating the NED frame to the body-fixed frame: yaw, pitch, roll
    ``[19:22] rho, c_s, mu``
        density, speed of sound, and viscosity outputs of the atmosphere model as a
        function of time and position in planet-fixed frame
    ``[22:25] V_T, alpha, beta``
        true air speed, angle of attack, and angle of sideslip used for aerodynamics
        (includes wind)
    ``[25:28] p_B, q_B, r_B``
        angular velocity components of body relative to NED used for aerodynamic damping
        derivatives
    ``[28:31] V_N, V_E, V_D``
        velocity of the vehicle center of mass in the planet-fixed frame expressed in
        NED frame
    ``[31:34] W_N, W_E, W_D``
        winds acting on the vehicle expressed in NED frame

    The **output** components are:

    ``[0:3] A_X, A_Y, A_Z``
        translational acceleration of vehicle center of mass due to non-gravitational
        forces in the inertial coordinate system expressed in body-fixed
        Forward-Right-Down (FRD) coordinates, computed assuming constant mass and total
        non-gravitational forces due to aerodynamics (base model and extra aerodynamic
        coefficients input components) and the extra force input components.
    ``[3:6] alpha_X, alpha_Y, alpha_Z``
        angular acceleration of vehicle coordinate system in the inertial coordinate
        system expressed in body-fixed FRD coordinates, about the center of mass
        computed assuming constant inertia and total moments due to aerodynamics (base
        model and extra aerodynamic coefficients input components) and the extra moment
        input components.

    Additional class attributes:

    ``dim_state, dim_output, dim_input, num_events``
        Used for the SimuPy interface; must be updated for any sub-classes
    ``<variable_name>_idx``
        Provides a named index value for indexing the state or output data columns
    ``output_column_names``
        A list of strings of the output names, useful for constructing data structures
    ``output_column_names_latex``
        A list of strings of the LaTeX output names, useful for labeling plots
    """

    dim_state = 0
    dim_output = 6
    dim_input = 34
    num_events = 0

    initial_condition = np.empty(0)

    output_column_names = ["A_X", "A_Y", "A_Z", "alpha_X", "alpha_Y", "alpha_Z"]

    output_column_names_latex = [
        "$A_{X}$",
        "$A_{Y}$",
        "$A_{Z}$",
        "$\\alpha_{X}$",
        "$\\alpha_{Y}$",
        "$\\alpha_{Z}$",
    ]

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

    def __init__(
        self,
        m,
        I_xx,
        I_yy,
        I_zz,
        I_xy,
        I_yz,
        I_xz,
        x_com,
        y_com,
        z_com,
        base_aero_coeffs,
        x_mrc,
        y_mrc,
        z_mrc,
        S_A,
        a_l,
        b_l,
        c_l,
        d_l,
        input_aero_coeffs=0.0,
        input_force_moment=0.0,
        input_aero_coeffs_idx=None,
        input_force_moment_idx=None,
        dim_additional_input=0,
    ):

        self.dim_input = Vehicle.dim_input + dim_additional_input

        if callable(base_aero_coeffs):
            pass
        else:
            base_aero_coeff_vals = np.array(base_aero_coeffs, dtype=np.float_).reshape(
                -1
            )
            if base_aero_coeff_vals.size == 1:
                base_aero_coeff_vals = np.tile(base_aero_coeff_vals, 9)
            base_aero_coeffs = get_constant_aero(*base_aero_coeff_vals)
        self.base_aero_coeffs = base_aero_coeffs

        if callable(input_aero_coeffs):
            pass
        else:
            input_aero_coeff_vals = np.array(
                input_aero_coeffs, dtype=np.float_
            ).reshape(-1)
            if input_aero_coeff_vals.size == 1:
                input_aero_coeff_vals = np.tile(input_aero_coeff_vals, 9)
            input_aero_coeffs = get_constant_aero(*input_aero_coeff_vals)
        self.input_aero_coeffs = input_aero_coeffs

        if callable(input_force_moment):
            pass
        else:
            input_force_moment_vals = np.array(
                input_force_moment, dtype=np.float_
            ).reshape(-1)
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

        (
            self.m,
            self.I_xx,
            self.I_yy,
            self.I_zz,
            self.I_xy,
            self.I_yz,
            self.I_xz,
            self.x_com,
            self.y_com,
            self.z_com,
            self.x_mrc,
            self.y_mrc,
            self.z_mrc,
            self.S_A,
            self.a_l,
            self.b_l,
            self.c_l,
            self.d_l,
        ) = (
            m,
            I_xx,
            I_yy,
            I_zz,
            I_xy,
            I_yz,
            I_xz,
            x_com,
            y_com,
            z_com,
            x_mrc,
            y_mrc,
            z_mrc,
            S_A,
            a_l,
            b_l,
            c_l,
            d_l,
        )

    def prepare_to_integrate(self, *args, **kwargs):
        """
        SimuPy calls each Block's ``prepare_to_integrate`` function prior to simulation
        Returns the first output vector.
        """
        return self.output_equation_function(*args, **kwargs)

    def _input_force_moment(
        self,
        t,
        p_x,
        p_y,
        p_z,
        q_0,
        q_1,
        q_2,
        q_3,
        v_x,
        v_y,
        v_z,
        omega_X,
        omega_Y,
        omega_Z,
        lamda_D,
        phi_D,
        h_D,
        psi,
        theta,
        phi,
        rho,
        c_s,
        mu,
        V_T,
        alpha,
        beta,
        p_B,
        q_B,
        r_B,
        V_N,
        V_E,
        V_D,
        W_N,
        W_E,
        W_D,
        qbar,
        Ma,
        Re,
        *args
    ):
        filtered_args = np.array(args)[self.input_force_moment_idx]
        return self.input_force_moment(
            t,
            p_x,
            p_y,
            p_z,
            v_x,
            v_y,
            v_z,
            q_0,
            q_1,
            q_2,
            q_3,
            omega_X,
            omega_Y,
            omega_Z,
            lamda_D,
            phi_D,
            h_D,
            psi,
            theta,
            phi,
            rho,
            c_s,
            mu,
            V_T,
            alpha,
            beta,
            p_B,
            q_B,
            r_B,
            V_N,
            V_E,
            V_D,
            W_N,
            W_E,
            W_D,
            qbar,
            Ma,
            Re,
            *filtered_args
        )

    def output_equation_function(self, t, u):
        """
        Compute the dynamic outputs used for simulation.
        """
        uu = u[..., : Vehicle.dim_input]
        u_extra = u[..., Vehicle.dim_input :]
        return dynamics.dynamics_output_function(self, t, *uu, *u_extra)

    def _tot_aero_forces_moments(
        self, qbar, Ma, Re, V_T, alpha, beta, p_B, q_B, r_B, *args
    ):
        filtered_args = np.array(args)[self.input_aero_coeffs_idx]
        return self.tot_aero_forces_moments(
            qbar, Ma, Re, V_T, alpha, beta, p_B, q_B, r_B, *filtered_args
        )

    def tot_aero_forces_moments(
        self, qbar, Ma, Re, V_T, alpha, beta, p_B, q_B, r_B, *args
    ):
        """
        Helper function to compute the total aerodynamic forces and moments for a given
        flight  and vehicle condition. Used to evaluate the aerodynamic forces and
        moments in calculating the dynamics output.

        Parameters
        ----------
        qbar
            dynamic pressure in Pascals
        Mach
            Mach number (unitless)
        Re
            Reynolds number (unitless)
        V_T, alpha, beta
            true air speed, angle of attack, and angle of sideslip used for aerodynamics
            (includes wind)
        p_B, q_B, r_B
            angular velocity components of body relative to NED used for aerodynamic
            damping derivatives
        *args
            Captures any additional control inputs for the controlled aerodynamics model
        """
        return dynamics.tot_aero_forces_moments(
            self, qbar, Ma, Re, V_T, alpha, beta, p_B, q_B, r_B, *args
        )

    def mrc_to_com_cpm(self):
        """
        Construct a skew symmetric matrix that can perform the right-sided cross product
        of the position vector from the moment reference center to the center of mass,
        used for accounting for the translational force contribution to the moment.
        """
        return dynamics.mrc_to_com_cpm(self)


def body_to_wind_dcm(alpha, beta):
    """
    Construct a direction cosine matrix (DCM) relating the body-fixed Forward-Right-Down
    (FRD) coordinate system to the wind coordinate system.

    Parameters
    ----------
    alpha
        angle of attack
    beta
        angle of sideslip
    """
    return dynamics.body_to_wind_dcm(alpha, beta)
