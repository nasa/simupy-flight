Search.setIndex({"docnames": ["api", "index", "nesc_test_cases/index", "nesc_test_cases/nesc_case01", "nesc_test_cases/nesc_case02", "nesc_test_cases/nesc_case03", "nesc_test_cases/nesc_case04", "nesc_test_cases/nesc_case05", "nesc_test_cases/nesc_case06", "nesc_test_cases/nesc_case07", "nesc_test_cases/nesc_case08", "nesc_test_cases/nesc_case09", "nesc_test_cases/nesc_case10", "nesc_test_cases/nesc_case11", "nesc_test_cases/nesc_case13p1", "nesc_test_cases/nesc_case13p2", "nesc_test_cases/nesc_case13p3", "nesc_test_cases/nesc_case13p4", "nesc_test_cases/nesc_case15", "nesc_test_cases/nesc_case16", "nesc_test_cases/sg_execution_times"], "filenames": ["api.rst", "index.rst", "nesc_test_cases/index.rst", "nesc_test_cases/nesc_case01.rst", "nesc_test_cases/nesc_case02.rst", "nesc_test_cases/nesc_case03.rst", "nesc_test_cases/nesc_case04.rst", "nesc_test_cases/nesc_case05.rst", "nesc_test_cases/nesc_case06.rst", "nesc_test_cases/nesc_case07.rst", "nesc_test_cases/nesc_case08.rst", "nesc_test_cases/nesc_case09.rst", "nesc_test_cases/nesc_case10.rst", "nesc_test_cases/nesc_case11.rst", "nesc_test_cases/nesc_case13p1.rst", "nesc_test_cases/nesc_case13p2.rst", "nesc_test_cases/nesc_case13p3.rst", "nesc_test_cases/nesc_case13p4.rst", "nesc_test_cases/nesc_case15.rst", "nesc_test_cases/nesc_case16.rst", "nesc_test_cases/sg_execution_times.rst"], "titles": ["SimuPy Flight API", "SimuPy Flight Vehicle Toolkit", "NESC Test Cases", "Case 1: Dropped sphere with no drag", "Case 2: Tumbling brick with no damping or drag", "Case 3: Tumbling brick with dynamic damping, no drag", "Case 4: Dropped sphere over non-rotating, spherical Earth", "Case 5: Dropped sphere over rotating, spherical Earth", "Case 6: Dropped sphere over rotating, ellipsoidal Earth", "Case 7: Sphere dropped through a steady wind field", "Case 8: Sphere dropped through a varying wind field", "Case 9: Sphere launched ballistically eastward along Equator", "Case 10: Sphere launched ballistically northward along Prime Meridian", "Case 11: Subsonic F-16 trimmed flight across earth", "Case 13.1: Altitude change of a subsonic aircraft", "Case 13.2: Velocity change of a subsonic aircraft", "Case 13.3: Course change of a subsonic aircraft", "Case 13.4: Lateral offset manuever of a subsonic aircraft", "Case 15: Circular F-16 flight around North Pole", "Case 16: Circular F-16 flight around Equator/dateline intersection", "Computation times"], "terms": {"class": [0, 1, 13], "simupy_flight": [0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "planet": [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "graviti": [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "wind": [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], "atmospher": [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "planetodet": [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "The": [0, 1, 3, 4, 13], "state": [0, 1], "dynam": [0, 1, 2, 3, 4, 6, 7, 8, 13, 20], "block": [0, 1, 13], "provid": [0, 1], "kinemat": [0, 1, 3, 4, 5, 6, 7, 8, 13], "equat": [0, 1, 2, 20], "motion": [0, 1], "integr": [0, 1, 6, 13], "an": [0, 1], "output": [0, 1, 13, 17, 18, 19], "commonli": 0, "us": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "variabl": [0, 1, 13], "vehicl": [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "accord": [0, 3], "model": [0, 1, 3, 4, 5, 13], "parameter": 0, "base": [0, 1, 13], "follow": [0, 1], "compon": [0, 1], "translat": [0, 1, 3, 11], "acceler": [0, 1, 13, 14, 18, 19], "due": [0, 1, 4, 9], "function": [0, 1, 13, 14, 18, 19], "fix": [0, 1, 6], "posit": [0, 1, 13], "rectangular": 0, "coordin": [0, 3], "see": [0, 1], "earth_j2_grav": [0, 3, 4, 5, 8, 9, 10, 11, 12, 13], "exampl": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "local": 0, "ned": 0, "frame": 0, "time": [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "indic": 0, "specifi": 0, "direct": 0, "so": [0, 1, 3], "from": [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "west": [0, 9], "w_e": 0, "get_constant_wind": [0, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13], "gener": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "simpl": [0, 14, 15, 16, 17], "densiti": [0, 13], "speed": [0, 13], "sound": 0, "viscos": 0, "i": [0, 18], "e": 0, "stochast": 0, "express": [0, 1], "geodet": [0, 18], "longitud": [0, 13, 18, 19], "latitud": [0, 13, 18, 19], "altitud": [0, 2, 10, 20], "get_constant_atmospher": [0, 3, 4], "atmosphere_1976": [0, 5, 6, 7, 8, 9, 10, 11, 12, 13], "standard": [0, 1], "1976": [0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "must": [0, 1], "rotat": [0, 2, 3, 4, 5, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], "rate": [0, 1, 4, 11, 12, 13], "rad": 0, "s": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "omega_p": [0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "pd2pcf": 0, "pcf2pd": 0, "convert": 0, "between": [0, 1], "planetocentr": 0, "spheric": [0, 2, 20], "futur": [0, 1], "version": [0, 1], "mai": [0, 1], "support": 0, "nutat": 0, "precess": 0, "should": [0, 1], "assum": 0, "pc2pd": 0, "exclud": 0, "sider": 0, "which": [0, 1, 13], "account": 0, "ar": [0, 1, 2, 13], "0": [0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], "3": [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 20], "p_x": 0, "p_y": 0, "p_z": 0, "center": [0, 1, 2], "mass": 0, "rel": [0, 4, 11, 12], "inerti": [0, 1, 3, 4, 5], "origin": [0, 1], "system": [0, 1, 3, 13, 14, 15, 16, 17, 18, 19], "7": [0, 2, 3, 8, 10, 11, 12, 20], "q_0": 0, "q_1": 0, "q_2": 0, "q_3": 0, "quaternion": 0, "repres": [0, 1], "bodi": [0, 1], "10": [0, 2, 4, 5, 6, 7, 13, 14, 18, 20], "v_x": 0, "v_y": 0, "v_z": 0, "veloc": [0, 1, 2, 3, 11, 12, 20], "13": [0, 2, 13, 18, 19, 20], "omega_x": 0, "omega_i": 0, "omega_z": 0, "angular": [0, 1, 4, 11, 12], "input": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18, 19], "a_x": 0, "a_i": 0, "a_z": 0, "non": [0, 2, 20], "gravit": [0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "forc": [0, 4], "forward": 0, "right": [0, 1, 16], "down": 0, "frd": 0, "6": [0, 1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 18, 19, 20], "alpha_x": 0, "alpha_i": 0, "alpha_z": 0, "c_q": 0, "scale": 0, "paramet": 0, "unit": [0, 1], "thi": [0, 1, 3, 13], "control": [0, 1, 13, 14, 15, 16, 17], "numer": 0, "behavior": 0, "can": [0, 1, 2], "ignor": 0, "most": 0, "circumst": 0, "list": [0, 3], "abov": [0, 3, 13], "16": [0, 1, 2, 14, 15, 16, 17, 20], "lamda_d": [0, 12, 13, 18, 19], "phi_d": [0, 12, 13, 18, 19], "h": [0, 1, 12, 13, 14, 18, 19], "19": 0, "psi": [0, 12, 13, 16, 18, 19], "theta": [0, 12, 13, 18, 19], "phi": [0, 12, 13, 18, 19], "euler": [0, 13], "angl": [0, 13, 18], "relat": 0, "yaw": 0, "pitch": [0, 13, 14, 18, 19], "roll": [0, 13, 14, 18, 19], "22": [0, 4, 5, 6, 7, 8, 20], "rho": [0, 13, 18], "c_": 0, "mu": 0, "25": [0, 18], "v_t": [0, 4, 5, 6, 7, 8, 13, 18], "alpha": [0, 13, 18, 19], "beta": [0, 13, 18], "true": [0, 13, 14, 15, 16, 17, 18, 19], "air": [0, 1, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19], "attack": 0, "sideslip": 0, "aerodynam": [0, 1, 3, 13], "includ": [0, 1], "28": 0, "p_b": [0, 12, 13, 18, 19], "q_b": [0, 12, 13, 18, 19], "r_b": [0, 12, 13, 18, 19], "damp": [0, 2, 20], "deriv": [0, 1], "31": [0, 15, 20], "v_n": [0, 12, 13, 18, 19], "v_e": [0, 12, 13, 18, 19], "v_d": [0, 12, 13, 18, 19], "34": 0, "w_n": 0, "w_d": 0, "addit": [0, 1, 2, 13], "attribut": 0, "dim_stat": [0, 13, 17], "dim_output": [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18, 19], "dim_input": [0, 13, 17], "num_ev": [0, 13], "interfac": [0, 1], "updat": 0, "ani": [0, 1], "sub": [0, 13], "variable_nam": 0, "_idx": 0, "name": [0, 1], "index": 0, "valu": [0, 3, 4, 5, 6, 7, 8, 13, 14, 18, 19], "data": [0, 1, 13], "column": 0, "output_column_nam": 0, "A": [0, 1], "string": 0, "construct": [0, 1, 3], "structur": 0, "output_column_names_latex": 0, "latex": 0, "label": [0, 18, 19], "plot": [0, 1, 3, 18, 19], "ic_from_planetodet": [0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18, 19], "helper": [0, 13], "defin": [0, 3], "initi": [0, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19], "condit": [0, 1, 3, 13], "definit": [0, 1], "north": [0, 2, 12, 20], "east": [0, 11], "inertial_to_ned_dcm": 0, "t": [0, 1, 10, 13, 14, 15, 16, 17, 18], "cosin": 0, "matrix": 0, "dcm": 0, "local_translational_trim_residu": [0, 13], "comput": [0, 1, 13], "trim": [0, 1, 2, 20], "residu": [0, 13], "optim": [0, 13, 14, 18, 19], "zero": [0, 1, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13], "achiev": [0, 13], "first": [0, 13], "nine": 0, "might": 0, "come": [0, 1], "calcul": [0, 13, 18], "It": [0, 1, 13], "ha": [0, 1, 13], "last": 0, "three": 0, "good": 0, "candid": 0, "free": [0, 1, 13], "output_equation_funct": [0, 13, 17, 18, 19], "x": [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18, 19], "simul": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "prepare_to_integr": [0, 13], "arg": [0, 13], "kwarg": 0, "call": [0, 1], "each": [0, 1], "prior": [0, 1], "return": [0, 10, 13, 17, 18], "vector": 0, "state_equation_funct": [0, 13, 17], "u": [0, 13, 18], "f": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 20], "namespac": [0, 1], "float": 0, "equatori": 0, "radiu": 0, "flatten": 0, "planetari": 0, "flat": 0, "m": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "i_xx": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "i_yi": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "i_zz": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "i_xi": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "i_yz": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "i_xz": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "x_com": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "y_com": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "z_com": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "base_aero_coeff": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "x_mrc": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "y_mrc": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "z_mrc": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "s_a": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "a_l": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "b_l": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "c_l": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "d_l": [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], "input_aero_coeff": [0, 1], "input_force_mo": [0, 1], "input_aero_coeffs_idx": [0, 1], "none": [0, 18, 19], "input_force_moment_idx": [0, 1], "dim_additional_input": 0, "less": [0, 1], "common": [0, 1], "properti": 0, "inertia": [0, 1, 13], "about": [0, 1, 2], "arbitrari": 0, "refer": [0, 1], "cad": 0, "mach": 0, "reynold": 0, "number": [0, 1], "signatur": 0, "ma": 0, "re": [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "arrai": [0, 13, 18], "coeffici": [0, 3, 4, 5], "order": 0, "drag": [0, 2, 20], "sideforc": 0, "lift": 0, "static": 0, "moment": 0, "also": [0, 1], "To": [0, 1, 13], "databas": 0, "transpos": 0, "body_to_wind_dcm": 0, "transform": [0, 3, 18], "get_constant_aero": [0, 4, 5, 6, 7, 8, 9, 10, 11, 12], "same": [0, 1], "area": 0, "length": 0, "doe": [0, 1], "usual": 0, "hold": [0, 1], "recommend": [0, 1], "set": [0, 1, 3, 4, 5, 13], "locat": [0, 1], "g": 0, "propuls": [0, 13], "process": [0, 1], "logic": 0, "child": 0, "over": [0, 2, 20], "write": [0, 1], "callabl": 0, "directli": [0, 1], "have": [0, 1], "correct": 0, "callback": [0, 1], "take": 0, "argument": [0, 13], "befor": [0, 1], "rout": [0, 1], "increment": 0, "onli": 0, "being": [0, 1], "els": [0, 13, 17, 18, 19], "try": 0, "build": [0, 13], "constant": [0, 5, 6, 7, 8, 9, 10, 11, 12], "extra": [0, 1], "aero": [0, 1, 5], "foce": 0, "respect": [0, 1], "all": [0, 1, 2], "integ": 0, "particular": [0, 1], "empti": 0, "__init__": [0, 13], "alloc": 0, "channel": [0, 13], "diagram": [0, 13], "act": 0, "total": [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], "mrc_to_com_cpm": 0, "skew": 0, "symmetr": 0, "perform": [0, 1, 3, 13], "side": 0, "cross": 0, "product": [0, 1], "contribut": 0, "tot_aero_forces_mo": 0, "qbar": 0, "given": [0, 13], "evalu": [0, 13, 14, 18, 19], "pressur": 0, "pascal": 0, "unitless": 0, "captur": 0, "body_to_ned_dcm": 0, "px": 0, "py": [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], "pz": 0, "j2": [0, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "cd_b": [0, 6, 7, 8, 9, 10, 11, 12], "cs_b": 0, "cl_b": 0, "clcal_b": 0, "cmcal_b": 0, "cncal_b": 0, "cp_b": [0, 5], "cq_b": [0, 5], "cr_b": [0, 5], "regardless": 0, "density_v": 0, "speed_of_sound_v": 0, "1": [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20], "viscocity_v": 0, "get_constant_force_mo": 0, "fx": [0, 10], "fy": [0, 10], "fz": [0, 10], "mx": 0, "my": 0, "mz": 0, "wx": 0, "wy": [0, 9], "wz": 0, "get_flat_pc2pd": 0, "pi": [0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 18, 19], "2": [0, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20], "type": 0, "get_flat_pd2pc": 0, "90": [0, 11, 18, 19], "degre": [0, 1, 2, 11, 12, 13, 16], "get_nonflat_pc2pd": 0, "ellipsoid": [0, 2, 20], "meter": 0, "consist": 0, "other": [0, 1], "erfa": 0, "implement": [0, 1, 2, 13], "fukushima": 0, "method": [0, 3, 13], "get_nonflat_pd2pc": 0, "get_spherical_grav": [0, 6, 7], "gravitational_const": 0, "inertial_to_body_dcm": 0, "import": [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "part": 1, "innov": 1, "aerospac": 1, "technolog": 1, "nasa": [1, 2], "engin": [1, 2], "safeti": [1, 2], "identifi": 1, "address": 1, "need": [1, 3, 4], "verifi": [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "through": [1, 2, 13, 20], "work": [1, 3, 4, 5, 6, 7, 8], "six": [1, 2], "freedom": [1, 2], "dof": [1, 2], "author": 1, "wa": 1, "prompt": 1, "develop": 1, "tool": 1, "would": 1, "allow": 1, "rapid": 1, "novel": 1, "concept": 1, "hyperson": 1, "entri": 1, "urban": 1, "mobil": 1, "softwar": 1, "leverag": 1, "open": 1, "sourc": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "scientif": 1, "effici": 1, "framework": 1, "python": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "compos": [1, 13], "altern": 1, "simulink": 1, "solv": 1, "result": [1, 3, 13], "differenti": 1, "scipi": [1, 13, 14, 15, 16, 17, 18, 19], "wrapper": 1, "fortran": 1, "orient": 1, "correspond": 1, "sympi": 1, "symbol": 1, "code": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "well": 1, "solut": [1, 13], "invers": 1, "geodesi": [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "problem": 1, "dynamicalsystem": [1, 13, 17], "api": 1, "priorit": 1, "explicit": 1, "modular": 1, "facilit": 1, "adher": 1, "best": 1, "practic": 1, "commun": 1, "principl": 1, "separ": 1, "reflect": 1, "distinct": 1, "two": 1, "could": 1, "been": 1, "wai": 1, "inde": 1, "within": 1, "howev": [1, 3], "often": 1, "conveni": 1, "break": 1, "when": 1, "intermedi": 1, "signal": [1, 13, 14, 15, 16, 17], "multipl": 1, "purpos": 1, "In": [1, 13], "easi": 1, "access": 1, "sens": 1, "measur": 1, "particularli": 1, "disciplin": 1, "approach": 1, "out": [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "where": 1, "depend": 1, "whose": 1, "chang": [1, 2, 19, 20], "without": 1, "blockdiagram": [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "manipul": 1, "like": 1, "connect": [1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18, 19], "guidanc": 1, "navig": 1, "reduc": 1, "learn": 1, "curv": 1, "new": 1, "user": 1, "dim_extra_input": 1, "drive": [1, 13], "For": [1, 14, 15, 16, 17], "more": 1, "detail": 1, "docstr": 1, "via": 1, "help": 1, "clone": 1, "repositori": 1, "pip": 1, "git": 1, "http": 1, "github": 1, "com": 1, "cd": 1, "These": [1, 2], "nesc_test_cas": [1, 20], "directori": 1, "view": 1, "galleri": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "document": 1, "inform": [1, 2], "found": [1, 2], "appendic": 1, "report": 1, "run": [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "requir": [1, 13], "few": 1, "thei": 1, "simpli": 1, "execut": [1, 20], "nesc_cas": 1, "file": [1, 20], "run_nesc_cas": 1, "script": [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "iter": [1, 13, 14, 18, 19], "load": 1, "nesc_data": 1, "them": [1, 13], "along": [1, 2, 20], "obtain": 1, "here": [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "divid": [1, 3, 4, 5, 6, 7, 8], "runtimewarn": [1, 3, 4, 5, 6, 7, 8], "check": [1, 3], "handl": [1, 3], "correctli": [1, 3], "By": 1, "default": [1, 4, 5, 6, 7, 8], "regress": 1, "against": 1, "pass": [1, 13], "flag": 1, "option": [1, 13], "everi": 1, "annot": 1, "least": 1, "basic": 1, "descript": 1, "adapt": [1, 13], "moder": 1, "highlight": 1, "usag": 1, "11": [1, 2, 14, 15, 16, 17, 19, 20], "demonstr": 1, "straight": [1, 14, 15, 16, 17, 18, 19], "level": [1, 13, 14, 15, 16, 17, 18, 19], "thoroughli": 1, "illustr": 1, "how": 1, "sophist": 1, "itself": 1, "becaus": 1, "one": 1, "american": 1, "institut": 1, "aeronaut": 1, "astronaut": 1, "aiaa": 1, "xml": 1, "exchang": 1, "format": 1, "aircraft": [1, 2, 20], "markup": 1, "languag": 1, "dave": 1, "ml": 1, "parse_daveml": 1, "submodul": 1, "parser": 1, "valid": 1, "processdaveml": 1, "filenam": 1, "creat": [1, 13], "replac": [1, 14, 15, 16, 17], "extens": 1, "featur": 1, "process_nesc_daveml": 1, "specif": [1, 4], "element": 1, "assist": 1, "verif": 1, "debug": 1, "add": 1, "run_check": 1, "f16": [1, 13, 18, 19], "themselv": 1, "f16_aero": 1, "pleas": 1, "feel": 1, "share": 1, "thought": 1, "opinion": 1, "issu": 1, "feedback": [1, 13, 18], "welcom": 1, "appreci": 1, "bug": 1, "pull": 1, "request": 1, "alwai": 1, "etc": 1, "discuss": 1, "If": 1, "isn": 1, "your": 1, "rational": 1, "spend": 1, "signific": 1, "prepar": 1, "ideal": 1, "you": 1, "perfectli": 1, "fine": 1, "progress": 1, "review": 1, "accept": 1, "contributor": 1, "agreement": 1, "we": [1, 3, 13], "dure": 1, "releas": 1, "under": 1, "copyright": 1, "2021": 1, "govern": 1, "administr": 1, "nation": 1, "space": 1, "reserv": 1, "No": 1, "warranti": 1, "THE": 1, "subject": 1, "IS": 1, "AS": 1, "OF": 1, "kind": 1, "either": 1, "impli": 1, "OR": 1, "statutori": 1, "BUT": 1, "NOT": 1, "limit": 1, "TO": 1, "THAT": 1, "WILL": 1, "conform": 1, "merchant": 1, "fit": 1, "FOR": 1, "infring": 1, "BE": 1, "error": [1, 13], "IF": 1, "IN": 1, "manner": 1, "constitut": 1, "endors": 1, "BY": 1, "agenc": 1, "recipi": 1, "hardwar": 1, "applic": 1, "further": 1, "AND": 1, "liabil": 1, "regard": 1, "third": 1, "parti": 1, "present": 1, "distribut": 1, "IT": 1, "waiver": 1, "indemn": 1, "agre": 1, "waiv": 1, "claim": 1, "ITS": 1, "contractor": 1, "subcontractor": 1, "demand": 1, "damag": 1, "expens": 1, "loss": 1, "aris": 1, "SUCH": 1, "ON": 1, "shall": 1, "indemnifi": 1, "harmless": 1, "extent": 1, "permit": 1, "law": 1, "sole": 1, "remedi": 1, "matter": 1, "immedi": 1, "unilater": 1, "termin": [1, 13, 14, 18, 19], "simupi": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "flight": [2, 3, 4, 5, 6, 7, 8, 20], "toolkit": 2, "drop": [2, 20], "sphere": [2, 20], "tumbl": [2, 20], "brick": [2, 20], "4": [2, 13, 14, 15, 18, 19, 20], "earth": [2, 17, 18, 19, 20], "5": [2, 14, 15, 18, 19, 20], "steadi": [2, 13, 20], "field": [2, 20], "8": [2, 13, 19, 20], "vari": [2, 20], "9": [2, 19, 20], "launch": [2, 20], "ballist": [2, 20], "eastward": [2, 20], "northward": [2, 20], "prime": [2, 20], "meridian": [2, 20], "subson": [2, 20], "across": [2, 20], "cours": [2, 13, 17, 20], "later": [2, 13, 20], "offset": [2, 13, 20], "manuev": [2, 14, 15, 16, 20], "15": [2, 13, 16, 20], "circular": [2, 20], "around": [2, 20], "pole": [2, 20], "datelin": [2, 20], "intersect": [2, 20], "download": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "nesc_test_cases_python": 2, "zip": 2, "jupyt": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "notebook": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "nesc_test_cases_jupyt": 2, "sphinx": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "click": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "full": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "eom": [3, 4, 11], "wg": [3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "84": [3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "std": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "still": [3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19], "dragless": [3, 4, 5], "note": [3, 4, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19], "block_diagram": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18, 19], "numpi": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 17, 18, 19], "np": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 17, 18, 19], "nesc_testcase_help": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "plot_nesc_comparison": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "int_opt": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "benchmark": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "ft_per_m": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19], "configur": [3, 5, 13, 14, 15, 16, 17], "earth_equitorial_radiu": [3, 4, 5, 8, 9, 10, 11, 12, 13], "earth_rotation_r": [3, 4, 5, 7, 8, 9, 10, 11, 12, 13], "earth_f": [3, 4, 5, 8, 9, 10, 11, 12, 13], "sinc": 3, "just": 3, "object": 3, "bd": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "lat_ic": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "180": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 18, 19], "long_ic": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "h_ic": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "30_000": [3, 4, 5, 6, 7, 8, 9, 10], "v_n_ic": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "v_e_ic": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "v_d_ic": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "psi_ic": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "theta_": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "phi_ic": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "omega_x_": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "omega_y_": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "omega_z_": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "initial_condit": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18, 19], "b": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "30": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17, 18, 20], "integrator_opt": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "home": [3, 4, 5, 6, 7, 8], "runner": [3, 4, 5, 6, 7, 8], "198": [3, 4, 5, 6, 7, 8], "invalid": [3, 4, 5, 6, 7, 8], "encount": [3, 4, 5, 6, 7, 8], "double_scalar": [3, 4, 5, 6, 7, 8], "x60": [3, 4, 5, 6, 7, 8], "x6": [3, 4, 5, 6, 7, 8], "x26": [3, 4, 5, 6, 7, 8], "x54": [3, 4, 5, 6, 7, 8], "x58": [3, 4, 5, 6, 7, 8], "x13": [3, 4, 5, 6, 7, 8], "x14": [3, 4, 5, 6, 7, 8], "x59": [3, 4, 5, 6, 7, 8], "x33": [3, 4, 5, 6, 7, 8], "x55": [3, 4, 5, 6, 7, 8], "x53": [3, 4, 5, 6, 7, 8], "321": 3, "01": [3, 20], "minut": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "081": [3, 20], "second": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "nesc_case01": [3, 20], "ipynb": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "kg_per_slug": [4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "ixx": [4, 5, 6, 7, 8, 9, 10, 11, 12], "001894220": [4, 5], "slug": [4, 5, 6, 7, 8, 9, 10, 11, 12], "ft2": [4, 5, 6, 7, 8, 9, 10, 11, 12], "iyi": [4, 5, 6, 7, 8, 9, 10, 11, 12], "006211019": [4, 5], "izz": [4, 5, 6, 7, 8, 9, 10, 11, 12], "007194665": [4, 5], "ixi": [4, 5, 6, 7, 8, 9, 10, 11, 12], "iyz": [4, 5, 6, 7, 8, 9, 10, 11, 12], "izx": [4, 5, 6, 7, 8, 9, 10, 11, 12], "155404754": [4, 5], "y": [4, 5, 6, 7, 8, 9, 10, 11, 12, 18, 19], "z": [4, 5, 6, 7, 8, 9, 10, 11, 12], "arang": [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18, 19], "20": [4, 5, 6, 7, 9, 10, 14, 15, 17], "02": [4, 18], "x15": [4, 5, 6, 7, 8], "select": [4, 5, 6, 7, 8, 13], "greater": [4, 5, 6, 7, 8], "955": 4, "830": [4, 20], "nesc_case02": [4, 20], "coupl": 5, "22222": 5, "aero_model": 5, "03": [5, 13, 14], "027": 5, "799": [5, 20], "nesc_case03": [5, 20], "r": [6, 7], "round": [6, 7], "c_d": [6, 7, 8, 9, 10, 11, 12], "earth_spherical_gravity_const": [6, 7], "20902255": [6, 7], "199": [6, 7], "1963495": [6, 7, 8, 9, 10, 11, 12], "008": 6, "04": [6, 18], "692": [6, 20], "nesc_case04": [6, 20], "011": [7, 10, 20], "05": 7, "685": [7, 20], "nesc_case05": [7, 20], "030": 8, "06": [8, 18, 20], "020": [8, 20], "nesc_case06": [8, 20], "effect": 9, "ft": [9, 11, 12, 13, 14, 17], "025": 9, "07": [9, 18, 20], "077": [9, 20], "nesc_case07": [9, 20], "2d": 10, "shear": 10, "linearli": 10, "wy_1": 10, "70": 10, "x_1": 10, "wy_2": 10, "x_2": 10, "def": [10, 13, 17, 18], "linear_wind": 10, "winds_out": 10, "039": 10, "08": [10, 20], "nesc_case08": [10, 20], "sqrt": [11, 12, 13, 18, 19], "2000": [11, 12, 17], "align": [11, 12], "45": [11, 12, 13, 18], "vertic": [11, 12], "head": [11, 12, 13, 16], "platform": [11, 12], "000": [11, 12, 13], "1000": [11, 12], "004178073": [11, 12], "034": 11, "09": [11, 13, 14, 18, 19], "016": [11, 20], "nesc_case09": [11, 20], "corioli": 12, "p_b_ic": 12, "q_b_ic": 12, "r_b_ic": 12, "026": 12, "043": [12, 20], "nesc_case10": [12, 20], "unaug": 13, "kffa": 13, "airport": 13, "335": 13, "kta": 13, "stabil": 13, "augment": 13, "off": 13, "test": 13, "find": 13, "sonic": 13, "sight": 13, "modul": 13, "minim": 13, "f16_model": [13, 18, 19], "f16_control": 13, "f16_vehicl": [13, 18, 19], "dictionari": 13, "keyword": 13, "spec_ic_arg": [13, 14, 16, 17, 18, 19], "dict": [13, 18, 19], "36": 13, "01916667": 13, "75": 13, "67444444": 13, "10_013": 13, "400": 13, "653814": 13, "knots_per_mp": [13, 18, 19], "94384": 13, "surfac": 13, "equival": [13, 18], "rho_0": [13, 18, 19], "controller_feedback_indic": 13, "h_d_idx": [13, 18], "v_t_idx": [13, 18, 19], "alpha_idx": [13, 18], "beta_idx": [13, 18], "psi_idx": [13, 18], "theta_idx": [13, 18], "phi_idx": [13, 18], "p_b_idx": [13, 18], "q_b_idx": [13, 18], "r_b_idx": [13, 18], "rho_idx": [13, 18, 19], "dim_feedback": [13, 14, 15, 16, 17, 18, 19], "len": [13, 18], "auto": [13, 14, 15, 16, 17, 18, 19], "pilot": [13, 14, 15, 16, 17, 18, 19], "daveml": 13, "sever": 13, "f16controllerblock": [13, 14, 15, 16, 17], "self": 13, "throttletrim": [13, 18], "longstktrim": [13, 18], "sason": [13, 14, 15, 16, 17, 18], "fals": [13, 18, 19], "apon": [13, 14, 15, 16, 17, 18], "event_t": [13, 14, 15, 16, 17], "update_equation_funct": 13, "event_channel": 13, "__call__": 13, "event_equation_funct": 13, "throttl": [13, 14, 18, 19], "longstk": [13, 14, 18, 19], "latstk": [13, 18], "pedal": [13, 18], "command": [13, 14, 15, 16, 17, 18], "alt": [13, 18], "pb": [13, 18], "qb": [13, 18], "rb": [13, 18], "airspe": [13, 18], "keascmd": [13, 18, 19], "altcmd": [13, 18], "latoffset": [13, 17], "basechicmd": 13, "vequiv": [13, 18], "control_eart": [13, 18], "experienc": 13, "longitudin": 13, "stick": 13, "eval_trim": [13, 18, 19], "flight_condit": 13, "kin_out": 13, "controller_func": 13, "aero_plus_prop_acceler": 13, "dynamics_output_funct": 13, "gen_accel": 13, "squeez": 13, "redisu": 13, "run_trimm": [13, 18, 19], "flight_ic_arg": 13, "throttle_": [13, 18, 19], "longstk_ic": [13, 18, 19], "allow_rol": [13, 18, 19], "len_var": 13, "initial_guess": 13, "extra_index": 13, "parse_x": 13, "weighting_matrix": 13, "ey": 13, "aileron": 13, "rudder": 13, "trim_opt_func": 13, "eval_arg": 13, "copi": 13, "linalg": 13, "norm": 13, "ord": 13, "opt_r": 13, "tol": 13, "1e": 13, "12": [13, 14, 19], "disp": 13, "fatol": 13, "maxit": 13, "20_000": 13, "xatol": 13, "nelder": 13, "mead": 13, "opt_theta": 13, "opt_phi": 13, "opt_longstk": 13, "opt_throttl": 13, "opt_result": 13, "opt_arg": [13, 18, 19], "opt_flight_condit": 13, "print": 13, "4e": 13, "4f": 13, "100": [13, 14, 18, 19], "n": 13, "reshap": 13, "trimmer": 13, "determin": 13, "opt_ctrl": [13, 14, 15, 16, 17, 18, 19], "trimmed_flight_condit": [13, 18, 19], "trimmed_kea": [13, 15, 18, 19], "successfulli": [13, 14, 18, 19], "current": [13, 14, 18, 19], "006500": [13, 14], "299": [13, 14], "562": [13, 14], "6351e": [13, 14], "00": [13, 14, 18, 19, 20], "0000e": [13, 14, 18, 19], "9236": [13, 14], "7561": [13, 14], "59610112e": [13, 14], "59610020e": [13, 14], "01317177e": [13, 14], "98437663e": [13, 14], "11189587e": [13, 14], "12618714e": [13, 14], "controller_block": 13, "keascmdoutput": 13, "keascmdblock": [13, 15], "systemfromcal": [13, 14, 15, 16, 18, 19], "lambda": 13, "altcmdoutput": 13, "altcmdblock": [13, 14], "basechicmdoutput": 13, "basechicmdblock": [13, 16, 17], "latoffsetstateequ": [13, 17], "chi_cmd": 13, "v_ground_magnitud": 13, "v_ground_head": 13, "arctan2": 13, "sin": [13, 18], "latoffsetoutputequ": [13, 17], "latoffsetblock": [13, 17], "variou": 13, "appropri": [13, 14, 15, 16, 17], "assess": 13, "nstep": [13, 18], "5_000": 13, "__name__": [13, 18], "__main__": [13, 18], "63": 13, "108": 13, "739": [13, 20], "nesc_case11": [13, 14, 15, 16, 17, 18, 19, 20], "multidimension": [14, 15, 16, 17], "tabl": [14, 15, 16, 17], "look": [14, 15, 16, 17], "up": [14, 15, 16, 17], "5s": [14, 15], "increas": 14, "modifi": [14, 15, 16, 17], "interpol": [14, 15, 16, 17, 18, 19], "plot_f16_control": [14, 15, 16, 17, 18, 19], "make_interp_splin": [14, 15, 16], "k": [14, 15, 16, 18, 19], "26": [14, 15], "182": 14, "13p1": 14, "32": [14, 20], "005": [14, 20], "nesc_case13p1": [14, 20], "decreas": 15, "kea": 15, "150": 15, "13p2": 15, "136": [15, 20], "nesc_case13p2": [15, 20], "41": [16, 20], "584": 16, "13p3": 16, "46": [16, 20], "866": [16, 20], "nesc_case13p3": [16, 20], "latoffsetoutputequationshift": 17, "v_n_idx": 17, "v_e_idx": 17, "60": [17, 18], "85": 17, "191": 17, "13p4": 17, "536": [17, 20], "nesc_case13p4": [17, 20], "propag": 18, "circumnavig": [18, 19], "engag": [18, 19], "os": [18, 19], "matplotlib": [18, 19], "pyplot": [18, 19], "plt": [18, 19], "get_baselin": [18, 19], "nesc_opt": [18, 19], "nesc_color": [18, 19], "f16_gnc": [18, 19], "trimmedkea": [18, 19], "89": 18, "95": [18, 19], "10_000": [18, 19], "563": [18, 19], "643": [18, 19], "earth_output_for_gnc_select": [18, 19], "phi_d_idx": [18, 19], "lamda_d_idx": [18, 19], "get_gnc_funct": [18, 19], "circlepolesw": [18, 19], "gnc_function": 18, "100_000": 18, "gnc_block": [18, 19], "y_idx_offset": [18, 19], "xy_for_north_pole_ground_track": 18, "lat": 18, "long": 18, "xx": 18, "co": 18, "yy": 18, "subplot": [18, 19], "constrained_layout": [18, 19], "axi": [18, 19], "equal": [18, 19], "sim_lat": [18, 19], "sim_long": [18, 19], "refence_lat": 18, "reference_long": 18, "360": 18, "ref_step": 18, "ref_lag_long_list": 18, "rang": 18, "append": 18, "j": 18, "baseline_pd": [18, 19], "baseline_pd_label": [18, 19], "o": [18, 19], "markerfacecolor": [18, 19], "markeredgecolor": [18, 19], "baseline_idx": [18, 19], "enumer": [18, 19], "latitude_deg": [18, 19], "longitude_deg": [18, 19], "nesc": [18, 19], "iloc": 18, "legend": 18, "interactive_mod": [18, 19], "show": [18, 19], "savefig": [18, 19], "path": [18, 19], "join": [18, 19], "save_relative_path": [18, 19], "15_groundtrack": 18, "pdf": [18, 19], "307430": 18, "263": 18, "545": 18, "6742e": 18, "0083": 18, "8357": 18, "30742973e": 18, "36137054e": 18, "05938035e": 18, "38762619e": 18, "03070072e": 18, "56269071e": 18, "647": 18, "205": 18, "54": [18, 20], "124": [18, 20], "nesc_case15": [18, 19, 20], "sign": 19, "179": 19, "000000": 19, "279": 19, "518": 19, "6659e": 19, "9899": 19, "8207": 19, "16304744e": 19, "40508088e": 19, "68081724e": 19, "84840687e": 19, "49043324e": 19, "45023664e": 19, "362": 19, "baseline_long": 19, "unwrap": 19, "baseline_lat": 19, "xlabel": 19, "deg": 19, "ylabel": 19, "grid": 19, "16_groundtrack": 19, "370": [19, 20], "nesc_case16": [19, 20], "031": 20, "case": 20, "mb": 20}, "objects": {"": [[0, 0, 0, "-", "simupy_flight"]], "simupy_flight": [[0, 1, 1, "", "Planet"], [0, 1, 1, "", "Planetodetic"], [0, 1, 1, "", "Vehicle"], [0, 3, 1, "", "atmosphere_1976"], [0, 3, 1, "", "body_to_NED_dcm"], [0, 3, 1, "", "body_to_wind_dcm"], [0, 3, 1, "", "earth_J2_gravity"], [0, 3, 1, "", "get_constant_aero"], [0, 3, 1, "", "get_constant_atmosphere"], [0, 3, 1, "", "get_constant_force_moments"], [0, 3, 1, "", "get_constant_winds"], [0, 3, 1, "", "get_flat_pc2pd"], [0, 3, 1, "", "get_flat_pd2pc"], [0, 3, 1, "", "get_nonflat_pc2pd"], [0, 3, 1, "", "get_nonflat_pd2pc"], [0, 3, 1, "", "get_spherical_gravity"], [0, 3, 1, "", "inertial_to_body_dcm"]], "simupy_flight.Planet": [[0, 2, 1, "", "ic_from_planetodetic"], [0, 2, 1, "", "inertial_to_NED_dcm"], [0, 2, 1, "", "local_translational_trim_residual"], [0, 2, 1, "", "output_equation_function"], [0, 2, 1, "", "prepare_to_integrate"], [0, 2, 1, "", "state_equation_function"]], "simupy_flight.Vehicle": [[0, 2, 1, "", "mrc_to_com_cpm"], [0, 2, 1, "", "output_equation_function"], [0, 2, 1, "", "prepare_to_integrate"], [0, 2, 1, "", "tot_aero_forces_moments"]]}, "objtypes": {"0": "py:module", "1": "py:class", "2": "py:method", "3": "py:function"}, "objnames": {"0": ["py", "module", "Python module"], "1": ["py", "class", "Python class"], "2": ["py", "method", "Python method"], "3": ["py", "function", "Python function"]}, "titleterms": {"simupi": [0, 1], "flight": [0, 1, 13, 18, 19], "api": 0, "vehicl": 1, "toolkit": 1, "librari": 1, "design": 1, "instal": 1, "nesc": [1, 2], "test": [1, 2], "case": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "daveml": 1, "pars": 1, "contribut": 1, "licens": 1, "notic": 1, "disclaim": 1, "content": 1, "1": [3, 14], "drop": [3, 6, 7, 8, 9, 10], "sphere": [3, 6, 7, 8, 9, 10, 11, 12], "drag": [3, 4, 5], "2": [4, 15], "tumbl": [4, 5], "brick": [4, 5], "damp": [4, 5], "3": [5, 16], "dynam": 5, "4": [6, 17], "over": [6, 7, 8], "non": 6, "rotat": [6, 7, 8], "spheric": [6, 7], "earth": [6, 7, 8, 13], "5": 7, "6": 8, "ellipsoid": 8, "7": 9, "through": [9, 10], "steadi": 9, "wind": [9, 10], "field": [9, 10], "8": 10, "vari": 10, "9": 11, "launch": [11, 12], "ballist": [11, 12], "eastward": 11, "along": [11, 12], "equat": [11, 19], "10": 12, "northward": 12, "prime": 12, "meridian": 12, "11": 13, "subson": [13, 14, 15, 16, 17], "f": [13, 18, 19], "16": [13, 18, 19], "trim": 13, "across": 13, "13": [14, 15, 16, 17], "altitud": 14, "chang": [14, 15, 16], "aircraft": [14, 15, 16, 17], "veloc": 15, "cours": 16, "later": 17, "offset": 17, "manuev": 17, "15": 18, "circular": [18, 19], "around": [18, 19], "north": 18, "pole": 18, "datelin": 19, "intersect": 19, "comput": 20, "time": 20}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx": 56}})