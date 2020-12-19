from simupy import systems
import os
from scipy import interpolate
from nesc_testcase_helper import plot_nesc_comparisons, plot_F16_controls, nesc_options, ft_per_m
from nesc_case11 import int_opts, get_controller_function, BD, spec_ic_args, opt_ctrl, dim_feedback

int_opts['max_step'] = 0.0 #2**-4
int_opts['name'] = 'dop853'

glob_path = os.path.join(data_relative_path, 'Atmospheric_checkcases', 'Atmos_13p1_SubsonicAltitudeChangeF16', 'Atmos_13p1_sim_*.csv')

altCmdBlock = systems.SystemFromCallable(interpolate.make_interp_spline([0, 5], [spec_ic_args['h']*ft_per_m, spec_ic_args['h']*ft_per_m+100.0], k=0), 0, 1)
    
BD.systems[-3] = altCmdBlock
BD.systems[2] = systems.SystemFromCallable(get_controller_function(*opt_ctrl, sasOn=True, apOn=True), dim_feedback + 4, 4)
import time
tstart = time.time()
res = BD.simulate(20, integrator_options=int_opts)
tend = time.time()
tdelta = tend - tstart
print("time to simulate: %f    eval time to run time: %f" % (tdelta, res.t[-1]/tdelta))
    

plot_nesc_comparisons(res, glob_path, '13p1')
plot_F16_controls(res, '13p1')
