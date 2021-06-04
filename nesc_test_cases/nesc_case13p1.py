from simupy import systems
from scipy import interpolate

from nesc_testcase_helper import plot_nesc_comparisons, plot_F16_controls, benchmark
from nesc_testcase_helper import ft_per_m
from nesc_case11 import int_opts, get_controller_function, BD, spec_ic_args, opt_ctrl, dim_feedback

int_opts['max_step'] = 0.0 #2**-4
int_opts['name'] = 'dop853'

altCmdBlock = systems.SystemFromCallable(interpolate.make_interp_spline([0, 5], [spec_ic_args['h']*ft_per_m, spec_ic_args['h']*ft_per_m+100.0], k=0), 0, 1)

BD.systems[-3] = altCmdBlock
BD.systems[2] = systems.SystemFromCallable(get_controller_function(*opt_ctrl, sasOn=True, apOn=True), dim_feedback + 4, 4)

with benchmark() as b:
    res = BD.simulate(20, integrator_options=int_opts)
    b.tfinal = res.t[-1]

plot_nesc_comparisons(res, '13p1')
plot_F16_controls(res, '13p1')
