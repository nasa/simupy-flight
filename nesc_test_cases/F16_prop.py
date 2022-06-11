import numpy
import ndsplines


MIL_PWR = 5.000000e+01
FEY = 0.000000e+00
FEZ = 0.000000e+00
TEL = 0.000000e+00
TEM = 0.000000e+00
TEN = 0.000000e+00


T_IDLE_fn_defn_knots_0 = numpy.array([0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0])
T_IDLE_fn_defn_knots_1 = numpy.array([0.0, 0.0, 10000.0, 20000.0, 30000.0, 40000.0, 50000.0, 50000.0])
T_IDLE_fn_defn_coeffs = numpy.array([[1060.0, 670.0, 880.0, 1140.0, 1500.0, 1860.0],
 [635.0, 425.0, 690.0, 1010.0, 1330.0, 1700.0],
 [60.0, 25.0, 345.0, 755.0, 1130.0, 1525.0],
 [-1020.0, -710.0, -300.0, 350.0, 910.0, 1360.0],
 [-2700.0, -1900.0, -1300.0, -247.0, 600.0, 1100.0],
 [-3600.0, -1400.0, -595.0, -342.0, -200.0, 700.0]])
T_IDLE_fn_defn = ndsplines.NDSpline([T_IDLE_fn_defn_knots_0,T_IDLE_fn_defn_knots_1], T_IDLE_fn_defn_coeffs, 1)

T_MIL_fn_defn_knots_0 = numpy.array([0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0])
T_MIL_fn_defn_knots_1 = numpy.array([0.0, 0.0, 10000.0, 20000.0, 30000.0, 40000.0, 50000.0, 50000.0])
T_MIL_fn_defn_coeffs = numpy.array([[12680.0, 9150.0, 6200.0, 3950.0, 2450.0, 1400.0],
 [12680.0, 9150.0, 6313.0, 4040.0, 2470.0, 1400.0],
 [12610.0, 9312.0, 6610.0, 4290.0, 2600.0, 1560.0],
 [12640.0, 9839.0, 7090.0, 4660.0, 2840.0, 1660.0],
 [12390.0, 10176.0, 7750.0, 5320.0, 3250.0, 1930.0],
 [11680.0, 9848.0, 8050.0, 6100.0, 3800.0, 2310.0]])
T_MIL_fn_defn = ndsplines.NDSpline([T_MIL_fn_defn_knots_0,T_MIL_fn_defn_knots_1], T_MIL_fn_defn_coeffs, 1)

T_MAX_fn_defn_knots_0 = numpy.array([0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0])
T_MAX_fn_defn_knots_1 = numpy.array([0.0, 0.0, 10000.0, 20000.0, 30000.0, 40000.0, 50000.0, 50000.0])
T_MAX_fn_defn_coeffs = numpy.array([[20000.0, 15000.0, 10800.0, 7000.0, 4000.0, 2500.0],
 [21420.0, 15700.0, 11225.0, 7323.0, 4435.0, 2600.0],
 [22700.0, 16860.0, 12250.0, 8154.0, 5000.0, 2835.0],
 [24240.0, 18910.0, 13760.0, 9285.0, 5700.0, 3215.3],
 [28070.0, 21075.0, 15975.0, 11115.0, 6860.0, 3950.0],
 [28885.0, 23319.0, 18300.0, 13484.0, 8642.0, 5057.0]])
T_MAX_fn_defn = ndsplines.NDSpline([T_MAX_fn_defn_knots_0,T_MAX_fn_defn_knots_1], T_MAX_fn_defn_coeffs, 1)



def F16_prop(PWR, ALT, RMACH):
    [T_IDLE] = T_IDLE_fn_defn((RMACH, ALT)).ravel()
    [T_MIL] = T_MIL_fn_defn((RMACH, ALT)).ravel()
    [T_MAX] = T_MAX_fn_defn((RMACH, ALT)).ravel()
    FEX = numpy.select([numpy.greater_equal(MIL_PWR, PWR)], [T_IDLE + PWR*(-T_IDLE + T_MIL)/MIL_PWR], default=T_MIL + (-MIL_PWR + PWR)*(T_MAX - T_MIL)/(100.0 - MIL_PWR))
    return (numpy.array([FEX, FEY, FEZ, TEL, TEM, TEN]))

check_data = [
 [[0.0, 0.0, 0.0],
 [1060.0, 0.0, 0.0, 0.0, 0.0, 0.0]],

 [[50.0, 0.0, 0.0],
 [12680.0, 0.0, 0.0, 0.0, 0.0, 0.0]],

 [[100.0, 0.0, 0.0],
 [20000.0, 0.0, 0.0, 0.0, 0.0, 0.0]],

 [[100.0, 0.0, 1.0],
 [28885.0, 0.0, 0.0, 0.0, 0.0, 0.0]],

 [[0.0, 50000.0, 1.0],
 [700.0, 0.0, 0.0, 0.0, 0.0, 0.0]],

 [[50.0, 50000.0, 1.0],
 [2310.0, 0.0, 0.0, 0.0, 0.0, 0.0]],

 [[100.0, 50000.0, 1.0],
 [5057.0, 0.0, 0.0, 0.0, 0.0, 0.0]],

 [[42.3, 23507.0, 0.625],
 [5319.3491, 0.0, 0.0, 0.0, 0.0, 0.0]],

 [[88.3, 33537.0, 0.895],
 [9298.8926, 0.0, 0.0, 0.0, 0.0, 0.0]] 
]


def run_checks():
    for check_input, check_output in check_data:
        if not numpy.allclose( F16_prop(*check_input), check_output):
            raise ValueError("Check for F16_prop failed!")
    print("All checks for F16_prop passed.")

if __name__ == "__main__":
    run_checks()
