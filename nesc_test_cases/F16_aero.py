import numpy
import ndsplines


cbar = 1.132000e+01
bspan = 3.000000e+01
sref = 3.000000e+02


CX_fn_knots_0 = numpy.array([-24.0, -24.0, -12.0, 0.0, 12.0, 24.0, 24.0])
CX_fn_knots_1 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
CX_fn_coeffs = numpy.array([[-0.099, -0.081, -0.081, -0.063, -0.025, 0.044, 0.097, 0.113, 0.145, 0.167, 0.174, 0.166],
 [-0.048, -0.038, -0.04, -0.021, 0.016, 0.083, 0.127, 0.137, 0.162, 0.177, 0.179, 0.167],
 [-0.022, -0.02, -0.021, -0.004, 0.032, 0.094, 0.128, 0.13, 0.154, 0.161, 0.155, 0.138],
 [-0.04, -0.038, -0.039, -0.025, 0.006, 0.062, 0.087, 0.085, 0.1, 0.11, 0.104, 0.091],
 [-0.083, -0.073, -0.076, -0.072, -0.046, 0.012, 0.024, 0.025, 0.043, 0.053, 0.047, 0.04]])
CX_fn = ndsplines.NDSpline([CX_fn_knots_0,CX_fn_knots_1], CX_fn_coeffs, 1)

CZ0_fn_knots_0 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
CZ0_fn_coeffs = numpy.array([0.77,
 0.241,
 -0.1,
 -0.416,
 -0.731,
 -1.053,
 -1.366,
 -1.646,
 -1.917,
 -2.12,
 -2.248,
 -2.229])
CZ0_fn = ndsplines.NDSpline([CZ0_fn_knots_0], CZ0_fn_coeffs, 1)

Cm0_fn_knots_0 = numpy.array([-24.0, -24.0, -12.0, 0.0, 12.0, 24.0, 24.0])
Cm0_fn_knots_1 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
Cm0_fn_coeffs = numpy.array([[0.205, 0.168, 0.186, 0.196, 0.213, 0.251, 0.245, 0.238, 0.252, 0.231, 0.198, 0.192],
 [0.081, 0.077, 0.107, 0.11, 0.11, 0.141, 0.127, 0.119, 0.133, 0.108, 0.081, 0.093],
 [-0.046, -0.02, -0.009, -0.005, -0.006, 0.01, 0.006, -0.001, 0.014, 0.0, -0.013, 0.032],
 [-0.174, -0.145, -0.121, -0.127, -0.129, -0.102, -0.097, -0.113, -0.087, -0.084, -0.069, -0.006],
 [-0.259, -0.202, -0.184, -0.193, -0.199, -0.15, -0.16, -0.167, -0.104, -0.076, -0.041, -0.005]])
Cm0_fn = ndsplines.NDSpline([Cm0_fn_knots_0,Cm0_fn_knots_1], Cm0_fn_coeffs, 1)

Cl0_fn_knots_0 = numpy.array([0.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 30.0])
Cl0_fn_knots_1 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
Cl0_fn_coeffs = numpy.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
 [-0.001, -0.004, -0.008, -0.012, -0.016, -0.022, -0.022, -0.021, -0.015, -0.008, -0.013, -0.015],
 [-0.003, -0.009, -0.017, -0.024, -0.03, -0.041, -0.045, -0.04, -0.016, -0.002, -0.01, -0.019],
 [-0.001, -0.01, -0.02, -0.03, -0.039, -0.054, -0.057, -0.054, -0.023, -0.006, -0.014, -0.027],
 [0.0, -0.01, -0.022, -0.034, -0.047, -0.06, -0.069, -0.067, -0.033, -0.036, -0.035, -0.035],
 [0.007, -0.01, -0.023, -0.034, -0.049, -0.063, -0.081, -0.079, -0.06, -0.058, -0.062, -0.059],
 [0.009, -0.011, -0.023, -0.037, -0.05, -0.068, -0.089, -0.088, -0.091, -0.076, -0.077, -0.076]])
Cl0_fn = ndsplines.NDSpline([Cl0_fn_knots_0,Cl0_fn_knots_1], Cl0_fn_coeffs, 1)

Cn0_fn_knots_0 = numpy.array([0.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 30.0])
Cn0_fn_knots_1 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
Cn0_fn_coeffs = numpy.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
 [0.018, 0.019, 0.018, 0.019, 0.019, 0.018, 0.013, 0.007, 0.004, -0.014, -0.017, -0.033],
 [0.038, 0.042, 0.042, 0.042, 0.043, 0.039, 0.03, 0.017, 0.004, -0.035, -0.047, -0.057],
 [0.056, 0.057, 0.059, 0.058, 0.058, 0.053, 0.032, 0.012, 0.002, -0.046, -0.071, -0.073],
 [0.064, 0.077, 0.076, 0.074, 0.073, 0.057, 0.029, 0.007, 0.012, -0.034, -0.065, -0.041],
 [0.074, 0.086, 0.093, 0.089, 0.08, 0.062, 0.049, 0.022, 0.028, -0.012, -0.002, -0.013],
 [0.079, 0.09, 0.106, 0.106, 0.096, 0.08, 0.068, 0.03, 0.064, 0.015, 0.011, -0.001]])
Cn0_fn = ndsplines.NDSpline([Cn0_fn_knots_0,Cn0_fn_knots_1], Cn0_fn_coeffs, 1)

CXq_fn_knots_0 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
CXq_fn_coeffs = numpy.array([-0.267,
 -0.11,
 0.308,
 1.34,
 2.08,
 2.91,
 2.76,
 2.05,
 1.5,
 1.49,
 1.83,
 1.21])
CXq_fn = ndsplines.NDSpline([CXq_fn_knots_0], CXq_fn_coeffs, 1)

CYr_fn_knots_0 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
CYr_fn_coeffs = numpy.array([0.882,
 0.852,
 0.876,
 0.958,
 0.962,
 0.974,
 0.819,
 0.483,
 0.59,
 1.21,
 -0.493,
 -1.04])
CYr_fn = ndsplines.NDSpline([CYr_fn_knots_0], CYr_fn_coeffs, 1)

CYp_fn_knots_0 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
CYp_fn_coeffs = numpy.array([-0.108,
 -0.108,
 -0.188,
 0.11,
 0.258,
 0.226,
 0.344,
 0.362,
 0.611,
 0.529,
 0.298,
 -0.227])
CYp_fn = ndsplines.NDSpline([CYp_fn_knots_0], CYp_fn_coeffs, 1)

CZq_fn_knots_0 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
CZq_fn_coeffs = numpy.array([-8.8,
 -25.8,
 -28.9,
 -31.4,
 -31.2,
 -30.7,
 -27.7,
 -28.2,
 -29.0,
 -29.8,
 -38.3,
 -35.3])
CZq_fn = ndsplines.NDSpline([CZq_fn_knots_0], CZq_fn_coeffs, 1)

Clr_fn_knots_0 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
Clr_fn_coeffs = numpy.array([-0.126,
 -0.026,
 0.063,
 0.113,
 0.208,
 0.23,
 0.319,
 0.437,
 0.68,
 0.1,
 0.447,
 -0.33])
Clr_fn = ndsplines.NDSpline([Clr_fn_knots_0], Clr_fn_coeffs, 1)

Clp_fn_knots_0 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
Clp_fn_coeffs = numpy.array([-0.36,
 -0.359,
 -0.443,
 -0.42,
 -0.383,
 -0.375,
 -0.329,
 -0.294,
 -0.23,
 -0.21,
 -0.12,
 -0.1])
Clp_fn = ndsplines.NDSpline([Clp_fn_knots_0], Clp_fn_coeffs, 1)

Cmq_fn_knots_0 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
Cmq_fn_coeffs = numpy.array([-7.21,
 -5.4,
 -5.23,
 -5.26,
 -6.11,
 -6.64,
 -5.69,
 -6.0,
 -6.2,
 -6.4,
 -6.6,
 -6.0])
Cmq_fn = ndsplines.NDSpline([Cmq_fn_knots_0], Cmq_fn_coeffs, 1)

Cnr_fn_knots_0 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
Cnr_fn_coeffs = numpy.array([-0.38,
 -0.363,
 -0.378,
 -0.386,
 -0.37,
 -0.453,
 -0.55,
 -0.582,
 -0.595,
 -0.637,
 -1.02,
 -0.84])
Cnr_fn = ndsplines.NDSpline([Cnr_fn_knots_0], Cnr_fn_coeffs, 1)

Cnp_fn_knots_0 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
Cnp_fn_coeffs = numpy.array([0.061,
 0.052,
 0.052,
 -0.012,
 -0.013,
 -0.024,
 0.05,
 0.15,
 0.13,
 0.158,
 0.24,
 0.15])
Cnp_fn = ndsplines.NDSpline([Cnp_fn_knots_0], Cnp_fn_coeffs, 1)

dlda_fn_knots_0 = numpy.array([-30.0, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 30.0])
dlda_fn_knots_1 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
dlda_fn_coeffs = numpy.array([[-0.041, -0.052, -0.053, -0.056, -0.05, -0.056, -0.082, -0.059, -0.042, -0.038, -0.027, -0.017],
 [-0.041, -0.053, -0.053, -0.053, -0.05, -0.051, -0.066, -0.043, -0.038, -0.027, -0.023, -0.016],
 [-0.042, -0.053, -0.052, -0.051, -0.049, -0.049, -0.043, -0.035, -0.026, -0.016, -0.018, -0.014],
 [-0.04, -0.052, -0.051, -0.052, -0.048, -0.048, -0.042, -0.037, -0.031, -0.026, -0.017, -0.012],
 [-0.043, -0.049, -0.048, -0.049, -0.043, -0.042, -0.042, -0.036, -0.025, -0.021, -0.016, -0.011],
 [-0.044, -0.048, -0.048, -0.047, -0.042, -0.041, -0.02, -0.028, -0.013, -0.014, -0.011, -0.01],
 [-0.043, -0.049, -0.047, -0.045, -0.042, -0.037, -0.003, -0.013, -0.01, -0.003, -0.007, -0.008]])
dlda_fn = ndsplines.NDSpline([dlda_fn_knots_0,dlda_fn_knots_1], dlda_fn_coeffs, 1)

dldr_fn_knots_0 = numpy.array([-30.0, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 30.0])
dldr_fn_knots_1 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
dldr_fn_coeffs = numpy.array([[0.005, 0.017, 0.014, 0.01, -0.005, 0.009, 0.019, 0.005, 0.0, -0.005, -0.011, 0.008],
 [0.007, 0.016, 0.014, 0.014, 0.013, 0.009, 0.012, 0.005, 0.0, 0.004, 0.009, 0.007],
 [0.013, 0.013, 0.011, 0.012, 0.011, 0.009, 0.008, 0.005, 0.0, 0.005, 0.003, 0.005],
 [0.018, 0.015, 0.015, 0.014, 0.014, 0.014, 0.014, 0.015, 0.013, 0.011, 0.006, 0.001],
 [0.015, 0.014, 0.013, 0.013, 0.012, 0.011, 0.011, 0.01, 0.008, 0.008, 0.007, 0.003],
 [0.021, 0.011, 0.01, 0.011, 0.01, 0.009, 0.008, 0.01, 0.006, 0.005, 0.0, 0.001],
 [0.023, 0.01, 0.011, 0.011, 0.011, 0.01, 0.008, 0.01, 0.006, 0.014, 0.02, 0.0]])
dldr_fn = ndsplines.NDSpline([dldr_fn_knots_0,dldr_fn_knots_1], dldr_fn_coeffs, 1)

dnda_fn_knots_0 = numpy.array([-30.0, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 30.0])
dnda_fn_knots_1 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
dnda_fn_coeffs = numpy.array([[0.001, -0.027, -0.017, -0.013, -0.012, -0.016, 0.001, 0.017, 0.011, 0.017, 0.008, 0.016],
 [0.002, -0.014, -0.016, -0.016, -0.014, -0.019, -0.021, 0.002, 0.012, 0.016, 0.015, 0.011],
 [-0.006, -0.008, -0.006, -0.006, -0.005, -0.008, -0.005, 0.007, 0.004, 0.007, 0.006, 0.006],
 [-0.011, -0.011, -0.01, -0.009, -0.008, -0.006, 0.0, 0.004, 0.007, 0.01, 0.004, 0.01],
 [-0.015, -0.015, -0.014, -0.012, -0.011, -0.008, -0.002, 0.002, 0.006, 0.012, 0.011, 0.011],
 [-0.024, -0.01, -0.004, -0.002, -0.001, 0.003, 0.014, 0.006, -0.001, 0.004, 0.004, 0.006],
 [-0.022, 0.002, -0.003, -0.005, -0.003, -0.001, -0.009, -0.009, -0.001, 0.003, -0.002, 0.001]])
dnda_fn = ndsplines.NDSpline([dnda_fn_knots_0,dnda_fn_knots_1], dnda_fn_coeffs, 1)

dndr_fn_knots_0 = numpy.array([-30.0, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 30.0])
dndr_fn_knots_1 = numpy.array([-10.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 45.0])
dndr_fn_coeffs = numpy.array([[-0.018, -0.052, -0.052, -0.052, -0.054, -0.049, -0.059, -0.051, -0.03, -0.037, -0.026, -0.013],
 [-0.028, -0.051, -0.043, -0.046, -0.045, -0.049, -0.057, -0.052, -0.03, -0.033, -0.03, -0.008],
 [-0.037, -0.041, -0.038, -0.04, -0.04, -0.038, -0.037, -0.03, -0.027, -0.024, -0.019, -0.013],
 [-0.048, -0.045, -0.045, -0.045, -0.044, -0.045, -0.047, -0.048, -0.049, -0.045, -0.033, -0.016],
 [-0.043, -0.044, -0.041, -0.041, -0.04, -0.038, -0.034, -0.035, -0.035, -0.029, -0.022, -0.009],
 [-0.052, -0.034, -0.036, -0.036, -0.035, -0.028, -0.024, -0.023, -0.02, -0.016, -0.01, -0.014],
 [-0.062, -0.034, -0.027, -0.028, -0.027, -0.027, -0.023, -0.023, -0.019, -0.009, -0.025, -0.01]])
dndr_fn = ndsplines.NDSpline([dndr_fn_knots_0,dndr_fn_knots_1], dndr_fn_coeffs, 1)



def F16_aero(vt, alpha, beta, p, q, r, el, ail, rdr):
    rtd = 57.2957795785523
    del_ = 0.04*el
    dail = 0.05*ail
    drdr = 0.0333333333333333*rdr
    [cxt] = CX_fn((el, alpha)).ravel()
    cy0 = -0.02*beta + 0.021*dail + 0.086*drdr
    [czt] = CZ0_fn((alpha,)).ravel()
    cz1 = czt*(1.0 - (beta/rtd)**2.0) - 0.19*del_
    tvt = 2.0*vt
    b2v = bspan/tvt
    cq2v = cbar*q/tvt
    absbeta = abs(beta)
    [absCl0] = Cl0_fn((absbeta, alpha)).ravel()
    [absCn0] = Cn0_fn((absbeta, alpha)).ravel()
    clt = numpy.select([numpy.less_equal(beta, 0.0)], [-absCl0], default=absCl0)
    [cmt] = Cm0_fn((el, alpha)).ravel()
    cnt = numpy.select([numpy.less_equal(beta, 0.0)], [-absCn0], default=absCn0)
    [cxq] = CXq_fn((alpha,)).ravel()
    [cyr] = CYr_fn((alpha,)).ravel()
    [cyp] = CYp_fn((alpha,)).ravel()
    [czq] = CZq_fn((alpha,)).ravel()
    [clr] = Clr_fn((alpha,)).ravel()
    [clp] = Clp_fn((alpha,)).ravel()
    [cmq] = Cmq_fn((alpha,)).ravel()
    [cnr] = Cnr_fn((alpha,)).ravel()
    [cnp] = Cnp_fn((alpha,)).ravel()
    [dclda] = dlda_fn((beta, alpha)).ravel()
    [dcldr] = dldr_fn((beta, alpha)).ravel()
    [dcnda] = dnda_fn((beta, alpha)).ravel()
    [dcndr] = dndr_fn((beta, alpha)).ravel()
    cl1 = clt + dail*dclda + dcldr*drdr
    cn1 = cnt + dail*dcnda + dcndr*drdr
    cx = cq2v*cxq + cxt
    cy = b2v*(cyp*p + cyr*r) + cy0
    cz = cq2v*czq + cz1
    cl = b2v*(clp*p + clr*r) + cl1
    cm = cmq*cq2v + cmt
    cn = b2v*(cnp*p + cnr*r) + cn1
    return (numpy.array([cbar, bspan, sref, cx, cy, cz, cl, cm, cn]))

check_data = [
 [[300.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
 [11.32, 30.0, 300.0, -0.004, 0.0, -0.416, 0.0, -0.005, 0.0]],

 [[300.0, 5.0, 2.34, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
 [11.32, 30.0, 300.0, -0.004, -0.0468, -0.41530612733219, -0.005616, -0.005, 0.008892]],

 [[300.0, 5.0, -2.34, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
 [11.32, 30.0, 300.0, -0.004, 0.0468, -0.41530612733219, 0.005616, -0.005, -0.008892]],

 [[300.0, 5.0, 0.0, 3.42, 0.0, 0.0, 0.0, 0.0, 0.0],
 [11.32, 30.0, 300.0, -0.004, 0.01881, -0.416, -0.07182, -0.005, -0.002052]],

 [[300.0, 5.0, 0.0, -3.42, 0.0, 0.0, 0.0, 0.0, 0.0],
 [11.32, 30.0, 300.0, -0.004, -0.01881, -0.416, 0.07182, -0.005, 0.002052]],

 [[300.0, 5.0, 0.0, 0.0, 0.98, 0.0, 0.0, 0.0, 0.0],
 [11.32, 30.0, 300.0, 0.02077570666667, 0.0, -0.99656506666667, 0.0, -0.10225389333333, 0.0]],

 [[300.0, 5.0, 0.0, 0.0, -0.98, 0.0, 0.0, 0.0, 0.0],
 [11.32, 30.0, 300.0, -0.02877570666667, 0.0, 0.16456506666667, 0.0, 0.09225389333333, 0.0]],

 [[300.0, 5.0, 0.0, 0.0, 0.0, 2.92, 0.0, 0.0, 0.0],
 [11.32, 30.0, 300.0, -0.004, 0.139868, -0.416, 0.016498, -0.005, -0.056356]],

 [[300.0, 5.0, 0.0, 0.0, 0.0, -2.92, 0.0, 0.0, 0.0],
 [11.32, 30.0, 300.0, -0.004, -0.139868, -0.416, -0.016498, -0.005, 0.056356]],

 [[300.0, 5.0, 0.0, 0.0, 0.0, 0.0, 12.92, 0.0, 0.0],
 [11.32, 30.0, 300.0, -0.02860333333333, 0.0, -0.514192, 0.0, -0.13206, 0.0]],

 [[300.0, 5.0, 0.0, 0.0, 0.0, 0.0, -12.92, 0.0, 0.0],
 [11.32, 30.0, 300.0, -0.02422, 0.0, -0.317808, 0.0, 0.11659333333333, 0.0]],

 [[300.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 24.1, 0.0],
 [11.32, 30.0, 300.0, -0.004, 0.025305, -0.416, -0.06266, -0.005, -0.010845]],

 [[300.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, -24.1, 0.0],
 [11.32, 30.0, 300.0, -0.004, -0.025305, -0.416, 0.06266, -0.005, 0.010845]],

 [[300.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.03],
 [11.32, 30.0, 300.0, -0.004, 0.034486, -0.416, 0.005614, -0.005, -0.018045]],

 [[300.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.03],
 [11.32, 30.0, 300.0, -0.004, -0.034486, -0.416, -0.005614, -0.005, 0.018045]],

 [[300.0, 16.2, -3.24, 0.56, -0.76, -0.94, 4.567, 7.654, -2.991],
 [11.32, 30.0, 300.0, 0.04794994533333, 0.02735386, -0.72934852554344, -0.026917840128, 0.05917625733333, 0.013526640528]] 
]


def run_checks():
    for check_input, check_output in check_data:
        if not numpy.allclose( F16_aero(*check_input), check_output):
            raise ValueError("Check for F16_aero failed!")
    print("All checks for F16_aero passed.")