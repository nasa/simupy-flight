import numpy


thetaCmdAltErrGain = -5.000000e-02
chiCmdYErrGain = -1.000000e-02
phiCmdChiErrGain = -1.000000e+01
trimmedAlpha = 2.653814e+00
trimmedTheta = 2.653814e+00
trimmedKEAS = 2.878089e+02
longLQR11 = -6.300907e-02
longLQR12 = 1.132304e-01
longLQR13 = 1.011343e+01
longLQR14 = 3.154983e+00
longLQR21 = 9.972606e-01
longLQR22 = -2.546771e-02
longLQR23 = 1.213308e+00
longLQR24 = 2.087444e-01
latdLQR11 = 3.078044e+00
latdLQR12 = 3.236586e-02
latdLQR13 = 4.557859e+00
latdLQR14 = 5.894432e-01
latdLQR21 = -7.058175e-01
latdLQR22 = -2.563629e-01
latdLQR23 = -1.073666e+00
latdLQR24 = 8.221146e-01


def F16_control(throttle, longStk, latStk, pedal, sasOn, apOn, keasCmd, altCmd, latOffset, baseChiCmd, altMsl, Vequiv, alpha, beta, phi, theta, psi, pb, qb, rb, throttleTrim, longStkTrim):
    fsasOn = apOn + sasOn
    altErr = -altCmd + altMsl
    deltaThetaCmd = numpy.clip(altErr*thetaCmdAltErrGain,-5.0,5.0)
    thetaCmd = deltaThetaCmd + trimmedTheta
    chiEst = beta + psi
    deltaChiCmd = numpy.clip(chiCmdYErrGain*latOffset,-30.0,30.0)
    chiCmd = baseChiCmd + deltaChiCmd
    chiErrUnwrapped = -chiCmd + chiEst
    chiErr = numpy.select([numpy.greater(abs(chiErrUnwrapped), 180.0)], [numpy.select([numpy.greater(chiErrUnwrapped, 0.0)], [chiErrUnwrapped - 360.0], default=chiErrUnwrapped + 360.0)], default=chiErrUnwrapped)
    phiCmd = numpy.clip(chiErr*phiCmdChiErrGain,-30.0,30.0)
    keasCmdSw = numpy.select([numpy.greater(apOn, 0.5)], [keasCmd], default=trimmedKEAS)
    thetaCmdSw = numpy.select([numpy.greater(apOn, 0.5)], [thetaCmd], default=trimmedTheta)
    phiCmdSw = numpy.select([numpy.greater(apOn, 0.5)], [phiCmd], default=0.0)
    deltaVequiv = Vequiv - keasCmdSw
    deltaAlpha = alpha - trimmedAlpha
    deltaPhi = phi - phiCmdSw
    deltaTheta = theta - thetaCmdSw
    longLQR = -deltaAlpha*longLQR12 - deltaTheta*longLQR14 - deltaVequiv*longLQR11 - longLQR13*qb
    throttleLQR = -deltaAlpha*longLQR22 - deltaTheta*longLQR24 - deltaVequiv*longLQR21 - longLQR23*qb
    latLQR = -beta*latdLQR12 - deltaPhi*latdLQR11 - latdLQR13*pb - latdLQR14*rb
    dirLQR = -beta*latdLQR22 - deltaPhi*latdLQR21 - latdLQR23*pb - latdLQR24*rb
    longLQRsw = numpy.select([numpy.greater(fsasOn, 0.5)], [longLQR], default=0.0)
    throttleLQRsw = numpy.select([numpy.greater(fsasOn, 0.5)], [throttleLQR], default=0.0)
    latLQRsw = numpy.select([numpy.greater(fsasOn, 0.5)], [latLQR], default=0.0)
    dirLQRsw = numpy.select([numpy.greater(fsasOn, 0.5)], [dirLQR], default=0.0)
    longStkSw = numpy.select([numpy.less_equal(apOn, 0.5)], [longStk], default=0.0)
    throttleSw = numpy.select([numpy.less_equal(apOn, 0.5)], [throttle], default=0.0)
    latStkSw = numpy.select([numpy.less_equal(apOn, 0.5)], [latStk], default=0.0)
    pedalSw = numpy.select([numpy.less_equal(apOn, 0.5)], [pedal], default=0.0)
    totLongStk = numpy.clip(longLQRsw + longStkSw + longStkTrim,-1.0,1.0)
    totLatStk = numpy.clip(latLQRsw + latStkSw,-1.0,1.0)
    totPedal = numpy.clip(dirLQRsw + pedalSw,-1.0,1.0)
    totThrottle = numpy.clip(throttleLQRsw + throttleSw + throttleTrim,0.0,1.0)
    el = -25.0*totLongStk
    ail = -21.5*totLatStk
    rdr = 0.008*ail - 30.0*totPedal
    PWR = 100.0*totThrottle
    return (numpy.array([el, ail, rdr, PWR]))
