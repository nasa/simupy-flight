def plot_nesc_comparisons(res, baseline_pds):
    long_lat_deg = res.y[:,13:15]*180/np.pi
    alt_ft = res.y[:,15]*ft_per_m
    
    psi_theta_phi_deg = res.y[:,16:19]*180/np.pi
    fig, ax = plt.subplots(3, sharex=True, constrained_layout=True)
    for idx, baseline_pd in enumerate(baseline_pds):
        try:
            baseline_pd[['longitude_deg', 'latitude_deg', 'altitudeMsl_ft']].plot(subplots=True, color='C%d'%idx, alpha=0.5, legend=False, ax=ax, label='')
        except KeyError:
            print('missing sim %d' % idx)
            pass
        
    ax[0].plot(res.t, long_lat_deg[:, 0], 'k--', label='SimuPy-based')
    ax[0].set_ylabel('longitude, deg')
    ax[0].legend(['NESC %d' % (idx+1) for idx in range(len(baseline_pds))] + ['SimuPy-based'])
    ax[1].plot(res.t, long_lat_deg[:, 1], 'k--')
    ax[1].set_ylabel('latitude, deg')
    ax[2].plot(res.t, alt_ft, 'k--')
    ax[2].set_ylabel('altitude, ft')
    ax[2].set_xlabel('time, s')
    
    for axis in ax:
        axis.grid(True)
    plt.show()
    
    fig, ax = plt.subplots(3, sharex=True, constrained_layout=True)
    for idx, baseline_pd in enumerate(baseline_pds):
        try:
            baseline_pd[['eulerAngle_deg_Yaw', 'eulerAngle_deg_Pitch', 'eulerAngle_deg_Roll']].plot(subplots=True, color='C%d'%idx, alpha=0.5, legend=False, ax=ax)
        except KeyError:
            print('missing sim %d' % idx)
            pass
        
    ax[0].plot(res.t, psi_theta_phi_deg[:, 0], 'k--')
    ax[1].plot(res.t, psi_theta_phi_deg[:, 1], 'k--')
    ax[2].plot(res.t, psi_theta_phi_deg[:, 2], 'k--')
    
    ax[0].set_ylabel('yaw, deg')
    ax[1].set_ylabel('pitch, deg')
    ax[2].set_ylabel('roll, deg')
    
    ax[0].set_ylabel('$\\psi$, deg')
    ax[1].set_ylabel('$\\theta$, deg')
    ax[2].set_ylabel('$\\phi$, deg')
    
    ax[2].set_xlabel('time, s')
    for axis in ax:
        axis.grid(True)
    plt.show()
    
    fig, ax = plt.subplots(3, sharex=True, constrained_layout=True)
    for idx, baseline_pd in enumerate(baseline_pds):
        try:
            baseline_pd[['bodyAngularRateWrtEi_deg_s_Roll', 'bodyAngularRateWrtEi_deg_s_Pitch', 'bodyAngularRateWrtEi_deg_s_Yaw']].plot(subplots=True, color='C%d'%idx, alpha=0.5, legend=False, ax=ax)
        except KeyError:
            print('missing sim %d' % idx)
            pass
    
    ax[0].plot(res.t, res.x[:, 10]*180/np.pi, 'k--')
    ax[0].set_ylabel('$p$, deg/s')
    ax[1].plot(res.t, res.x[:, 11]*180/np.pi, 'k--')
    ax[1].set_ylabel('$q$, deg/s')
    ax[2].plot(res.t, res.x[:, 12]*180/np.pi, 'k--')
    ax[2].set_ylabel('$r$, deg/s')
    ax[2].set_xlabel('time, s')
    for axis in ax:
        axis.grid(True)
    plt.show()