

function sys_sim = zero_eject(sys_nodelay)
    Ts = sys_nodelay.Ts;
    [wp_real_x, wz_real_x] = w_zp_real(sys_nodelay);
    rho_1 = wz_real_x(1)/wp_real_x(1);
    g_eject = zpk(exp(-wp_real_x(1)*Ts), exp(-wz_real_x(1)*Ts), 1, Ts);
    g_eject = g_eject/dcgain(g_eject);
    sys_eject = minreal(sys_nodelay*g_eject);
    sys_sim = SSTools.deltaUkSys(sys_eject);
end
