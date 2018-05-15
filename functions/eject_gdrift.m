% [sys_vibr, g_eject] = Gdrift_eject(sys_nodelay)
%
% Removes the low frequency pole zero pair from the stage model.
% This function is very non-robust, so if the model changes it may fail.
% Check the results.
%
% Returns:
%   sys_vibr: a state space model of just the high freq. vibrational
%   dynamics.
%   g_eject: a zpk model of the INVERSE low freq pole zero pair.

function [sys_vibr, g_eject] = eject_gdrift(sys_nodelay, normalize_dc)
    Ts = sys_nodelay.Ts;
    [wp_real_x, wz_real_x] = w_zp_real(sys_nodelay);
%     rho_1 = wz_real_x(1)/wp_real_x(1);
    g_eject = zpk(exp(-wp_real_x(1)*Ts), exp(-wz_real_x(1)*Ts), 1, Ts);
    if ~exist('normalize_dc', 'var') || normalize_dc
      g_eject = g_eject/dcgain(g_eject);
    end
    sys_vibr = minreal(sys_nodelay*g_eject);
end
