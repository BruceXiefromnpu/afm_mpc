% [sys_vibr, g_eject] = Gdrift_eject(sys, normalize_dc)
%
% Removes the low frequency pole zero pair from the stage model.
% This function is very non-robust, so if the model changes it may fail.
% Check the results.
%
% Inputs
% ------
%  sys : The base system from which to eject the first real pole/zero pair.
%        This should probably be a system without delay, but I havent tried to
%        leave it in. 
%  normalize_dc : (boolean, optional) If true then normalize g_drift such that
%                 it has a dc-gain of 1. Default is true (ie, if you do not 
%                 include the parameter normalize_dc, the dc-gain will be 
%                 normailized to 1).

% Returns:
%   sys_vibr: a state space model of just the high freq. vibrational
%   dynamics.
%   g_eject: a zpk model of the INVERSE low freq pole zero pair.

function [sys_vibr, g_eject] = eject_gdrift(sys, normalize_dc)
    Ts = sys.Ts;
    [wp_real_x, wz_real_x] = w_zp_real(sys);
%     rho_1 = wz_real_x(1)/wp_real_x(1);
    g_eject = zpk(exp(-wp_real_x(1)*Ts), exp(-wz_real_x(1)*Ts), 1, Ts);
    if ~exist('normalize_dc', 'var') || normalize_dc
      g_eject = g_eject/dcgain(g_eject);
    end
    sys_vibr = minreal(sys*g_eject);
end
