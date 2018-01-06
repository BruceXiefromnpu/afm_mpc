% Get the frequency response of an observer in closed loop.
% 
%  H = -K*[zI - A + BK + LC]^-1 * L
%
% sys_frf is optional.

function [ kalman_frf, w_s ] = observer_frf(sys, K, L, sys_frf, w_s)

if ~exist('sys_frf', 'var')
    [sys_frf, w_s] = getFRF(sys);
end


A_tilde = sys.a - sys.b*K - L*sys.c;

sys_kalman = ss(A_tilde, L, K, 0, sys.Ts);

kalman_frf = getFRF(sys_kalman, w_s);




end

