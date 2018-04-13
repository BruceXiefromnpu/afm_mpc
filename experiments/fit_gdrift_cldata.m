function [ yerr ] = fit_gdrift_cldata(theta, gvib, y_exp, u_exp, t_exp)
%FIT_GDRIFT_CLDATA Summary of this function goes here
%   Detailed explanation goes here


gdrift = zpk(theta(1), theta(2), theta(3), gvib.Ts);

yerr = lsim(gvib*gdrift, u_exp, t_exp) - y_exp;


end

