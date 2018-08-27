function [ yerr ] = fit_gdrift(theta, gvib, y_exp, u_exp, t_exp, np)
%FIT_GDRIFT_CLDATA Summary of this function goes here
%   Detailed explanation goes here


gdrift = zpk(theta(np+1:end-1), theta(1:np), theta(end), gvib.Ts);

yerr = lsim(gvib*gdrift, u_exp, t_exp) - y_exp;


end

