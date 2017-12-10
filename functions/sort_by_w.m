% Sorts a set of complex numbers by the equivalent s-plane natural
% frequency. This is useful if the complex numbers are the poles or zeros
% of a discrete time system.
%
% I sort them such that wn_i = | log(z_i) |is decreasing.
% 

function [ p ] = sort_by_w(p)

logP = abs(log(p));

[pp, inds] = sort(logP);

p = p(inds);

end

