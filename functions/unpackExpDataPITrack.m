% [y_exp, u_exp, xhat_exp xdelay_exp] = unpackExpData(data, Ns, )
% 
% Unpacks the experimental data that LabView spits out into pieces that
% make sense.
%
% Inputs: 
%   data: from labview. We expect the following format:
%       dataBoth = [y,  X1.....Xns,  u,  Xd]
%   Ns:   Number of states (excluding delay)
%
% Outputs:
%   y_exp, u_exp. Actual measured input/output to AFM.
%   xhat_exp: the estimated states from the observer.
%   xdelay_exp: Input delay states, from shift register. 
%
%   These are all structs, e.g. y_exp.Data, y_exp.Time.


function [y_exp, u_exp] = unpackExpDataPITrack(data,  Ts)
    % Make Data Accesible
    % dataBoth = [y,  X1.....Xns,  u,  Xd]
    y_exp      = timeseries();
    u_exp      = timeseries();
    
   
    
    y_exp.Data = data(:,1);
    ty = (0:Ts:(length(y_exp.Data)-1)*Ts)';
    y_exp.Time = ty;
    
%     keyboard
    u_exp.Data = data(:, end);
    tu = (0:Ts:(length(u_exp.Data) - 1)*Ts)';
    u_exp.Time = tu;
    
    
    
end