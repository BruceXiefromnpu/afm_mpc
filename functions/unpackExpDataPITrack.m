% [y_exp, u_exp, ypow_exp] = unpackExpData(data, Ts)
% 
% Unpacks the experimental data that LabView spits out into pieces that
% make sense.
%
% Inputs: 
%   data: from labview. We expect the following format:
%       dataBoth = [y,  X1.....Xns,  u,  Xd]
%   Ts: sample time
% Outputs:
%   y_exp, u_exp, ypow_exp. Actual measured input/output to AFM.
% 
%
%   These are all timeseries objects.


function [y_exp, u_exp, ypow_exp] = unpackExpDataPITrack(data,  Ts)
    % Make Data Accesible
    % dataBoth = [y,  X1.....Xns,  u,  Xd]
    y_exp      = timeseries();
    u_exp      = timeseries();
    
   
    
    y = data(:,1);
    t = (0:Ts:(length(y)-1)*Ts)';
    y_exp = timeseries(y, t);
    
%     keyboard
    u = data(:, 2);
    u_exp = timeseries(u, t);
    
    ypow = data(:,3);
    ypow_exp = timeseries(ypow, t);

    
    
    
    
end