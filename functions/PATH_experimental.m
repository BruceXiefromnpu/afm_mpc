function [ PATH ] = PATH_experimental()
% PATH constant to where all experimental data is stored for the 
% MPC journal paper.

if ispc
    PATH = 'Z:\mpc-journal';
else
    PATH = '/media/labserver/mpc-journal';
end

end

