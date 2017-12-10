function [ PATH ] = PATH_experimental()
% PATH constant to where all experimental data is stored for the 
% MPC journal paper.

if ispc
    PATH = 'C:\Users\arnold\Documents\Matlab';
else
    PATH = '/home/arnold/matlab/afm_mpc_journal';
end

end

