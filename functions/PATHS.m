classdef PATHS
    %PATHS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function [ PATH ] = exp()
            % PATH constant to where all experimental data is stored for the
            % MPC journal paper.
            if ispc
                PATH = 'Z:\mpc-journal';
            else
                PATH = '/media/labserver/mpc-journal';
            end
        end
        
        function [ PATH ] = MPCJ_root()
            % PATH constant to where all experimental data is stored for the 
            % MPC journal paper.

            PATH = fullfile(getMatPath, 'afm_mpc_journal');
        end
        
        function [ PATH ] = sim_models()
            % PATH constant to where all experimental data is stored for the
            % MPC journal paper.
            
            PATH = fullfile(PATH_MPCJ_root, 'models');
        end
        
        function [ PATH ] = step_exp()
            %PATH_EXP Summary of this function goes here
            %   Detailed explanation goes here
            
            PATH = fullfile(PATH_experimental, 'step-exps');
            
            
        end
        
        function [ PATH ] = sysid()
            % PATH constant to where the system ID data is stored
            
            PATH = fullfile(PATH_experimental, 'sysID');
        end

    end
    
end

