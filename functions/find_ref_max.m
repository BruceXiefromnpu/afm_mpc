% ref_max = find_ref_max(sim_struct, ref_range, ref_step)
%
% This function runs a (saturated) LQR or MPC simlation with the same
% control parameters over many different setpoints to determine the largest
% achievable setpoint (before instability results due to saturation).
% 
% It is designed to run the MPC_fp simulink model.
%
% Inputs
% ------
%   sim_struct: a struct containing the data simulink needs to run. It must
%               contain the following fields
%
%     sim_struct.PLANT :  the LTI plant model
%     sim_struct.Nx    :  normalized Xss
%     sim_struct.K_lqr :  LQR gain (equivilent gain if doing to MPC
%     simulation)
%     sim_struct.mpcProb1 : an instance of sparseMPCprob
%     sim_struct.trun : how long the simulation should for
%     sim_struct.du_max : the slew rate saturation limit
%     sim_struct.mpc_on : boolean flag of weather to run linear or MPC sim.
%
%  ref_range: an array with two elements 
%              ref_range = [ref_f_lower, ref_f_upper]
%  ref_step: the increment of the reference. All values in 
%  ref_f_lower:ref_step:ref_f_upper will be tested.
%
% Option Inputs
% ------------
%   find_ref_max(..., 'verbose', int) if int=0, no plotting is down. If
%   int>0, will plot the y(t), u(t), and du(t) for each tested setpoint
%   into the same set of three figures. Default is verbose=0.
% 
%  find_ref_max(..., 'fig_base', int). Only used if verbose >=1. Allows the
%  figures generated to be incremented. The scheme is 
%  figure(1 + 10*fig_base) etc. Default is fig_base=10;
%
% Outputs
% -------
%   ref_max:  The maximim stable setpoint that was found. If the very first
%   iteration is unstable, ref_max=0.
%
%
% See Also
% --------
%   sparseMPCprob 

function [ref_max, t_settle_s] = find_ref_max(sim_struct, ref_s, varargin)
    % Auxiliary functions for checking parameter values:
    isnonnegint = @(x) (isnumeric(x) && mod(x, 1) == 0 && x >= 0);
    
    % Parse options
    defaultVerbose = 0;
    defaultFigbase = 10;
    p = inputParser;
    addParameter(p, 'verbose', defaultVerbose, isnonnegint)
    addParameter(p, 'fig_base', defaultFigbase, isnonnegint)
    parse(p, varargin{:})
    verbose = p.Results.verbose;
    fig_base = p.Results.fig_base;
    
    % ------------------------------------------------------------------- %
    % Pull out all the data stored in sim_struct to expose it to simulink.
    PLANT = sim_struct.PLANT;
    Ts = PLANT.Ts;
    Nx = sim_struct.Nx;
    
    K_lqr = sim_struct.K_lqr;
    mpcProb1 = sim_struct.mpcProb1;
    trun     = sim_struct.trun;
    du_max   = sim_struct.du_max;
    mpc_on   = sim_struct.mpc_on;


    x0 = sim_struct.xss*0;
    uss_0 = 0;

    figs = [];
    if verbose > 0
         F1 = figure(1 + 10*fig_base);
         figs = [figs, F1];
    end
    if verbose > 1
        F2 = figure(2 + 10*fig_base);
        F3 = figure(3 + 10*fig_base);
        figs = [figs, F2, F3];
    end
    
    % We have to set the simulink workspace to this function, since
    % default is the main scripts workspace.
    options = simset('SrcWorkspace','current');
    
    t_settle_s = zeros(1, length(ref_s));
    for iter = 1:length(ref_s)
        ref_f = ref_s(iter);
        sim('MPC_fp', [], options)
        y1 = y_mpcDist; 
        [t_settle, k_s] = settle_time(y1.time, y1.Data, ref_f, 0.01*ref_f,...
                                      [], [], 30);

        if verbose
            figs = plot_local(y1, u_mpcDist, du_mpcDist, figs, fig_base, verbose);
        end
        
        if isnan(t_settle)
            break
        else
            t_settle_s(iter) = t_settle;
        end


    end 
    
    if iter > 1
        ref_max = ref_s(iter-1);
    else
        ref_max = 0;
    end

end


function figs = plot_local(y1, u1, du1, figs, fig_base, verbose)

    % Figure (1)
    if verbose > 0
        try
            change_current_figure(figs(1))
        catch
            figs(1) = figure(1 + 10*fig_base);
        end
        hold on
        plot(y1.Time, y1.Data)
        title('y(t)')
    end
    if verbose > 1
        % Figure 2
        try
            change_current_figure(figs(2))
        catch
            figs(2) = figure(2 + 10*fig_base); 
        end
        hold on
        plot(u1.Time, u1.Data)
        title('u(t)')

        % Figure 3
        try
            change_current_figure(figs(3))
        catch
            figs(3) = figure(3 + 10*fig_base); 
        end
        hold on
        plot(du1.Time, du1.Data)
        title('du(t)')
    end
    drawnow()

        
    end








