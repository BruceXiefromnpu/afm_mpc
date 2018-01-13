% ref_max = find_ref_max(sim_struct, ref_s, varargin.)
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
%     ref_s : a vector of references to check. Should be in
%     increasing order, as they will be tested from smallest to largest.
%
% Option Inputs
% ------------
%   find_ref_max(..., 'verbose', int) if int<=2, no plotting is down. If
%   int>2, will plot the y(t). If int > 3, will also plot for
%   u(t), and du(t) for each tested setpoint into the same set of
%   three figures. Default is verbose=0.
% 
%  find_ref_max(..., 'fig_base', int). Only used if verbose >2. Allows the
%  figures generated to be incremented. The scheme is 
%  figure(1 + fig_base) etc. Default is fig_base=10;
%
% Outputs
% -------
%   result.The maximim stable setpoint that was found. If the very first
%   iteration is unstable, ref_max=0.
%   result.ref_max : The maximum achievable reference.
% 
%   result.y_traj_s : A vector of time series trajectories from
%   ref_s(1) to ref_s(ref_max_ind).
% 
%   result.t_settle_s : A vector of all the settle times from
%   ref_s(1) to ref_s(ref_max_ind).
%
% See Also
% --------
%   sparseMPCprob, build_max_setpoints, StepData, StepParamsLin,
%   StepParamsMPC
%
% Depends on the simulink model 'MPC_fp'.

function result = find_ref_max(sim_struct, ref_s, varargin)
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
    

    figs = [];
    if verbose > 2
         F1 = figure(1 + fig_base);
         figs = [figs, F1];
    end
    if verbose > 3
        F2 = figure(2 + fig_base);
        F3 = figure(3 + fig_base);
        figs = [figs, F2, F3];
    end
    
    % We have to set the simulink workspace to this function, since
    % default is the main scripts workspace.
    % options = simset('SrcWorkspace','current');
    y_traj_s = repmat(timeseries(0,0), 1, length(ref_s));
    t_settle_s = zeros(1, length(ref_s));
    for iter = 1:length(ref_s)
        ref_f = ref_s(iter);
        % sim('MPC_fp', [], options)
        % y1 = y_mpcDist; 
        [Y, U, dU] = sim_local(sim_struct, ref_f);
        [t_settle, k_s] = settle_time(Y.time, Y.Data, ref_f, 0.01*ref_f,...
                                      [], [], 30);

        if verbose
            figs = plot_local(Y, U, dU, figs, fig_base, verbose);
        end
        
        if isnan(t_settle)
            break
        else
            t_settle_s(iter) = t_settle;
            y_traj_s(iter) = Y;
        end    
    end 
    
    if iter > 1
        ref_max = ref_s(iter-1);
        result.ref_max = ref_max;
        result.y_traj_s = y_traj_s;
        result.t_settle_s = t_settle_s;
    else
        result.ref_max = NaN;
        result.t_settle_s = [];
        result.y_traj_s = Y;
    end

end

function [Y, U, dU] = sim_local(sim_struct, ref_f)
% ------------------------------------------------------------------- %
% Pull out all the data stored in sim_struct to expose it to
% simulink. There must be a better way...

    options = simset('SrcWorkspace','current');
    % Expose the sime struct to simulink.
    K_lqr = sim_struct.K_lqr;
    PLANT = sim_struct.PLANT;
    trun = sim_struct.trun;
    mpcProb1 = sim_struct.mpcProb1;
    du_max = sim_struct.du_max;
    mpc_on = sim_struct.mpc_on;
    xss = sim_struct.xss;
    Nx = sim_struct.Nx;
    
    x0 = sim_struct.xss*0;
    uss_0 = 0;
    Ts = PLANT.Ts;
    sim('MPC_fp', [], options)
    
    Y = y_mpcDist;
    U = u_mpcDist;
    dU = du_mpcDist;
end


function figs = plot_local(y1, u1, du1, figs, fig_base, verbose)

    % Figure (1)
    if verbose > 2
        try
            change_current_figure(figs(1))
        catch
            figs(1) = figure(1 + fig_base);
        end
        hold on
        plot(y1.Time, y1.Data)
        title('y(t)')
    end
    if verbose > 3
        % Figure 2
        try
            change_current_figure(figs(2))
        catch
            figs(2) = figure(2 + fig_base); 
        end
        hold on
        plot(u1.Time, u1.Data)
        title('u(t)')

        % Figure 3
        try
            change_current_figure(figs(3))
        catch
            figs(3) = figure(3 + fig_base); 
        end
        hold on
        plot(du1.Time, du1.Data)
        title('du(t)')
    end
    drawnow()

        
end







