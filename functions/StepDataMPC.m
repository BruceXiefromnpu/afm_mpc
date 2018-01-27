classdef StepDataMPC < StepDataQuad
    
    properties
        % params;
        % file;
        % fig_files;
        % savedata;
        % results;
    end
    
    methods
        function self = StepDataMPC(Params, varargin)
            self = self@StepDataQuad(Params, varargin{:})
        end
        function self = build_max_setpoints(self, varargin)
        % self = build_max_setpoints(self, varargin) 
        %    
        % Builds a set of maximum setpoints for a family of LQR based
        % linear or MPC based controllers. The controllers are
        % parameterized by gamma. 
        %
        % Outputs
        % -------
        %   step_data : the function will populate the StepData.results
        %   property upon exit. In addition, if the StepData.file field has
        %   been set, will also save the resulting instance into that file
        %   location.
        % 
        %   Upon (succseful exit), the results property will contain the
        %   following fields for the linear case:
        %   results.max_setpoints.
        %          results.data{idx}.t_settle_s
        %          results.data{idx}.y_traje_s
        %          results.data{idx}.ref_max
        %   where idx is the index for gamma.
        % 
        %   For the mpc case, results is itself a cell array, with size=length(N_mpc).
        % 
        % Optional Inputs
        % --------------
        %   build_max_setpoints(..., 'force', (true|false)) if true,  will
        %   force the simulations to be re-run regardless of if the
        %   underlying parameters are the same as the saved ones.
        % 
        %   build_max_setpoints(..., 'fid', fid) Default fid=1. Pass in a
        %   file id from fopen to write all logging info to a file. Will
        %   only be used if step_data.verbose >0.
        %
        %   verbose = 0 --> only initial and final status printed to console.
        %   verbose =  1 print progress bar to console
        %   verbose =  2 plot gamma vs max in Fig 1, and settle-times vs
        %   ref in Fig (all gammas).
        %   verbose = 3 plot y-trajs into a figure. A new figure is
        %   generated for each gamma.
        %   verbose = 4 also plot u/du-trajs into 2 figures. A new figure is
        %   generated for each gamma. This generates a lot of figures, so
        %   currently, they are all deleted after every 10th gamma.
        % 
        % Algorithm/Methods
        % ----------------
        %    Iterate over a bunch of gammas and try to find the maximum setpoint for
        %    each one. This code is pretty niave. We start at a very low setpoint and
        %    slowly increase the setpoint until the settling time is reported as NaN.
        %    A bisection search would probably be faster. This is complicated though
        %    by the fact that instability seems to occur in the "middle": if the
        %    setpoint is large enough, we dont have the stability problem, which is
        %    weird.
        %
        %  MPC Case 
        %  --------
        %    This works similarly, but now we also allow the simulations to
        %    be done over a family of control horizons. Thus, there are two
        %    loops: the outer loop loops of N_mpc_s, the inner loop loops
        %    over gammas, like in the linear case.
            
            
        % pull data out of data variable:
            defaultForce = 0;
            p = inputParser;
            p.addParameter('force', defaultForce);
            p.addParameter('verbose', 1);
            parse(p, varargin{:});
            force = p.Results.force;
            
            verbose = p.Results.verbose;

            self.logger('%s LOG (build_max_setpoints, mpc_on=%d)\n', ...
                        datestr(now), self.params.sim_struct.mpc_on);
            if force
                self.logger('Force flag = true: re-building max setpoints.\n');
            elseif self.stepdata_struct_unchanged()
                other = load(self.file);
                self.results = other.step_data.results;
                self.logger(['Data appears to be the same. Loading data ',...
                        'without re-calculation.\n']);
                return
            else
                self.logger(['data has changed: re-building max ' ...
                             'setpoints.\n']);
            end
            
            
            % ------------------------------------------------------------------ %
            % The parameter data has changed, so we re-run. Setup
            % and run the simulations 

            %        -Where max_setpoints is a vector (list) of the maximum
            %         setpoint achievable for each supplied gamma.
            %        - data is a vector of structs containing:
            %            data(idx).t_settle_s
            %            data(idx).y_traj_s
            %            data(idx).ref_max
            
            N_mpc_s = self.params.N_mpc_s;
            gam_s = self.params.gam_s;
            result_s = cell(1, length(self.params.N_mpc_s));
            sim_struct.mpc_on = true;
            ProgBar = @self.ProgBar;
            warning('OFF', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
            parfor mpc_iter = 1:length(N_mpc_s)
                warning('OFF', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
                N_mpc_iter = N_mpc_s(mpc_iter);
                if verbose > 0
                    start_str = sprintf('MPC, N=%.0f', N_mpc_iter);
                    PB = ProgBar(length(gam_s), 'start_str', ...
                                 start_str);
                end
%                 self.params.sim_struct.N_mpc = N_mpc_iter;
                result_s_iter = self.build_max_sp_local(PB, verbose, ...
                                                   N_mpc_iter);
                result_s{mpc_iter} = result_s_iter;
            end
            % warning('ON', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
            % results.max_setpoints = max_setpoints;
            self.results = result_s;

            if self.savedata
                try
                    step_data = self;
                    save(step_data.file, 'step_data', '-v7.3')
                catch ME
                    errMsg = getReport(ME,  'extended','hyperlinks', 'off');
                    self.logger('*********** DATA NOT SAVED ************* \n');
                    self.logger('Matlab said:\n %s', errMsg);
                end
            end
            
        end % END MAIN FUNCTION
        
    end
    
end


% -----------------------------------------------------------------------%
% ----------------------- Utility Functions -----------------------------%

function plot_local(t_settle_s, gam_s, ref_s, max_setpoints, gam_iter)
% expose for easy access while plotting:
    
    if mod(gam_iter, 10)==0
        self.logger('\ndeleting current figures\n\n')
        for i=1:gam_iter
            close(figure(10*i + 1))
            close(figure(10*i + 2))
            close(figure(10*i + 3))
        end
    end
    F1 = figure(1); hold on;
    xlabel('$\gamma$', 'interpreter', 'latex', 'FontSize', 16)
    ylabel('Max Ref', 'interpreter', 'latex', 'FontSize', 16)
    
    F2 = figure(2);  hold on;
    ylabel('settle time [ms]')
    xlabel('ref_f')

    change_current_figure(F2);
    % colrs = get(gca, 'colororder');
    k_max = find(t_settle_s ~= 0, 1, 'last');
    
    plot(ref_s(1:k_max), t_settle_s(1:k_max)*1e3)
    % ylim([0, 10])
    hold on
    drawnow()

    change_current_figure(F1)
    hlin = plot(gam_s(1:gam_iter), max_setpoints(1:gam_iter),...
                '-o',  'LineWidth', 2);

    drawnow()
    hold on
end



