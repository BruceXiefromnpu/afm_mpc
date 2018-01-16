classdef StepDataCLQR < StepData
    
    properties
        % params;
        % file;
        % fig_files;
        % savedata;
        % results;
        plot_ref_vs_settle_fcn;
    end
    
    methods
        function self = StepDataCLQR(Params, varargin)
        % % First, parse the inputs we care about here:
        %     p=inputParser();
        %     p.KeepUnmatched = true;
        %     p.addParameter('logger', @fprintf);
        %     p.addParameter('ProgBar', @ProgressBar);
        %     parse(p, varargin{:});
        %     % Now, pull out everything that doesn't match, and
        %     % re-package it into a cell array (varargin like) to
        %     % pass into StepData. p.Unmatched is a struct. 
        %     unm  = p.Unmatched;
        %     unm_flds = fieldnames(unm);
        %     varg = {};
        %     for k = 1:length(unm_flds)
        %         varg{end+1} = unm_flds{k};
        %         varg{end+1} = unm.(unm_flds{k});
        %     end
            self = self@StepData(Params, varargin{:});
            

            
            self.plot_ref_vs_settle_fcn = @plot_ref_vs_settle;

            % p = inputParser;
            % p.addParameter('savedata', true)
            % p.addParameter('file', '')
            % p.addParameter('fig_files', '')
            % parse(p, varargin{:});
            
            % self.params = Params;
            % self.file = p.Results.file;
            % self.fig_files = p.Results.fig_files;
            % self.savedata =  p.Results.savedata;
            % self.results = [];
        end
        function [h, ax] = plot_ref_vs_settle(self, ax, gamma, varargin)
            [h, ax]= self.plot_ref_vs_settle_fcn(self, ax, gamma, varargin{:});
        end
        
        function [h, ax] = plot_single_ytraj(self, index, ax, varargin)
        % plot the y-trajectory held at self.results.opt_trajs_cell{1}.Y_vec_s{index}
        %   -- If ax is empty, will plot to gca().
        %   -- varargin is passed straight to matlabs plot function. 

            if ~exist('ax', 'var')
                ax = gca();
            elseif isempty(ax)
                ax = gca();
            end
            traj_y = self.results.opt_trajs_cell{1}.Y_vec_s{index};
            hy = plot(ax, traj_y.Time, traj_y.Data, varargin{:});
            
        end

        function [h, ax] = plot_single_utraj(self, index, ax, varargin)
        % Plot the u-trajectory held at self.results.opt_trajs_cell{1}.U_vec_s{index}
        %   -- If ax is empty, will plot to gca().
        %   -- varargin is passed straight to matlabs plot function. 
            if ~exist('ax', 'var')
                ax = gca();
            elseif isempty(ax)
                ax = gca();
            end

            traj_u = self.results.opt_trajs_cell{1}.U_vec_s{index};
            hu = plot(ax, traj_u.Time, traj_u.Data, varargin{:});
            
        end
        
        function self = build_clqr_trajs(self, varargin)
        % clqr_data = build_clqr_trajs(step_data, varargin)
        %
        % Builds CLQR trajectories parameteriezed by gamma. Will run a sequence
        % of simulations over a list of references, for each gamma. Will record
        % both the resulting trajectories as well as the associated settling
        % times. 
        %
        % Inputs
        % -----
        %   step_data : a class instance of StepData 
        % 
        % Optional Inputs
        % --------------
        %   build_timeopt_trajs(..., 'force', (true|false)) Force a
        %   rerun of the simulations, regardless of wheather or not
        %   the parameters have changed.
        % 
        %   build_timeopt_trajs(..., 'max_iter', 20) Maximum number
        %   of iterations allowed to find an upper bound in the
        %   bisection search.
        % 
        %   build_timeopt_trajs(..., 'fid' 1) A file id to write
        %   logging data to.
        % 
        %   build_timeopt_trajs(..., 'verbose', 1) 1 --> write
        %   logging info to file id. 2 --> make plot to figure 200 also.
        % 
        %   build_timeopt_trajs(..., 'savedata', 1) save the
        %   results and the whole object to self.file? Will be saved
        %   as variable 'step_data'
        %
        % Outputs
        %  -----
        %  step_data : the same class that was passed in, but with the
        %  results field populated. The results field shall have the form:
        %  
        %    N.B. The results fields are cell arrays, to account for the
        %    occasion that we gam_s is a vector.
        % 
        %    results.settle_times_opt_cell : settle times for each reference.
        % 
        %    results.opt_trajs_cell        : the trajectories associated
        %    with each reference (Y, U, X).
        %
        %
        % The big function of this function is to check if step_data.file
        % exists. If it does, then the function compares the data saved in that
        % .mat file to the data provided in data_struct. If the data is the same,
        % the function simply loads the .mat file and returns the structure
        % contained therin. If the data has changed, the function will perform
        % the simulations again, saving the new data into the .mat file specified
        % in data_struct.file.
        %
        %
        % See Also: opt_traj_gen, StepData

            defaultForce = 0;
            p = inputParser;
            p.addParameter('force', defaultForce);
            p.addParameter('verbose', 1);
            parse(p, varargin{:});
            force = p.Results.force;
            verbose = p.Results.verbose;
            
            if self.stepdata_struct_unchanged() && ~force
                other=load(self.file); % should provide step_data
                self.logger('LOG: (build_clqr_trajs)\n');
                self.logger(['Data appears to be the same. Loading data ',...
                              'without re-calculation.\n\n']);
                self = other.step_data;
                return
            end
            self.logger('LOG (build_clqr_trajs\n');
            self.logger(['Data has changed or force=true: re-building ' ...
                          'CLQR trajectories.\n']);
            
            gam_s = self.params.gam_s;
            % Pre-allocate
            % settle_times_opt_save = cell(1, length(gam_s));
            self.results.settle_times_opt_cell = cell(1, length(gam_s));
            % opt_trajs_save = cell(1, length(gam_s)); 
            self.results.opt_trajs_cell = cell(1, length(gam_s));

            
            if verbose >= 2
                Fig = figure(200);
                ax = gca();
                hands = [];
                change_current_figure(Fig);
            end
            ref_0 = 0;
            for iter = 1:length(gam_s)
                gamma_k = gam_s(iter);

                [traj_s, settle_times_opt] = self.opt_traj_gen(gamma_k, ref_0, 'verbose', verbose,...
                                                          'mpc_mode', self.params.mpc_mode);
                % settle_times_opt_save{iter} = settle_times_opt;
                % opt_trajs_save{iter} = traj_s;
                self.results.settle_times_opt_cell{iter} = settle_times_opt;
                self.results.opt_trajs_cell{iter} = traj_s;
                if verbose >=2
                    self.plot_ref_vs_settle(ax, gamma_k, 'LineWidth', 2);
                end
            end
            % results.settle_times_opt_cell =  settle_times_opt_save;
            % results.opt_trajs_cell = opt_trajs_save;
            % self.results = results;
            
            step_data = self;
            if self.savedata
                save(self.file, 'step_data')
            end
            % end build_CLQR_trajs
        end 

        function [traj_s, settle_times] = opt_traj_gen(self, R, ref_0_s, varargin )
        % [traj_s, settle_times] = opt_traj_gen(Q, R, N_traj, sys, ref_s, ref_0_s, varargin )
        % The goal of this function is to generate a set of optimal trajectories
        % over a LONG horizon (for a sequence of setpoints). The immediate goal
        % here is to look at what is the minimum acheivable settling time for a set
        % of setpoints, given a (Q,R) pair, laying aside all issues about the MPC
        % control horizon.
        % 
        % Required Inputs
        % ------
        %  Q, R: state and control weighting matrices
        %  
        %  N_traj : the TOTAL length of the optimal trajectory. Eg., 400 is
        %  usuall reasonable for the AFM
        % 
        %  sys : discrete time dynamical system. The optimal control problem will
        %  be generated based on sys.
        %
        %  ref_s : a list of reference inputs to generate trajectories for.
        %  
        %  ref_0_s : a list of starting points. We will assume the system is at
        %  steady for y = ref_0(iter). Can be either a single point, in which case
        %  the same ref_0 is used for all trajectories, or can be a list the same
        %  length as ref_s.
        %
        %  Optional Inputs
        %  ---------------
        %    opt_traj_gen(..., 'S', S) : cross weighting matrix for the optimal
        %    control problem
        %
        %    opt_traj_gen(..., verbose, 1) (default 0) flag indiciating what plots to make
        %        verbose = 1, plot output
        %        verbose = 2, plot output and control
        %        verbose = 3, plot output, control and accumulated control (useful
        %        if doing optimization over incremental form)
        %
        %   opt_traj_gen(..., 'mpc_mode', {condensed|sparse}) specifies if
        %   we should use the sparse or condensed CLQR (mpc) formulation.
        %   If 'sparse', will use the sparseMPCprob class, if 'condensed',
        %   will use the condensedMPCprob class. 
        
            % ----------- Parse and Sanitize inputs ------------- %
                
            sys = self.params.sys;
            p = inputParser;
            
            defaultS = zeros(size(sys.b, 1), size(sys.b, 2));
            
            addParameter(p, 'S', defaultS);
            addParameter(p, 'verbose', 0);
            addParameter(p, 'mpc_mode', 'condensed');
            
            parse(p, varargin{:})
            S = p.Results.S;
            verbose = p.Results.verbose;
            mpc_mode = p.Results.mpc_mode;
            
            % Expose parameters:
            params = self.params;
   
            gam_s = params.gam_s;
            ref_s = params.ref_s;
            du_max = params.du_max;
            N_traj = params.N_traj;
            mpc_mode = params.mpc_mode;
            Q = params.Q;
            
            if length(ref_0_s) == 1
                ref_0_s = 0*ref_s + ref_0_s;
            elseif length(ref_0_s) ~= length(ref_s)
                error('length(ref_0_s) must be either one or the same as length(ref_s)')
            end
            
            Qp = dare(sys.a, sys.b, Q, R);
            if strcmpi(mpc_mode, 'sparse')
                NLQR_prob = sparseMPCprob(sys, N_traj, Q, Qp, R, S);
            elseif strcmpi(mpc_mode, 'condensed')
                NLQR_prob = condensedMPCprob(sys, N_traj, Q, Qp, R, S);
            end

            if du_max ~= 0 
                NLQR_prob.add_U_constraint('slew', du_max);
            end

            if verbose >= 3
                %     mkfig(verbose)
                Fig = figure();
            else
                Fig = 0;
            end

            Nx = SSTools.getNxNu(sys);
            
            N_refs = length(ref_s);
            Y_vec_s = cell(1, N_refs);
            U_vec_s = cell(1, N_refs);
            X_vec_s = cell(1, N_refs);
            
            tvec = (0:1:N_traj-1)*sys.Ts; 
            settle_times = zeros(1, length(ref_s));
            start_str = sprintf('CLQR, gamma: %.0f', R);
            upd = self.ProgBar(length(ref_s), 'start_str', start_str);
            
            upd.upd(0);
            for iter = 1:length(ref_s)
                ref_f = ref_s(iter);
                ref_0 = ref_0_s(iter);
                
                x0_err = ref_0*Nx - ref_f*Nx;
                
                [U, Xerr] = NLQR_prob.solve(x0_err, 'getX', true); 
                
                % Convert the solution trajectory from error coords to
                % regular coords.
                % keyboard
                X = Xerr + ref_f*Nx;
                Y = sys.c*X;
                t_settle = settle_time(tvec, Y, ref_f, 0.01*ref_f,...
                                             [], [], 30);
                settle_times(iter) = t_settle;
                % We have to transpose to columns, otherwise matlab will create a 3d
                % matrix out of this stuff, for some bizare reason.
                Uvec = timeseries(U', tvec);
                
                % We have one extra element N+1 states, but only N controls. So drop
                % the last one. 
                Xvec = timeseries(X(:,1:end-1)', tvec);
                Yvec = timeseries(Y(:,1:end-1)', tvec);
                
                U_vec_s{iter} = Uvec;
                X_vec_s{iter} = Xvec;
                Y_vec_s{iter} = Yvec;
                
                if verbose >= 3
                    Fig = plot_local(Yvec, Uvec, verbose, Fig);
                end
                
                upd.upd(iter);
                
            end
            % self.logger('PUT\n')
            traj_s.X_vec_s = X_vec_s;
            traj_s.Y_vec_s = Y_vec_s;
            traj_s.U_vec_s = U_vec_s;

        end            
        
    end
    % end methods    
end
% end classdef


function Fig= plot_local(Yvec, Uvec, verbose, Fig)
% Local plotting 
    
    if ~isvalid(Fig)
        Fig = figure();
        
    end % if ~isvalid(Fig)
        
    change_current_figure(Fig);
        
    if verbose == 1
        subplot(111); hold on
        plot(Yvec.time, Yvec.data)
        ylabel('y(k)')
        xlabel('time')
        drawnow()
        
    elseif verbose ==2
        subplot(211); hold on
        plot(Yvec.time, Yvec.data)
        ylabel('y(k)')
        
        subplot(212); hold on
        plot(Uvec.time, Uvec.data)
        ylabel('u(k)')
        xlabel('time')
        drawnow()
    elseif verbose == 3
        subplot(311); hold on
        plot(Yvec.time, Yvec.data)
        ylabel('y(k)')
        
        subplot(312); hold on
        plot(Uvec.time, Uvec.data)
        ylabel('u(k)')
        
        subplot(313); hold on
        plot(Uvec.time, cumsum(Uvec.data))
        ylabel('accum(u(k))')
        xlabel('time')
        drawnow()    
    end
end


function [h, ax] = plot_ref_vs_settle(self,ax,gam, varargin)
% plot_ref_vs_settle(self,ax, varargin)
% plot reference vs settle time for data contained in
% self.data. 
%   -- If ax is empty, will plot to gca().
%   -- varargin is passed straight to matlabs plot function. 
    if ~exist('ax', 'var') && isvalid(ax)
        ax = gca();
    elseif isempty(ax)
        figure()
        ax = gca();
    end
    ref_s = self.params.ref_s;
    clqr_settletime_s = self.results.settle_times_opt_cell{1};
    
    h = plot(ax, ref_s, clqr_settletime_s*1000, varargin{:});
    set(h, 'DisplayName', sprintf('$\\gamma = %.0f$', gam));
                     
     leg = legend(h);
     set(leg, 'interpreter', 'latex');
    ylabel('settle time [ms]', 'FontSize', 16)
    xlabel('setpoint', 'FontSize', 16)
    drawnow()
    grid on
    
end