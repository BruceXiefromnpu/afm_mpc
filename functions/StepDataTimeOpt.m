classdef StepDataTimeOpt < StepData
    
    properties
        % params;
        % file;
        % fig_files;
        % savedata;
        % results;
    end
    
    methods
        function self = StepDataTimeOpt(Params, varargin)
            self = self@StepData(Params, varargin{:})
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
        
        function[h, ax] = plot_ref_vs_settle(self,ax, varargin)
        % plot_ref_vs_settle(self,ax, varargin)
        % plot reference vs settle time for data contained in
        % self.data. 
        %   -- If ax is empty, will plot to gca().
        %   -- varargin is passed straight to matlabs plot function. 
            if ~exist('ax', 'var')
                ax = gca();
            elseif isempty(ax)
                ax = gca();
            end
            ref_s = self.params.ref_s;
            topt_settletime_s = self.results.settle_times_opt_cell{1};
            h = plot(ax, ref_s, topt_settletime_s*1000, varargin{:});
            set(h, 'DisplayName', 'TimeOpt')
        end
        
        function [hy, ax] = plot_single_ytraj(self, index, ax, ...
                                             varargin)
        % plot the y-trajectory held at 
        % self.results.opt_trajs_cell{1}.Y_vec_s{index}.
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

        function [hu, ax] = plot_single_utraj(self, index, ax, varargin)
        % Plot the u-trajectory held at 
        % self.results.opt_trajs_cell{1}.U_vec_s{index}.
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

        function [self, status] = build_timeopt_trajs(self, varargin)
        % Builds Time-Optimal trajectories Will run a sequence
        % of simulations over a list of references, performing a
        % bisection search for the time optimal trajectory on each
        % reference Will record both the resulting trajectories as 
        % well as the associated settling times. 
        %
        % The big function of this function is to check if data_struct.file
        % exists. If it does, then the function compares the data saved in that
        % .mat file to the data provided in data_struct. If the data is the same,
        % the function simply loads the .mat file and returns the structure
        % contained therin. If the data has changed, the function will perform
        % the simulations again, saving the new data into the .mat file specified
        % in data_struct.file.
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
        %   build_timeopt_trajs(..., 'do_eject' (true|false))
        %   Wheather or not to eject the real pole zero pair at ~100
        %   Hz. I should augment this to also allow inputing an overloading
        %   function handle that will do this. 
        % 
        %   build_timeopt_trajs(..., 'verbose', 1) 1 --> write
        %   logging info to self.logger. 2 --> make plot of settle times at
        %   figure 200 also. 3--> also plot all the resulting Y and U 
        %   trajectories to figure 500.
        % 
        % Outputs
        % ------
        %   On exit, this function will populate self.results.

            defaultForce = 0;
            p = inputParser;
            p.addParameter('force', defaultForce);
            p.addParameter('max_iter', 20);
            p.addParameter('do_eject', true);
            p.addParameter('verbose', 1);
            p.addParameter('savedata', 1);

            parse(p, varargin{:});

            force = p.Results.force;
            max_iter = p.Results.max_iter;
            do_eject = p.Results.do_eject;
            verbose = p.Results.verbose;
            savedata = p.Results.savedata;
            
            status = 0;
            self.logger('%s LOG (build_timeopt_trajs):\n', ...
                        datestr(now));
            if force
                self.logger(['Force = true, rebuilding time-optimal ' ...
                             'trajectories\n']);
            elseif self.stepdata_struct_unchanged() 
                other = load(self.file);
                self.logger(['Data appears to be the same. Loading data ',...
                              'without re-calculation.\n\n']);
                self = other.step_data;
                return
            else
            self.logger(['data has changed: re-building max ' ...
                         'setpoints.\n']);
            end
            
            params = self.params;
            
            % Expose parameters:
            ref_s = params.ref_s;
            du_max = params.du_max;
            
            sys_nodelay = params.sys_nodelay;
            Nd = params.Nd;
            
            
            if verbose >=2
                Fig = figure(200);
                hands = [];
                change_current_figure(Fig);
                colrs = get(gca, 'colororder');
            end
            if verbose >=2
                Fig2 = figure(500);
            end

            if do_eject
                % Pull out open-loop pole-zero information.
                sys_sim = eject_gdrift(params.sys_nodelay);
                sys_sim = SSTools.deltaUkSys(sys_sim);
            else
                sys_sim = SSTools.deltaUkSys(sys_nodelay);
            end
            
            [Nx_sim, Nu_sim] = SSTools.getNxNu(sys_sim);
            x0_sim = Nx_sim*0;
            
            % Pre-allocate
            opt_trajs_cell = cell(1,1);
            settle_times_opt = ref_s*0;
            
            N_refs = length(ref_s);
            Y_vec_s = cell(1, N_refs);
            U_vec_s = cell(1, N_refs);
            X_vec_s = cell(1, N_refs);
            
            
            warning('off', 'MATLAB:nargchk:deprecated')
            toBisect = TimeOptBisect(sys_sim, du_max);
            toBisect.max_iter = max_iter;
            k0 = size(sys_sim.B,1);
            
            for iter = 1:length(ref_s)
                ref_f = ref_s(iter);
                xf = Nx_sim*ref_s(iter);
                [X, U, status]=toBisect.time_opt_bisect(x0_sim, xf, 'k0', k0);
                if status
                    self.logger(['Time-optimal bisection failed at ref = ' ...
                                  '%0.3f. Exiting...'], ref_f);
                    break
                else
                    self.logger(['Time Optimal bisection for ref=%.3f, ref_iter = %.0f ',...
                              'of %.0f, N*=%.0f\n\n'], ref_f,...
                              iter, length(ref_s), length(U.Time));
                end
                k0 = length(U.Time);
                % Extend the trajectory with steady state u values
                u = [U.Data; ones(2*k0,1)*ref_s(iter)*Nu_sim]; % zero for deltaUk
                t = [0:1:length(u)-1]'*sys_sim.Ts;
                sys_sim.InputDelay = self.params.sys.InputDelay;
                [y, t, x] = lsim(self.params.sys, u, t, x0_sim);
                sys_sim.InputDelay = 0;
                
                U_vec_s{iter} = timeseries(u, t);
                X_vec_s{iter} = timeseries(x, t);
                Y_vec_s{iter} = timeseries(y, t);
                settle_times_opt(iter) = settle_time(t, y, ref_f, 0.01*ref_f,...
                                             [], [], 30);
                
                %U.Time(end) + Nd*sys_sim.Ts;

                if verbose >=2
                    change_current_figure(Fig); 
                    hold on
                    hands(iter) = plot(ref_s(iter), settle_times_opt(iter)*1000, '-xk'); 
                    drawnow()
                end
                if verbose >=3
                    change_current_figure(Fig2);
                    if mod(iter,20) ==0
                        clf;
                    end
                    subplot(2,1,1); 
                    plot(t, y); hold on
                    xlim([0, 0.0320])
                    subplot(2,1,2)
                    plot(t, u); hold on
                    xlim([0, 0.0320])
                    drawnow
                end
            end
            
            % If the bisection failed, truncate the results, so that the
            % simulation will run again.
            if status
                self.params.ref_s = ref_s(1:iter-1);
                settle_times_opt = settle_times_opt(1:iter-1);
                Y_vec_s = Y_vec_s{1:iter};
                U_vec_s = U_vec_s{1:iter};
                X_vec_s = X_vec_s{1:iter};
            end
            opt_trajs_cell{1}.Y_vec_s = Y_vec_s;
            opt_trajs_cell{1}.X_vec_s = X_vec_s;
            opt_trajs_cell{1}.U_vec_s = U_vec_s;
            
            self.results.opt_trajs_cell = opt_trajs_cell;
            self.results.settle_times_opt_cell{1} =  settle_times_opt;
            
            if savedata
                step_data = self;
                save(self.file, 'step_data')
            end
            
            
        end % end main function
    end %METHODS
end % CLASSDEF




function sys_vibrational = eject_gdrift(sys_nodelay)
    Ts = sys_nodelay.Ts;
    [wp_real_x, wz_real_x] = w_zp_real(sys_nodelay);
%     rho_1 = wz_real_x(1)/wp_real_x(1);
    g_eject = zpk(exp(-wp_real_x(1)*Ts), exp(-wz_real_x(1)*Ts), 1, Ts);
    g_eject = g_eject/dcgain(g_eject);
    sys_vibrational = minreal(sys_nodelay*g_eject);
end
