classdef StepDataLin < StepDataQuad
    
    properties
        % params;
        % file;
        % fig_files;
        % savedata;
        % results;
    end
    
    methods
        function self = StepDataLin(Params, varargin)
            self = self@StepDataQuad(Params, varargin{:})
        end
        
        function [hy, ax] = plot_single_ytraj(self, idx_ref, idx_gam, ax, varargin)
        % plot the y-trajectory held at self.results{1}.data{idx_gam}.y_traj_s(idx_ref);
        %   -- If ax is empty, will plot to gca().
        %   -- varargin is passed straight to matlabs plot function. 
            if ~exist('ax', 'var')
                ax = gca();
            elseif isempty(ax)
                ax = gca();
            end
            traj_y = self.results{1}.data{idx_gam}.y_traj_s(idx_ref);
            hy = plot(ax, traj_y.Time, traj_y.Data, varargin{:});
        end
        function [hu, ax] = plot_single_utraj(self, idx_ref, idx_gam, ax, varargin)
        % plot the u-trajectory held at self.results{1}.data{idx_gam}.u_traj_s(idx_ref);
        %   -- If ax is empty, will plot to gca().
        %   -- varargin is passed straight to matlabs plot function. 
            if ~exist('ax', 'var')
                ax = gca();
            elseif isempty(ax)
                ax = gca();
            end

            traj_du = self.results{1}.data{idx_gam}.du_traj_s(idx_ref);
            hu = plot(ax, traj_du.Time, traj_du.Data, varargin{:});
        end
        
        function plot_single_traj(self, idx_ref, idx_gam, ax1, ax2, varargin)
        % plot_single_traj(self, ref_idx, ax1, ax2, varargin)
        %
        % plot the trajectory for gam_s(idx_gam) and ref_s(idx_gam)
        % to ax1 (for y(k)) and ax2 for u(k).
            if ~exist('ax1', 'var') || isempty(ax1) || ~exist('ax2', 'var') || isempty(ax2)
                fig = gcf();
                ax1 = subplot(211);
                ax2 = subplot(212);
            end
            
            self.plot_single_ytraj(idx_ref, idx_gam, ax1, varargin{:});
            self.plot_single_utraj(idx_ref, idx_gam, ax2, varargin{:});
        end
        
        function self = build_max_setpoints(self, varargin)
        % step_data = build_max_setpoints(step_data, varargin) 
        %    
        % Builds a set of maximum setpoints for a family of LQR based
        % linear or MPC based controllers. The controllers are
        % parameterized by gamma. 
        %
        % Inputs
        % -----
        %   step_data : a class instance of StepData 
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
            
                    
        % This function just checks checks wheather we need to
        % re-run the simulation, then calls
        % build_max_setpoints_local in the parent class, StepDataQuad. 
            defaultForce = 0;
            p = inputParser;
            p.addParameter('force', defaultForce);
            p.addParameter('verbose', 1);
            parse(p, varargin{:});
            force = p.Results.force;
            
            verbose = p.Results.verbose;

            self.logger('%s LOG (build_max_setpoints, linear)\n', datestr(now));
            if force
                self.logger('Force flag = true: re-building max setpoints.\n');
            elseif stepdata_struct_unchanged(self)
                other = load(self.file);
                self.results = other.step_data.results;
                self.logger(['Data appears to be the same. Loading data ',...
                        'without re-calculation.\n']);
                return
            else
                self.logger(['data has changed or self.file does not',...
                             'exist or force flag=true : re-building',...
                             'max setpoints.\n']);
            end
            
            % The parameter data has changed, so we re-run. Setup
            % and run the simulations 

            start_str = 'Linear: ';
            PB = self.ProgBar(length(self.params.gam_s), 'start_str', start_str);

            result_s = cell(1, 1);
            result_s{1} = self.build_max_sp_local(PB, verbose);

            self.results = result_s;

            if self.savedata
                step_data = self;
                save(step_data.file, 'step_data')
            end
            
        end % END MAIN FUNCTION
        
    end
    
end
