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
        function step_data = build_max_setpoints(self, varargin)
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
