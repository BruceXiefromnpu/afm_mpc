classdef TsByMaxRefParams < handle
    properties
        StepData;
        rmax_s;
        file;
        results;
    end
    methods
        function self = TsByMaxRefParams(StepData, rmax_s, file)
            self.StepData = StepData;
            self.rmax_s = rmax_s;
            self.file = file;
        end
        
        function self = run_ts_by_refs_mpc(self, varargin)
        % This method will call self.run_ts_by_refs for each N_mpc
        % in self.StepData.params.N_mpc_s.
            p = inputParser;
            defaultForce = 0;
            p.addParameter('force', defaultForce);
            p.addParameter('fid', 1);
            p.addParameter('verbose', 0);
            parse(p, varargin{:});
            opts.force = p.Results.force;
            opts.fid = p.Results.fid;
            opts.verbose = p.Results.verbose;
            N_mpc_s = self.StepData.params.N_mpc_s;
            results = cell(1, length(N_mpc_s));
            
            for k = 1:length(N_mpc_s)
                % idx_N_mpc = k
                results{k} = run_ts_by_refs_local(self, k, opts);
                
            end
            self.results = results;
            ts_by_rmax_data = self;
            save(self.file, 'ts_by_rmax_data');
        end

        function self = run_ts_by_refs(self, varargin)
            p = inputParser;
            defaultForce = 0;
            p.addParameter('force', defaultForce);
            p.addParameter('fid', 1);
            p.addParameter('verbose', 1);
            parse(p, varargin{:});
            opts.force = p.Results.force;
            opts.fid = p.Results.fid;
            opts.verbose = p.Results.verbose;

            
            result = run_ts_by_refs_local(self, [], opts);
            self.results = result;
            ts_by_rmax_data = self;
            save(self.file, 'ts_by_rmax_data');
        end
        
        function [ax, hands, leg] = plot_ts_v_r2max(self, ax)
            if isempty(self.results)
                error('No Results to plot!')
            end
            
            if isempty(ax)
                ax = gca();
            end
            hold on
            hands = [];
            % keyboard
            for k=1:length(self.results)
                rmax = self.results{k}.ref_s_to_rmax(end);
                gamma = self.results{k}.gamma;
                leg_str = sprintf('$\\gamma = %.0f$, $r_{max} = %.1f$', gamma, rmax);
                h = plot(ax, self.results{k}.ref_s_to_rmax, ...
                         self.results{k}.t_settle_s*1000);
                set(h, 'DisplayName', leg_str);
                xlabel('reference', 'FontSize', 14, 'interpreter', ...
                       'latex')
                ylabel('settle time [ms]', 'FontSize', 14, 'interpreter', ...
                       'latex')
                hands(k) = h;
            end
            leg = legend(hands);
            set(leg, 'interpreter', 'latex', 'FontSize', 14, ...
                     'Location', 'NorthWest');
            
            
        end
        

        function eqq = isequal(self, other, varargin)
            if length(varargin > 0)
                strict = 1;
            else
                strict = 0;
            end
                        
            if strict
                eqq = isequal(self, other);
            else
                eqq1 = isequal(self.StepData.params, ...
                               other.StepData.params);
                eqq2 = isequal(self.rmax_s, other.rmax_s)
                eqq = eqq1 && eqq2;
            end
        end
    end
end


function results = run_ts_by_refs_local(self, idx_N_mpc, opts)
    
    if isempty(idx_N_mpc)
        % Its empty in the linear case. 
        idx_N_mpc = 1;
        N_mpc = [];
    else
        N_mpc = self.StepData.params.N_mpc_s(idx_N_mpc)
    end
    

    
    results = cell(1, length(self.rmax_s));
    for rmax_iter = 1:length(self.rmax_s)
        % Expose for easy access.
        rmax = self.rmax_s(rmax_iter);
        ref_s = self.StepData.params.ref_s;
        gam_s = self.StepData.params.gam_s;
        
        sim_struct = self.StepData.params.sim_struct;
        % keyboard
        max_setpoints_by_gamma = ...
            self.StepData.results{idx_N_mpc}.max_setpoints;
        % For this size reference, find the first gamma we
        % can use. 
        kk = find(max_setpoints_by_gamma >=rmax, 1, 'first');
        gamma = gam_s(kk);
        % if rmax == 4.8
        %     gamma = 9035;
        % end
        
        sim_struct = update_sim_struct(sim_struct, self.StepData.params, gamma, N_mpc);
        % Get the feedback gain that can do it.
        
        % Index to slice ref_s from the first, up the one
        % corresponding to our max setpoints. N.B.
        % max_setpoints and ref_s ARE NOT equal!.
        k_ref = find(ref_s == max_setpoints_by_gamma(kk), 1, ...
                     'first');
        t_settle_s = zeros(1, k_ref);
        % for each ref BELOW rmax, simulate and compute the
        % settling times:
        if opts.verbose > 0
            if isempty(N_mpc) % linear case 
                start_str = sprintf('Lin: rmax = %.2f', rmax);
                PB = ProgressBar(k_ref, 'start_str', ...
                                 start_str, 'fid', opts.fid);
            else
                start_str = sprintf('MPC, N=%.0f, rmax=%.2f', N_mpc, ...
                                    rmax);
                PB = ProgressBar(k_ref, 'start_str', ...
                                 start_str, 'fid', opts.fid);
            end
            PB.upd(0);
        end
        for k=1:k_ref
            mpc_on = 0;
            ref_f = ref_s(k);
            y1 = sim_local(sim_struct, ref_f);
            
            [t_settle, k_s] = settle_time(y1.time, y1.Data, ref_f, 0.01*ref_f,...
                                                   [], [], 30);
            t_settle_s(k) = t_settle;
            if opts.verbose > 0
                PB.upd(k);
            end
        end
        result_rmx_iter.t_settle_s = t_settle_s;
        result_rmx_iter.ref_s_to_rmax = ref_s(1:k);
        result_rmx_iter.gamma = gamma;
        results{rmax_iter} = result_rmx_iter;
    end
    
    results = results;
end

function [Y, U, dU] = sim_local(sim_struct, ref_f)
    options = simset('SrcWorkspace','current');
    % Expose the sime struct to simulink.
    PLANT = sim_struct.PLANT;
    trun = sim_struct.trun;
    mpcProb1 = sim_struct.mpcProb1;
    du_max = sim_struct.du_max;
    mpc_on = sim_struct.mpc_on;
    xss = sim_struct.xss;
    Nx = sim_struct.Nx;
    K_lqr = sim_struct.K_lqr;
    x0 = sim_struct.xss*0;
    uss_0 = 0;
    Ts = PLANT.Ts;
    sim('MPC_fp', [], options)
    
    Y = y_mpcDist;
    U = u_mpcDist;
    dU = du_mpcDist;
end

function sim_struct = update_sim_struct(sim_struct, params, gamma, N_mpc_iter)
   if isempty(N_mpc_iter) % we are in linear mode then.
    sim_struct.K_lqr = dlqr(params.sys.a, params.sys.b,...
                                params.Q, gamma);
   else % if N_MPC_iter not empty --> in MPC mode.
       du_max = sim_struct.du_max;
       N_mpc  = N_mpc_iter;
       Qp = dare(params.sys.a, params.sys.b, params.Q, gamma); 
       sim_struct.mpcProb1 = condensedMPCprob(params.sys, N_mpc,...
                                              params.Q, Qp, gamma);
       sim_struct.mpcProb1.add_U_constraint('box', [-du_max, du_max]);
       sim_struct.K_lqr = params.sys.C*0;
   end
end
