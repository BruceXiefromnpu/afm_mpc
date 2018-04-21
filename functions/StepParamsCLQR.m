classdef StepParamsCLQR
  properties
    sys;
    ref_s;
    du_max;
    Q;
    S;
    gam_s;
    plant;
    N_traj;
    mpc_mode;
    
  end
  
  methods
    function self = StepParamsCLQR(sys, ref_s, du_max, Q, gam_s, ...
        plant, N_traj, mpc_mode, varargin)
      
      if ~strcmp(mpc_mode, 'sparse') && ~strcmp(mpc_mode, 'condensed')
        error('Expected mpc_mode = ("sparse"|"condensed"\n')
      end
      p = inputParser;
      p.addParameter('S', sys.b*0);
      p.parse(varargin{:});
      self.sys = sys;
      self.ref_s = ref_s;
      self.du_max = du_max;
      self.Q = Q;
      self.S = p.Results.S;
      self.gam_s = gam_s;
      self.plant = plant;
      self.N_traj = N_traj;
      self.mpc_mode = mpc_mode;
    end % end contstructor
    
    function [Y, U, dU] =  sim(self, ref, gam, figs)
      %           Qp = dare(self.sys.a, self.sys.b, self.Q, gam);
      %           mpcProb1 = condensedMPCprob(self.sys, 2, self.Q, Qp, 100);

      Qp = dare(self.sys.a, self.sys.b, self.Q, gam, self.S);
      if strcmpi(self.mpc_mode, 'sparse')
          NLQR_prob = sparseMPCprob(self.sys, self.N_traj, self.Q, Qp,...
            gam, self.S);
      elseif strcmpi(self.mpc_mode, 'condensed')
          NLQR_prob = condensedMPCprob_OA(self.sys, self.N_traj, self.Q,...
            Qp, gam, self.S);
      end

      if self.du_max ~= 0 
          CON = CondenCon([], [], NLQR_prob.N_mpc);
          CON.add_input_con('box', self.du_max);
          NLQR_prob.CON = CON;
          %NLQR_prob.add_U_constraint('box', du_max);
      end

      Nx = SSTools.getNxNu(self.sys);
      x0_err = - ref*Nx;
      
      [dU, Xerr] = NLQR_prob.solve(x0_err, 'getX', true);
      U = cumsum(dU);
      X = Xerr + ref*Nx;
      Y = X'*self.sys.c';
      tvec = (0:1:self.N_traj-1)*self.sys.Ts; 

      Y = timeseries(Y(1:end-1), tvec);
      dU = timeseries(dU, tvec);
      
      if nargout == 0
        if ~exist('figs', 'var')
          figs(1) = figure;
          figs(2) = figure();
        end
        if isvalid(figs(1))
          change_current_figure(figs(1));
        else
          figure()
        end
        plot(Y)
        ylabel('y(k')
        if isvalid(figs(2))
          change_current_figure(figs(2));
        else
          figure()
        end
        plot(dU)
        ylabel('$\Delta u(k)$')
      end
    end
    
  end % end methods

end % end classdef
