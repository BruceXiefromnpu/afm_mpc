classdef SimAFM

  properties
    PLANT;
    controller;
    Nx;
    Nbar;
    useNbar;
    sys_obs;
    L;
    du_max;
    x0_obs;
    x0;
    sys_obs_fp;
    
    thenoise;
    step_amp; % disturbance
    r;
    w;
    rp;
    wp;
    d;
    ws;
    dp;
    wsp;
    gdrift_inv;
    gdrift;
   
    isfxp;
    nw;
    nf;
    
  end

  methods
    
    function self = SimAFM(PLANT, controller, Nx, sys_obs, L, du_max, isfxp, varargin)
    % SimAFM(PLANT, controller, Nx, sys_obs, L, du_max, isfxp, varargin)  
    % varargin is name-value pairs of 
    % 'thenoise', 'step_amp', 'r', 'w', 'rp', 'wp', 'gdrift',
    % 'gdrift_inv', 'nw', 'nf'.
    
      fi_zero = fi(0, 1, 16, 11);
      kd_sys =  tf(1,1, PLANT.Ts);
      
      p = inputParser();
      p.addParameter('thenoise', []);
      p.addParameter('step_amp', 0);
      p.addParameter('r', 0);
      p.addParameter('w', 0);
      p.addParameter('rp', 0);
      p.addParameter('wp', 0);
      p.addParameter('d', 0);
      p.addParameter('ws', 0);
      p.addParameter('dp', fi_zero);
      p.addParameter('wsp', fi_zero);
      p.addParameter('gdrift',  kd_sys);
      p.addParameter('gdrift_inv', kd_sys);
      p.addParameter('nw', 0);
      p.addParameter('nf', 0);
      p.addParameter('useNbar', false)
      p.parse(varargin{:});
      
      self.thenoise = p.Results.thenoise;
      self.step_amp =  p.Results.step_amp;
      self.r = p.Results.r;
      self.w = p.Results.w;
      self.rp = p.Results.rp;
      self.wp = p.Results.wp;
      self.d = p.Results.d;
      self.ws = p.Results.ws;
      self.dp = p.Results.dp;
      self.wsp = p.Results.wsp;
      self.useNbar = p.Results.useNbar;
      
      
      self.PLANT = PLANT;
      self.controller = controller;
      self.Nx = Nx;
      if self.useNbar
        if ~isnumeric(controller)
          error('useNbar only makes sense if controller is feedback gain, K')
        end
        self.Nbar = controller*Nx;
      else
        self.Nbar = [];
      end
      self.sys_obs = sys_obs;
      self.L = L;
      self.du_max = du_max;
      self.isfxp = isfxp;
      self.x0_obs = sys_obs.b*0;
      self.x0 = PLANT.b*0;
      
      
      self.gdrift = p.Results.gdrift;
      self.gdrift_inv = p.Results.gdrift_inv;
      self.nw = p.Results.nw;
      self.nf = p.Results.nf;
    end
    
    function [Y, U_full, U_nominal, dU, Xhat, Xerr] = sim(sim_obj, ref_traj, dist_traj)
    % NOTE: rather than self, call it sim obj so it is more
    % readable in simulink.  
    % [Y, U, dU] = sim_MPC_fp(self, ref_f)
      %
      % Inputs
      % ------
      %  self : a structure which must have the following feilds
      %
      %    self.K_lqr;
      %    self.PLANT;
      %    self.trun;
      %    self.mpcProb1;
      %    self.du_max;
      %    self.mpc_on;
      %    self.Nx;
      %
      %  ref_f : (scalar) the setpoint to track
      % Outputs
      % -------
      %  Y : (timeseries) Plant Output
      %  U : (timeseries) Accumulated deltaU
      %  dU : (timeseries) deltaU(k) Control (ie, output of MPC block)
      %
      % ------------------------------------------------------------------- %
      % Pull out all the data stored in self to expose it to
      % simulink. There must be a better way...
      VSS_LINEAR_CONTROL=Simulink.Variant('VSS_CONTROL_MODE==1');
      VSS_MPC_CONTROL=Simulink.Variant('VSS_CONTROL_MODE == 2');
      
 
      %#ok<*NASGU>
      VSS_LINEAR_CONTROL=Simulink.Variant('VSS_CONTROL_MODE==1');
      VSS_MPC_CONTROL=Simulink.Variant('VSS_CONTROL_MODE == 2');

      fprintf('-----------------------\n');  
      % Expose the sim struct to simulink.
      if ~isnumeric(sim_obj.controller) % fi is numeric too.
        mpcProb = sim_obj.controller;
        VSS_CONTROL_MODE=2;
        % keyboard
        mpcProb.reset_warm_start_data();
        fprintf('simulating as MPC\n')
      else
        K_lqr = sim_obj.controller;
        VSS_CONTROL_MODE=1;
        fprintf('simulating as LINEAR\n');
      end
      if ~exist('dist_traj', 'var')
        dist_traj = ref_traj;
        dist_traj.Data = dist_traj.Data*0;
      end
      PLANT = sim_obj.PLANT;
      trun = ref_traj.Time(end);
      
      x0 = SSTools.getNxNu(PLANT)*0;
      uss_0 = 0;
      
      Ts = PLANT.Ts;
      
      % Expose the sim struct to simulink.
      PLANT = sim_obj.PLANT;
      trun = ref_traj.Time(end);
      Ts = PLANT.Ts;
      x0 = sim_obj.x0;
      
      
      thenoise = sim_obj.thenoise;
      %thenoise = timeseries(ref_traj.Time*0, ref_traj.Time);
      
      sim_obj.step_amp = 0;

      r = sim_obj.r;
      w = sim_obj.w;
      d = sim_obj.d;
      ws = sim_obj.ws;
      
      rp = sim_obj.rp;
      wp = sim_obj.wp;
      dp = sim_obj.dp;
      wsp = sim_obj.wsp;
      
      % sim_obj.gdrift_inv = tf(1,1, PLANT.Ts);
      % sim_obj.gdrift = tf(1,1, PLANT.Ts);
      
      [ ndist, Ns_obs] = size(sim_obj.sys_obs.c);
      Ident_obs = eye(Ns_obs);
      
      % the last row
      C_ydist = Ident_obs(end-ndist+1:end, :);
      % all rows but the last 1
      Ident_obs = Ident_obs(1:end-ndist, :);
      
      
      % ------------------------  RUN THE SIM ---------------------------- %
      options = simset('SrcWorkspace','current');
      if ~sim_obj.isfxp
        fprintf('simulating as floating-point\n');
        sim('AFMss_fp_obshas_uk', [], options)
      else
        % Turn off warning about fxp states not getting logged.
        fprintf('simulating as FIXED-point\n');
        nw = sim_obj.nw;
        nf = sim_obj.nf;
        id = 'Simulink:Engine:BlkIgnoringUsedAsDStateFlag';
        warning('off', id);
        

        if ~isnumeric(sim_obj.controller)
        % Since we are in the fixed-point case, check if we are doing MPC. It
        % is too slow to use the s-fun, we have to go back to the matlab fun,
        % this sucks but is what it is.
          IH_fxp = sim_obj.controller.I_HL;
          ML_x0_fxp = sim_obj.controller.ML;
          beta_fxp =  sim_obj.controller.beta;
          maxIter =  sim_obj.controller.maxIter;
          N_mpc =  sim_obj.controller.N_mpc;
          z0 =  sim_obj.controller.warm_start_data(1:N_mpc);
          y0 =  sim_obj.controller.warm_start_data(N_mpc+1:end);
          du_max = sim_obj.controller.uMax;
          nw_fgm = sim_obj.controller.nw;
          nf_fgm = sim_obj.controller.nf;
          Ns = size(ML_x0_fxp,2);
        end
        %keyboard
        sim('FXP_AFMss_obshas_uk', [], options)
      end
      % provides Y, U, dU, Xhat, U_nominal      
      
      % hyst_sat_data = struct('u_gdinv', u_gdinv, 'u_sinv', u_sinv,...
      %   'u_hinv', u_hinv, 'dp', dp, 'wsp', wsp, 'rp', rp, 'wp', wp)

      % provides Y, U, dU
    end
   
    function write_control_data(self, data_path, ref, varargin)
    % write_control_data(self, data_path, ref, varargin)
      
      if ~isnumeric(ref)
        ref_traj_path = varargin{1};
        ref_f_1 = 0;
        fid = fopen(ref_traj_path, 'w+');
        fprintf(fid, '%.12f, ', ref.Data(:)'); % as a row.
        fprintf(fid, '\n');
        fclose(fid);
      else
        ref_f_1 = ref(1);
      end

      Ns = size(self.sys_obs.b,1);
      Ns_mpc = Ns; %*size(self.sys_obs.b,2);
      umax = 0;
      
      fid = fopen(data_path, 'w+');
      fprintf(fid, '%d, %.5f, 0, %f, %.12f,  %.12f\n', Ns, ref_f_1, umax, self.du_max, self.Nbar);
      
      
      if isnumeric(self.controller)
        MPC_mat = [];
        fprintf(fid, '0, 0, 0, 0\n');
        K = self.controller;
      else
        % Should I make the FGMProb class do this??
        Nmpc = self.controller.N_mpc;
        I_HL = self.controller.I_HL;
        ML = self.controller.ML;
        MPC_mat = [double(I_HL), zeros(Nmpc,2), double(ML), zeros(Nmpc, 2)];
        K = zeros(size(ML,2),1);
        
        fprintf(fid, '%.12f, %f, %f, %f\n', double(self.controller.beta),...
                self.controller.maxIter, self.controller.N_mpc, Ns_mpc);
      end
      
      if isempty(self.sys_obs_fp)
        error(['to write the controller data, set property ' ...
               'sys_obs_fp'])
      end
      
      if self.useNbar
        AllMatrix = packMatrixDistEst(self.sys_obs_fp,...
          double(self.L), double(K), []);
      else
        AllMatrix = packMatrixDistEst(self.sys_obs_fp,...
          double(self.L), double(K), double(self.Nx));
      end

      % These get written all as one column
      for k=1:size(AllMatrix,1)-1
        fprintf(fid, '%.12f, ', AllMatrix(k,:));
      end
      fprintf(fid, '%.12f\n', AllMatrix(end,:));
      %fprintf('\n');
      
      if ~isempty(MPC_mat)
        for k=1:size(MPC_mat, 1)-1
          fprintf(fid, '%.12f, ', MPC_mat(k,:));
        end
        fprintf(fid, '%.12f\n', MPC_mat(end,:));
        fclose(fid);
      end
    end % write control data
  
   end %methods
end

% allmatrix = packMatrix(sys_obs, L, K, Nx)
%
% Packs observer system matrices, control gain, and observer gain into a
% single column vector. To be read by labview for experiement.
% This is performs the same functionality but in addition to packing 
% [B;L;K';Arows;]
% it appends xss to the end s.t
% [B;L;K';Arows; xss]. To use this with MPC, provde K as empty or zero.
%
% ASSUMES A IS PROVIDED AS A~ = A-LC
function AllMatrix = packMatrixDistEst(sys_obs, L,K, Nx)
    A = sys_obs.a;
    B = sys_obs.b;
    % C = sys_obs.c;
    if isempty(K)
      K = B*0;
    else
      K_T = K';
    end
    
    Ns = size(B, 1);

    AllMatrix = [B(:); 
                 L(:);
                 K_T(:)]; % filler for K
    
    for i=1:Ns
       AllMatrix = [AllMatrix; A(i,:)']; 
    end
    
    AllMatrix = [AllMatrix;
                 Nx];

end