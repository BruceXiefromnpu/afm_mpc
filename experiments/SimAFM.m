classdef SimAFM

  properties
    PLANT;
    controller;
    Nx;
    sys_obs;
    L;
    du_max;
    x0_obs;
    x0;
    
    thenoise;
    step_amp; % disturbance
    r;
    w;
    rp;
    wp;
    gdrift_inv;
    gdrift;
   
    isfxp;
    nw;
    nf;
    
  end

  methods
    
    function self = SimAFM(PLANT, controller, Nx, sys_obs, L, du_max, isfxp, varargin)
      self.PLANT = PLANT;
      self.controller = controller;
      self.Nx = Nx;
      self.sys_obs = sys_obs;
      self.L = L;
      self.du_max = du_max;
      self.isfxp = isfxp;
      self.x0_obs = sys_obs.b*0;
      self.x0 = PLANT.b*0;
      kd_sys =  tf(1,1, PLANT.Ts);
      
      p = inputParser();
      p.addParameter('thenoise', []);
      p.addParameter('step_amp', 0);
      p.addParameter('r', 0);
      p.addParameter('w', 0);
      p.addParameter('rp', 0);
      p.addParameter('wp', 0);
      p.addParameter('gdrift',  kd_sys);
      p.addParameter('gdrift_inv', kd_sys);
      p.addParameter('nw', 0);
      p.addParameter('nf', 0);
      p.parse(varargin{:});
      
      self.thenoise = p.Results.thenoise;
      self.step_amp =  p.Results.step_amp;
      self.r = p.Results.r;
      self.w = p.Results.w;
      self.rp = p.Results.rp;
      self.wp = p.Results.wp;
      self.gdrift = p.Results.gdrift;
      self.gdrift_inv = p.Results.gdrift_inv;
      self.nw = p.Results.nw;
      self.nf = p.Results.nf;
    end
    
    function [Y, U, U_nominal, dU, Xhat, Xerr] = sim(sim_obj, ref_traj)
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
      
      % VSS_NO_HYST_NO_DRIFT=Simulink.Variant('VSS_EXT_DIST_MODE==1');
      % VSS_HYST_AND_DRIFT=Simulink.Variant('VSS_EXT_DIST_MODE==2');
      % VSS_HYST_ONLY=Simulink.Variant('VSS_EXT_DIST_MODE==3');
      % VSS_EXT_DIST_MODE = 1;
      
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
      
      PLANT = sim_obj.PLANT;
      trun = ref_traj.Time(end);
      
      x0 = SSTools.getNxNu(PLANT)*0;
      uss_0 = 0;
      
      Ts = PLANT.Ts;
      
      % Expose the sim struct to simulink.
      PLANT = sim_obj.PLANT;
      trun = ref_traj.Time(end);
      Ts = PLANT.Ts;
      %trun = 5*Ts;
      x0 = sim_obj.x0;
      % uss_0 = 0;
      
      thenoise = sim_obj.thenoise; 
      thenoise = timeseries(ref_traj.Time*0, ref_traj.Time);
      
      sim_obj.step_amp = 0;

      r = sim_obj.r;
      w = sim_obj.w;
      r = 0;
      w = 0;
      
      rp = sim_obj.rp;
      wp = sim_obj.wp;
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
        sim('FXP_AFMss_obshas_uk', [], options)
      end
      % provides Y, U, dU, Xhat, U_nominal      
      
      % provides Y, U, dU
    end
   
    function write_control_data(self, data_path, ref, varargin)
      if ~isnumeric(ref)
        ref_traj_path = varargin{1};
        ref_f_1 = 0;
      else
        ref_f_1 = ref(1);
      end

      Ns = size(self.sys_obs.b,1) - size(self.sys_obs.c,2);
      Ns_mpc = Ns + size(self.sys_obs.b,2);
      umax = 0;
      
      fid = fopen(data_path, 'w+');
      fprintf(fid, '%d, %.5f, 0, %f, %.12f\n', Ns, ref_f_1, umax, self.du_max);
      
      
      if isnumeric(self.controller)
        MPC_MAT = [];
        fprintf(fid, '0, 0, 0, 0\n');
      else
        % Should I make the FGMProb class do this??
        Nmpc = self.controller.N_mpc;
        I_HL = self.controller.I_HL;
        ML = self.controller.ML;
        MPC_MAT = [double(I_HL), zeros(Nmpc,2), double(ML), zeros(Nmpc, 2)];
        K = zeros(size(ML,2),1);
        
        fprintf(fid, '%.12f, %d, %d, %d\n', double(self.controller.beta),...
                self.controller.maxIter, self.controller.N_mpc, Ns_mpc);
      end
      
      % These get written all as one column
      for k=1:size(AllMatrix,1)
        fprintf(fid, '%.12f, ', AllMatrix(k,:));
      end

      for k=1:size(MPCMatrix, 1)
        fprintf(fid, '%.12f, ', MPCMatrix(k,:));  
      end
      fprintf(fid, '\n');
      fclose(fid);
      
      
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