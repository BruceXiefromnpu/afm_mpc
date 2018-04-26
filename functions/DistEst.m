classdef DistEst
% Tools for the creation of disturbance estimation observers and their update
% gains. Most of what is in this class is taken from the paper
%   "Offset-free Reference Tracking for Predictive Controllers,
%    Urban Maeder, and Manfred Morari, CDC 2017"
%
%  -- When working with the deltaUk system, use a state disturbance,
%      DistEst.state_dist_est, which allows us to deal with the pole at z= 1. 
%      The function allows you to pre-design the nominal estimator poles and
%      then add the (usually integrating) pole of the disturbance estimator
%      without changing the original pole locations. 
%  
  
  properties
  end
  
  methods (Static)
    function [L, sys_obs_dist, Ident_obs, eNs_1] = state_dist_est(sys,...
                                                  Lx, p_int, Bd, Ad)
    % [L, sys_obs_dist] = state_dist_est(sys, Lx, p_int, Bd, Ad)
    %   
    % When working with the deltaUk system, we have a pole at z=1, which
    % prevents us from modeling the disturbance as an output disturbance (as
    % far as I can tell. 
    % 
    % Instead, use this function, which models things as a state disturbance.
    %
    % Inputs
    % ------
    % sys   : the state space system under consideration. Pole at z=1 allowed.
    % Lx    : the nominal observer gain for sys.
    % p_int : the desired closed-loop pole location of the disturbance state.
    % Bd    : (optional) Disturbance input matrix. By default, Bd = sys.b. 
    % Ad    : (optional) Disturbance state transition matrix. By default, Ad = I.
    %
    % Outputs:
    % -------
    % L            : the full observer gain which should yield
    %                  \sigma(A_ - LC_) = \sigma(A - LxC) \union p_int
    % sys_obs_dist : the state space observer system with augmented state xd
    %                such that 
    %                A_ = [A, Bd
    %                      0, Ad]
    %                B_ = [B;
    %                      0]
    %                C_ = [C, 0]
    %  I_obs:        I_obs*x_extended will return the original part of the
    %                state. IE, I_obs is the first Ns rows of eye(Ns+ndist)
    %         
    %  eNs_1:        The last row of eye(Ns+ndist), ie the disturbance state
    %                xd can be obtained via xd = eNs_1*x_aug.    
    
      ns    = size(sys.b, 1);
      ndist = size(sys.c, 1);
      
      if ~exist('Ad', 'var') || ~exist('Bd', 'var')
        Ad = eye(ndist);
        Bd = zeros(ns, ndist);
        Bd(end) = 1;
        Cd = zeros(ndist, ndist);
      else
        ndist = size(Ad, 1);
      end      
      
      I_ns = eye(ns);
      I_nd = eye(ndist);
      Z = zeros(ns, 1);
      [sys_obs_dist, Ident_obs, eNs_1]  = distEstObserver(sys, Ad, Bd, Cd);
      
      % The transformations require Kx, but they seem to work fine without it.
      % (ie, with Kx = 0)I need to investigate this. So just make one for now. 
      Kx = dlqr(sys.a, sys.b, I_ns, 1);
      
      H_bar = I_nd + sys.c*((I_ns - sys.a + sys.b*Kx)\Lx)
      
      T = [I_ns, -(I_ns - sys.A + Lx*sys.c)\Bd;
        zeros(1, ns),  1];
      
      %A_bar = sys_dist_obs.a - [Lx; 0]*sys_dist_obs.c
      %At_bar = T*(A_bar/T);  % T*A*T^-1
      %Ct_bar = H_bar*[sys.c, sys.c*( (I_ns-sys.a + Lx*sys.c)\Bd)];
      
      c_ = H_bar*sys.c*( (I_ns - sys.A +Lx*sys.C)\Bd)
      
      Ld_bar = place(I_nd', c_', p_int);
      
      L =[Lx; 0] +  T\[Z;
                      Ld_bar*H_bar];
      
    end
     function [sys_obs, I_obs, eNs_1] = distEstObserver(PLANT, Ad, Bd, Cd)
     % [sys_obs, I_obs, eNs_1] = distEstObserver(PLANT, Ad, Bd, Cd)
     %
     % Adds a constant disturbance estimate to the system, PLANT, which is
     % discrete time.
     %
     % Inputs:
     %   PLANT: discrete time ss system
     %   Ad, Bd, Cd. Disturbances matrices. If these are neglected, the function
     %   assumes that the number of disturbances to model is equal to the number
     %   of system outputs, ie, ndist = size(PLANT.c, 1). It then sets
     %  Ad = eye(ndist), Bd = zeros(ndist,1), Cd = eye(ndist). This is the
     %  standard setup for output disturbance estimation I have been using. (as
     %  of jan 2017).
     %
     % Outputs:
     % sys_obs(PHI, Gam, H) where
     %  PHI = [A, Bd
     %         0, Ad]
     %  Gam = [B;
     %         0]
     %  H   = [C, Cd]
     %
     %  I_obs: I_obs*x_extended will return the original part of the state. IE,
     %         I_obs is the first Ns rows of eye(Ns+ndist)
     %  eNs_1: The last row of eye(Ns+ndist), ie the disturbance state xd can be
     %         obtained via xd = eNs_1*x_aug.
      if PLANT.Ts ==0
        error('Function only defined for discrete time plants')
      end
      
      %
      Ns    = size(PLANT.b, 1);
      ndist = size(PLANT.c, 1);
      
      if ~exist('Ad', 'var') || ~exist('Bd', 'var') || ~exist('Cd', 'var')
        Ad = eye(ndist);
        Bd = zeros(Ns, ndist);
        Cd = eye(ndist);
      else
        ndist = size(Ad, 1);
      end
      
      [A, B, C, D] = ssdata(PLANT);
      A_obs = [A,             Bd;
        zeros(ndist,Ns),  Ad];
      B_obs   = [B; zeros(ndist,ndist)];
      C_obs   = [C, Cd];
      sys_obs = ss(A_obs, B_obs, C_obs, 0, PLANT.Ts);
      
      Ns_obs  = Ns + ndist;
      II      = eye(Ns_obs);
      
      
      eNs_1 = II(end-ndist+1:end,:);
      I_obs = II(1:Ns,:);
      
    end
  end
  
end

