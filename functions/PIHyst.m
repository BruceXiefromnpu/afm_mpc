classdef PIHyst
  %PIHYST Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    r;
    w;
    rp
    wp;
    umax;
    n;
  end
  
  methods 
    function self = PIHyst(umax, n)
      self.r = linspace(0, umax, n)';
      self.w = ones(n, 1);
      self.n = n;
      self.umax = umax;
    end
    
    function invert_hyst_PI2(self)
      [rp, wp] = PIHyst.invert_hyst_PI(self.r, self.w);
    end
  end
  methods (Static)
    function [rp, wp, dp, wsp] = invert_hyst_sat_PI(r, w, d, ws)
      [rp, wp] = PIHyst.invert_hyst_PI(r, w);
      [dp, wsp] = PIHyst.invert_sat(d, ws);
    end
    function [dp, wsp] = invert_sat(d, ws)
      d = d(:);
      ws = ws(:);
      wsp = ws*0;
      dp = d*0;
      
      nd = length(d);
      if mod(nd, 2) ~= 1
        error('Need length(2) odd')
      end
      
      idx_neg = find(d <= 0);
      d_neg = d(idx_neg);
      
      idx_0 = find(d==0);
      d0 = 0;
      
      idx_pos = find(d >= 0);
      d_pos = d(idx_pos);
      
      dp = d*0;

      % negative branch
      for idx=idx_neg'
        dp(idx) = sum(ws(idx:idx_neg(end)).*(d(idx) - d(idx:idx_neg(end))));
      end
      % Positive branch
      for idx=idx_pos'
        dp(idx) = sum(ws(idx_pos(1):idx).*(d(idx) - d(idx_pos(1):idx)));
      end
      
      % now the weights
      wsp = ws*0;
      wsp0 = 1/ws(idx_0);
      wsp(idx_0) = wsp0;

      % negative branch
      for idx=idx_neg(1:end-1)'
        % s1 = sum(ws(1:i));
        % s2 = sum(ws(1:i-1));
        s1 = sum(ws(idx:idx_neg(end)) );
        s2 = sum(ws(idx+1:idx_neg(end)) );
        wsp(idx) = -ws(idx)/(s1*s2);
      end
      if length(idx_pos) > 1
        for idx=idx_pos(2:end)'
          % s1 = sum(ws(1:i));
          % s2 = sum(ws(1:i-1));
          s1 = ws(idx_0) + sum(ws(idx_pos(2):idx) );
          s2 = ws(idx_0) + sum(ws(idx_pos(2):idx-1 ));
          wsp(idx) = -ws(idx)/(s1*s2);
        end
      end
    end
    
    function [rp, wp] = invert_hyst_PI(r, w)
    % [rp, wp] = invert_hyst_PI(r, w)
    %
    % Given a set of PI hysteresis operator parameters r and w, computes the
    % paramters r_prime and w_prime of the inverse operator. 
    
      r = r(:);
      w = w(:);
      rp = r*0;
      wp = w*0;
      
      wp(1) = 1/w(1);
      
      for i=2:length(w)
        s1 = sum(w(1:i));
        s2 = sum(w(1:i-1));
        wp(i) = -w(i)/(s1*s2);
      end
      
      for i=1:length(r)
        rp(i) = sum( w(1:i).*(r(i)-r(1:i)));
      end
    end
    

    function [y] = inverse_hyst_play_sat_op(u, rp, wp, dp,wsp, y0)
      % [y, y_vec_k ] = hyst_play_sat_op(u, r, w, d, ws, y0)
      %
      % Given a control vector u, PI parameters r and w and a hysteresis
      % initial condition y0, computes the output vector y. Also computed are is
      % the internal state sequence, x_vec_k.
      
      n = length(rp);
      wsp = wsp(:);
      wp = wp(:);
      dp = dp(:);
      
     % x_vec_k = zeros(length(u), length(rp));
      
      z = PIHyst.sat_op(u, dp, wsp);
      y = PIHyst.hyst_play_op(z, rp, wp, y0);
    end

    
    function [y, x_mat, dH_dw] = hyst_play_op(u, r, w, y0)
    % [y, y_vec_k ] = hyst_play_op(u, r, w, y0)
    %
    % Given a control vector u, PI parameters r and w and a hysteresis
    % initial condition y0, computes the output vector y. Also computed are is
    % the internal state sequence, x_vec_k.
      
      w = w(:);
      x_mat = zeros(length(u), length(r));
      x_mat(1, :) = y0(:)';
      
      for k=2:length(u)
        uk = u(k);
        for j = 1:length(r)
          x_mat(k, j) = max(uk - r(j), min(uk+r(j), x_mat(k-1, j)));
        end
      end
      
      y = x_mat * w;
      if nargout == 3
        dH_dw = x_mat;
      end
    end
    
    function [y, y_mat, dS_dws, dS_dx] = sat_op(u_vec, d, ws)
    % [y_mat, y] = sat_op(u_vec, d, ws)
    % Implements the (two-sided) saturation operator described in 
    % "Modeling, Identification and Compensation of complex
    % Hysteretic Nonlinearities: A modified Prandtl-Ishlinkii
    % Approach", Klaus Kuhnen
    %
    % y_mat is a matrix of the unweighted values, such that
    % y_mat*ws(:)' = y, where each row corresponds to u_vec(k)
      if length(d) ~= length(ws)
        error('d and ws must of the sime size')
      end
      y = u_vec*0;
      ws = ws(:);
      y_mat = zeros(length(u_vec), length(d));
      for k=1:length(u_vec)
        u_k = u_vec(k);
        Sd_vec = 0*d;
        for i=1:length(d)
          if d(i) == 0
            Sd_vec(i) = u_k;
          elseif d(i) > 0
            Sd_vec(i) = max(u_k - d(i), 0);
          else
            Sd_vec(i) = min(u_k - d(i), 0);
          end
        end
        y_mat(k, :) = Sd_vec;
      end
      
      y = y_mat*ws;
      
      if nargout >= 3
        dS_dws = y_mat;
      end
      if nargout ==4
        dS_dx = y_mat*0;
        for k=1:length(d)
          dS_dx(:,k) = PIHyst.S_prime(u_vec, d(k));
        end
        dS_dx = dS_dx*ws;
      end
      
    end
    
    function sprime = S_prime(x, d)
      % would be faster to not make a copy. Be explicit for now.
      if d > 0
        sprime = x*0+1; 
        sprime(x - d < 0) =0;
      elseif d < 0
        sprime = x*0+1; 
        sprime(x - d > 0) = 0;
      else
        sprime = 0*x +1;
      end
    end
  
    
    function [y, dHS_dw_ws ] = hyst_play_sat_op(u, r, w, d,ws, y0)
      % [y, y_vec_k ] = hyst_play_sat_op(u, r, w, d, ws, y0)
      %
      % Given a control vector u, PI parameters r and w and a hysteresis
      % initial condition y0, computes the output vector y. Also computed are is
      % the internal state sequence, x_vec_k.
      
      n = length(r);
      ws = ws(:);
      w = w(:);
      d = d(:);
      
      [z_k_vec,~, dH_dw] = PIHyst.hyst_play_op(u, r, w, y0);
      [y, ~, dS_dws, dS_dx] = PIHyst.sat_op(z_k_vec, d, ws);
      
      dHS_dw_ws = [dS_dx.*dH_dw, dS_dws];
    end
    
    
    function [r, w, d, ws] = fit_hyst_sat_weights(u, y, Nhyst, Nsat, varargin)
      if length(varargin) >1
        eps_ = varargin{1}
      else
        eps_ = -0.011;
      end
      
      if mod(Nsat,2) ~= 1
        error('Requre Nsat to be odd, but Nsat=%d', Nsat)
      end
      
      n_d = (Nsat - 1)/2;
      
      umax = max(abs(u));
      ymax = max(abs(y)); %
      
      % Create r_s
      r = ([0:Nhyst-1]'./(Nhyst) )*umax;

      % Create d_prime
      id_plus = (1:n_d);
      id_neg = (-n_d:-1);
      dplus = ((id_plus - 0.5)/n_d ) * ymax;
      dmin = ( (id_neg + 0.5)/n_d ) *ymax; 
      dp = [dmin, 0, dplus]';
      
      % To do the fit as a quadratic program, we only the input run
      % through the hyst-op, but unweighted, and the output run
      % through the (inverse) sat-op, but unweighted.
      [~, HU_mat] = PIHyst.hyst_play_op(u, r, r*0, r*0);
      [~, Syp_mat] = PIHyst.sat_op(y, dp, dp*0);
      
      % Create the Hessian of the quadprog
      H = [HU_mat'; -Syp_mat']*[HU_mat, -Syp_mat];
      H = (H+H')/2;
      
      % Create the inequality constraints.
      UH = -eye(Nhyst);
      neg1 = -ones(1, n_d);
      N = 2*n_d+1;
      US = zeros(N, N);
      On = ones(n_d+1);
      Us1 = -triu(On);
      Us2 = -triu(On)';
      US = blkdiag(Us1(1:n_d, 1:n_d), Us2);
      US(1:n_d, n_d+1) = -1;

      Ainq = blkdiag(UH, US);
      binq = [r(:)*0; dp(:)*0] + eps_;
      
      % Create the equality constraint.
      Aeq = [ymax*ones(1, length(r)) - r', zeros(1, length(dp))];
      beq = ymax;      
      
      % Solve the QP
      opts = optimset('quadprog');
      opts.Display = 'off';
      [w_wsp, JVAL] = quadprog(H, H(:,1)*0, Ainq, binq, Aeq, beq, [], [], [], opts);
      %fprintf('Nhyst = %.0f,  Nsat = %.0f, JVAL = %g\n', Nhyst, Nsat, JVAL);
      
      % Split apart the decision variable into w (for hyst) and
      % ws_prime (for sat).
      w = w_wsp(1:length(r));
      wsp = w_wsp(length(r)+1:end);
      
      % We fit d_prime, and ws_prime. Invert to get d and ws.
      [d, ws] = PIHyst.invert_sat(dp, wsp);
    end
    
    function [r, w, d, ws, C, D] = fit_hyst_sat_drift_weights(u, y, Nhyst, Nsat, lams, Ts, varargin)
      if length(varargin) >1
        eps_ = varargin{1}
      else
        eps_ = -0.011;
      end
      
      if mod(Nsat,2) ~= 1
        error('Requre Nsat to be odd, but Nsat=%d', Nsat)
      end
      Ndrift = length(lams);
      A = diag(lams);
      B = lams(:)*0+1;
      C = eye(Ndrift);
      g = ss(A, B, C, zeros(Ndrift,1), Ts);
      t = (0:length(u)-1)'*Ts;
      X = lsim(g, u, t);
      One_ = ones(length(t), 1);
      n_d = (Nsat - 1)/2;
      
      umax = max(abs(u));
      ymax = max(abs(y)); %
      
      % Create r_s
      r = ([0:Nhyst-1]'./(Nhyst) )*umax;

      % Create d_prime
      id_plus = (1:n_d);
      id_neg = (-n_d:-1);
      dplus = ((id_plus - 0.5)/n_d ) * ymax;
      dmin = ( (id_neg + 0.5)/n_d ) *ymax; 
      dp = [dmin, 0, dplus]';
      
      % To do the fit as a quadratic program, we only the input run
      % through the hyst-op, but unweighted, and the output run
      % through the (inverse) sat-op, but unweighted.
      [~, HU_mat] = PIHyst.hyst_play_op(u, r, r*0, r*0);
      [~, Syp_mat] = PIHyst.sat_op(y, dp, dp*0);
      %keyboard
      % Create the Hessian of the quadprog
      
      H = [HU_mat'; -Syp_mat'; X'; One_']*[HU_mat, -Syp_mat, X, One_];
      H = (H+H')/2;
      
      % Create the inequality constraints.
      UH = -eye(Nhyst);
      neg1 = -ones(1, n_d);
      N = 2*n_d+1;
      US = zeros(N, N);
      On = ones(n_d+1);
      Us1 = -triu(On);
      Us2 = -triu(On)';
      US = blkdiag(Us1(1:n_d, 1:n_d), Us2);
      US(1:n_d, n_d+1) = -1;

      Ainq = blkdiag(UH, US, -eye(Ndrift+1));
      %Ainq = [Ainq, zeros(size(Ainq,1), Ndrift+1)];
      binq = [r(:)*0; dp(:)*0; zeros(Ndrift+1, 1)] + eps_;
      
      % Create the equality constraint.
      Aeq = [ymax*ones(1, length(r)) - r', zeros(1, Nsat+Ndrift+1)];
      beq = ymax;      
      
      % Solve the QP
      opts = optimset('quadprog');
      opts.Display = 'off';
      [w_wsp_cd, JVAL] = quadprog(H, H(:,1)*0, Ainq, binq, Aeq, beq, [], [], [], opts);
      %fprintf('Nhyst = %.0f,  Nsat = %.0f, JVAL = %g\n', Nhyst, Nsat, JVAL);
      
      % Split apart the decision variable into w (for hyst) and
      % ws_prime (for sat).
      w = w_wsp_cd(1:length(r));
      wsp = w_wsp_cd(length(r)+1:length(r)+length(dp));
      C = w_wsp_cd(length(r)+length(dp)+1:end-1);
      D = w_wsp_cd(end);
      % We fit d_prime, and ws_prime. Invert to get d and ws.
      [d, ws] = PIHyst.invert_sat(dp, wsp);
    end % fit_hyst_sat_drift_weights
    
    function [r, w] = fit_hyst_weights(u, y, Nhyst, varargin)
      if length(varargin) >1
        eps_ = varargin{1}
      else
        eps_ = -0.011;
      end
      umax = max(abs(u));
      %ymax = max(abs(y)); %
      r = ([0:Nhyst-1]'./(Nhyst) )*umax;
      [~, HU_mat] = PIHyst.hyst_play_op(u, r, r*0, r*0);
      
      H = HU_mat'*HU_mat;
      H = (H+H')/2;
      
      f = -y'*HU_mat; % times 2 but /2
      Ainq = -eye(Nhyst);
      binq = r(:)*0  + eps_;
      [w, JVAL] = quadprog(H, f, Ainq, binq);
      %fprintf('Nhyst = %.0f,  JVAL = %g\n', Nhyst, JVAL+ 0.5*y(:)'*y(:));
      
    end % fit_hyst_weights
    
    
    
    function u = gen_reset_u(t1, t_final, Ts, k1, umax, omega)
    % u = gen_reset_u(t1, t_final, Ts, k1, umax, omega)
    % Generates a control u(k) that is a decaying sinusoid modulated from 0 to
    % t1 by a decaying ramp, and from t1 to t_final by a decaying exponential.
    %
    % Inputs
    % ------
    %  t1, t_final:  double
    %
    %  Ts: sample rate
    %  k1: ramp rate, such that from 0 to t1, a(t) = umax - t*k1
    %  umax: max control amplitude. 
    %  omega: (optional) natural frequency of the sinusoid. Default is 
    %          omega = 1
    %  phi: (optional) phase of the sinusoid. Default is phi = 0;
    %
    % Outputs
    % -------
    %  u : a vector of the control inputs from 0 to t_final
    % 
    % More About
    % ----------
    %   The returned control vector u(k) is supposed to reset the hysteresis
    %   to the relaxed inititial state. This was taken from  
    %    
    %    "Hysteresis and creep modeling and compensation for a piezoelectric
    %    actuator using a fractional-order Maxwell resistive capacitor
    %    approach," Yangfang Liu et al., IOP Smart Materials and Structures,
    %    2013. 
    %
    % 
      n1 = floor(t1/Ts);
      n_final = floor(t_final/Ts);
      
      T1 = (0:n1)'*Ts;
      T2 = (n1+1:n_final)'*Ts;
      T = [T1; T2];
      
      if (umax - k1*t1 < 0)
        delta = 0.1; % arbitrary;
        k1_orig = k1;
        k1 = (umax - delta)/t1;
        warning('umax -k1*t1 <0. Reseting k1 from %f to %f\n', k1_orig, k1)
      end
      
      k2 = k1/(umax - k1*t1);
      if k2 < 0
        error('Need k2 positive, but with the chosen parameters, k2 <0')
      end
      
      a1 = umax - k1*T1;
      a2 = a1(end)*exp(k2*(t1-T2));
      a = [a1; a2];
      if ~exist('omega', 'var')
        omega = 1*2*pi;
      end
      u = a.*sin(omega*T);
    end
    

  end
  
end

