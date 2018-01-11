% [x, u, S] = TimeOptBisect(sys, umax)
classdef TimeOptBisect
    properties
        sys; % The discrete time LTI ss to solve the problem for.
        Ts;  % Sampling period.
        umax; % Upper control bound.
        umin; % Lower control bound.
        PHI;  % PHI = sys.a
        Gam;  % gam = sys.b
        bump; % How much to bump the upper bound by when searching. 
              % Default is 10.
        TOL;  % Tolerance on ||xf - X(N)*|| for deciding if we are
              % close enough.
        maxIter; % Maximum number of times to bump while searching
                 % for the upper bound. Default is 20.
    end
    methods
        function self = TimeOptBisect(sys, umax)
        % self = TimeOptBisect(sys, umax)
        % Main constructor method. 
            self.umax = umax;
            self.umin = -umax;
            self.sys = sys;
            self.Ts = sys.Ts;
            self.PHI = sys.a;
            self.Gam = sys.b;

            self.bump = 10;
            self.TOL = 5e-3;
            self.maxIter = 20;

        end
        
        % ---------------------------------------------------------------- %
        function [xx, uu, S] = time_opt_bisect(self, x0, xf, varargin)
        % Solve the minimum-time problem in Discrete Time using a bisection
        % search method.
            p = inputParser;
            p.addParameter('verbose', true)
            parse(p, varargin{:});
            verbose = p.Results.verbose;
            
            % Find upper and lower bounds on N.
            [lowerBound, upperBound] = self.find_bounds(x0, xf, 'verbose', verbose);
            % Start in the middle of the bounds. 
            k0 = floor(mean(lowerBound, upperBound));
            % Hueristic to determine max_iter.     
            maxIter_bisect =( upperBound - lowerBound + 1)*2;
            % Perform the Bisection search. 
            for i = 1:maxIter_bisect
                results = self.time_opt_cvx_prob(k0, x0, xf);
                Jval = results.cvx_optval;
                if (Jval >= self.TOL)
                    lowerBound = k0;
                    k0 = ceil( mean([upperBound, lowerBound]));
                elseif(Jval < self.TOL)
                    upperBound = k0;
                    k0 = ceil(mean([upperBound, lowerBound]));
                else
                    lowerBound = k0+1;
                    k0 = k0+1;
                end
                if( abs(upperBound - lowerBound) < 2)
                    k0 = upperBound;
                    break
                end
            end

            results = self.time_opt_cvx_prob(k0, x0, xf);
            u = results.U;
            Jval = results.cvx_optval;
            if isnan(Jval)
                keyboard
            end
            if verbose
                fprintf('Final CVX Jval = %d\n', Jval);
            end
            

            [~, t_DT, X_DT] = lsim(self.sys, u, (0:1:k0-1)'*self.Ts, x0);

            xx = timeseries(X_DT, t_DT);
            uu = timeseries(u, t_DT);
            S  = sprintf('Settle Time (discrete, all states): %.5f', t_DT(end));

        end

        % ---------------------------------------------------------------- %
        function [lb, ub] = find_bounds(self, x0, xf, varargin)
            p = inputParser;
            p.addParameter('verbose', true)
            parse(p, varargin{:});
            verbose = p.Results.verbose;
                % if length(varargin) == 1
                    
        % [lb, ub] = find_bounds(self, x0, xf)
        % Find lower and upper bounds on maximum trajectory length, k0.
            lowerBound = length(self.Gam);
            k0 = lowerBound;
            for i = 1:self.maxIter
                results = self.time_opt_cvx_prob(k0, x0, xf);
                Jval = results.cvx_optval;

                if Jval > self.TOL || isnan(Jval)
                    if i>1
                        % Jval still too big, so bump the lower bound up
                        lowerBound = k0;
                    end
                    % Add bump instead of multiply. Multiply make the bound
                    % too large. This step takes longer, but bisection is
                    % quick now.
                    k0 = k0 + self.bump;
                else
                    upperBound = k0;
                    if verbose
                        fprintf('Bounds found!\n')
                        fprintf('upperBound = %.0f\n', upperBound)
                        fprintf('lowerBound = %.0f \n', lowerBound)
                    end
                    break;
                end
            end
            if isempty(upperBound)
               error(['Exhausted maxIter but no Upper Bound found. ',...
                       'Possible solution is to increase maxIter']);
            end
                            
            lb = lowerBound;
            ub = upperBound;
        end % find_bounds
        
        % ---------------------------------------------------------------- %
        function results = time_opt_cvx_prob(self, k0, x0, xf)
        % results = time_opt_cvx_prob(self, k0, x0, xf)
        % Solve the minimization problem
        % min  ||X(N) - xf||
        % s.t. 
        %     |U|_{inf} < umax
        %     X(N) = dynamics
        % 
        % Inputs
        % ------
        %   k0 : control horizon to use. 
        %   x0 : Initial condition
        %   xf : Desired final state
        %    
        % Outputs
        % -------
        %   results : a cvx structure. Contains 
        %   results.U : the optimal control vector
        %   results.cvx_optval : value of optimum.
              
        
            % path for cvx
            addpath(genpath(fullfile(getMatPath, 'solvers/cvx')))
            C=[];
            for k=0:k0-1
                C = [self.PHI^k*self.Gam C];
            end
            e_k0 = zeros(1,k0);
            e_k0(end) = 1; % Kth unit vector. [0 0 .... 0 0 1]. To pick out last element of u.
            I = eye(length(self.Gam));

            cvx_begin quiet
            variable u(k0);
            variable t;
            L = self.PHI^k0*x0 + C*u;
            % minimize norm(L - xf)
            % minimize norm(u)
            minimize norm(t)
            subject to
            norm(u, Inf) <= self.umax;
            % Enforce steady state requirement
            % xss = Axss + Buss --> (I-A)xss = Buss
            % (I-A)xss - Buss = 0. Not sure this is necessary.
            % If minimization succeeds, should be the case anyway.
            (I-self.PHI)*L - self.Gam*(e_k0*u) ==0;
            L - xf == 0;
            cvx_end

            results.U = u;
            results.cvx_optval = cvx_optval;

        end
    end % methods
end %classdef
