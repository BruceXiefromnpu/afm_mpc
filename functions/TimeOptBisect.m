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
        max_iter; % Maximum number of times to bump while searching
                 % for the upper bound. Default is 20.
        logger;
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
            self.max_iter = 20;
            self.logger = @fprintf;

        end
        
        % ---------------------------------------------------------------- %
        function [xx, uu, status] = time_opt_bisect(self, x0, xf, varargin)
        % [xx, uu, status] = time_opt_bisect(self, x0, xf, varargin)
        % Solve the minimum-time problem in Discrete Time using a bisection
        % search method.
            p = inputParser;
            p.addParameter('verbose', true);
            p.addParameter('k0', []);
            parse(p, varargin{:});
            
            verbose = p.Results.verbose;
            status = 0;
            % Find upper and lower bounds on N.
            [lowerBound, upperBound, status] = self.find_bounds(x0, ...
                                                              xf, varargin{:});
            if status
                xx = [];
                uu = [];
                status = 1;
                return
            end
            
            % Start in the middle of the bounds. 
            k0 = floor(mean(lowerBound, upperBound));
            % Hueristic to determine max_iter.     
            max_iter_bisect =( upperBound - lowerBound + 1)*2;
            % Perform the Bisection search. 
            for i = 1:max_iter_bisect
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
                % keyboard 
                self.logger(['Final Jval is NaN. Exiting with status ' ...
                              '1']);
                status = 1;
            end
            if verbose
                self.logger('Final CVX Jval = %d\n', Jval);
            end
            

            [~, t_DT, X_DT] = lsim(self.sys, u, (0:1:k0-1)'*self.Ts, x0);

            xx = timeseries(X_DT, t_DT);
            uu = timeseries(u, t_DT);

        end

        % ---------------------------------------------------------------- %
        function [lb, ub, status] = find_bounds(self, x0, xf, varargin)
            
            p = inputParser;
            p.addParameter('verbose', true);
            p.addParameter('k0', []);
            parse(p, varargin{:});
            k0 = p.Results.k0;
             status = 0;
                
            % [lb, ub] = find_bounds(self, x0, xf)
            % Find lower and upper bounds on maximum trajectory length, k0.
            if isempty(k0)
                k0 = length(self.Gam);
            end
            results = self.time_opt_cvx_prob(k0, x0, xf);
            Jval = results.cvx_optval;
            if  Jval < self.TOL && ~isnan(Jval)
                ub = k0;
                [lb, ub,  status] = findBoundsDown(self, x0, xf, ub, self.bump, varargin{:});
            else % Jval > tol, 
                
                lb = k0;
                [lb, ub] = findBoundsUp(self, x0, xf,lb, self.bump, varargin{:});
            end
  
            if isempty(ub)
                str = sprintf(['MyWarning: Exhausted max_iter but no Upper Bound found. ',...
                       'Possible solution is to increase ',...
                        'max_iter=%0.0f'], self.max_iter);
               self.logger( '%s\n', str)
               status = 1;
            end
                            

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



function [lb, ub, status] = findBoundsUp(obj, x0, xf, lb, bump, varargin);
% go up, lb is too small. (the original way)
fprintf('inside Up\n');
p = inputParser;
p.addParameter('verbose', true);
p.addParameter('k0', []);

parse(p, varargin{:});
verbose = p.Results.verbose;
% k0 = p.Results.k0;
status = 0;

upperBound = [];
k0 = lb;
lowerBound = lb;
for i = 1:obj.max_iter
    results = obj.time_opt_cvx_prob(k0, x0, xf);
    Jval = results.cvx_optval;

    if Jval > obj.TOL || isnan(Jval)
        if i>1
            % Jval still too big, so bump the lower bound up
            lowerBound = k0;
        end
        % Add bump instead of multiply. Multiply make the bound
        % too large. This step takes longer, but bisection is
        % quick now.
        k0 = k0 + bump;
    else
        upperBound = k0;
        if verbose
            obj.logger('Bounds found!\n')
            obj.logger('lowerBound = %.0f \n', lowerBound)
            obj.logger('upperBound = %.0f\n', upperBound)            
        end
        break;
    end
end
if isempty(upperBound)
    str = sprintf(['MyWarning: Exhausted max_iter but no Upper Bound found. ',...
           'Have k0 = %d. Possible solution is to increase ',...
            'max_iter=%0.0f'], k0, obj.max_iter);
   obj.logger( '%s\n', str)
   status = 1;
end

lb = lowerBound;
ub = upperBound;
end % find_bounds


function [lb, ub, status] = findBoundsDown(obj, x0, xf,ub, bump, varargin);
p = inputParser;
p.addParameter('verbose', true);
p.addParameter('k0', []);
parse(p, varargin{:});
verbose = p.Results.verbose;
% k0 = p.Results.k0;
status = 0;

fprintf('inside Down\n');
lowerBound = [];
k0 = ub;
upperBound = k0;

for i = 1:obj.max_iter
    results = obj.time_opt_cvx_prob(k0, x0, xf);
    Jval = results.cvx_optval;

    if Jval < obj.TOL || isnan(Jval)
        if i>1
            % Jval still too big, so bump the lower bound up
            upperBound = k0;
        end
        % Add bump instead of multiply. Multiply make the bound
        % too large. This step takes longer, but bisection is
        % quick now.
        k0 = k0 - bump;
    else
        lowerBound = max(0, k0);
        if verbose
            obj.logger('Bounds found!\n')
            obj.logger('lowerBound = %.0f \n', lowerBound)
            obj.logger('upperBound = %.0f\n', upperBound)
        end
        break;
    end
end
if isempty(lowerBound)
    str = sprintf(['MyWarning: Exhausted max_iter but no Lower Bound found. ',...
           'Possible solution is to increase ',...
            'max_iter=%0.0f'], obj.max_iter);
   obj.logger( '%s\n', str)
   status = 1;
end

lb = lowerBound;
ub = upperBound;
end % find_bounds


