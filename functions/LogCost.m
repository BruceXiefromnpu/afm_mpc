classdef LogCost <handle
    %LOGCOST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        theta;
        pls;
        zrs;
        k;
        Ts;
        nz;
        np;
        G_exp;
        omegas;
    end
    
    methods
        function self = LogCost(Ts, G_exp_frf, omegas)
            self.theta = [];
            self.Ts = Ts;
            self.G_exp = G_exp_frf;
            self.omegas = omegas;
        end
        function add_complex(self, z_or_p, zeta, wn)
           
            s = -zeta*wn + 1j*wn*sqrt(1-zeta^2);
            
            if z_or_p == 'z'
                self.zrs = [self.zrs, exp(s*self.Ts), exp(conj(s)*self.Ts)];
            elseif z_or_p == 'p'
                self.pls = [self.pls, exp(s*self.Ts), exp(conj(s)*self.Ts)];
            else
                error('Expected z_or_p = ("z"|"p"), but recieved: %s\n', z_or_p);
            end
        end
        
        function add_real(self, z_or_p, w)
            s = -w;
            
            if z_or_p == 'z'
                self.zrs = [self.zrs, exp(s*self.Ts)];
            elseif z_or_p == 'p'
                self.pls = [self.pls, exp(s*self.Ts)];
            else
                error('Expected z_or_p = ("z"|"p"), but recieved: %s\n', z_or_p);
            end
        end
        
        function update_theta(self)
            if isempty(self.k)
                k = 1;
            else
                k = self.k;
            end
            g = zpk(self.zrs, self.pls, k, self.Ts);
            [num, den] = tfdata(g);
            num = num{1};
            den = den{1};
            self.nz = length(self.zrs);
            self.np = length(self.pls);
            
            self.theta = [num(end-self.nz:end), den(end-self.np+1:end)];
        end
        function solve_lsq(self)
            self.theta = lsqnonlin(@self.logcost, self.theta);
        end
        function G = realize(self)
           G = tf(self.theta(1:self.nz+1), [1 self.theta(self.nz+2:end)], self.Ts); 
            
        end
        function J_vec = logcost(self, theta)
            % We can weight the cost function by frequency spacing.
            % w_weight = [self.omega(1);self.omega;self.omega(end)];
            % logW = log(w_weight(3:end)) - log(w_weight(1:end-2));
            
            H_frf = H(self.omegas, theta, self.nz, self.Ts);
            J_vec = abs(log(H_frf./self.G_exp));
        end
        
    end
    
end




function H_frf = H(w, theta, nz, Ts)

num = theta(1:nz+1);
den = theta(nz+2:end);
% keyboard
g = tf(num, [1 den], Ts);
H_frf = squeeze(freqresp(g, w));

end
