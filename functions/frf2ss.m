classdef frf2ss

    properties
       frf;
       freqs_rad;
       frf_even
       freqs_rad_even;
       impulse_experimental;
       Nd; 
       opts;
       
       ss_full;        
    end
    
    methods
        function self = frf2ss(frf, freqs_rad, Nd, opts)
            self.opts = opts;
            
            [frf, freqs_rad] = monotonicFRF(frf, freqs_rad);
            
            self.freqs_rad = freqs_rad;
            self.frf = frf;
            self.Nd = Nd;
            
             t = [0:opts.Ts:opts.impulse_length_sec]';
             N = length(t);
             K = length(freqs_rad);

             % For the fourier transform, N we need evenly spaced frequencies, where N
             % is length of impulse response.

             freqs_even = [0:(1/(N*opts.Ts)):floor((1/opts.Ts)*0.5)]'; %Frequencies up to nyquist.
             ws_even    = freqs_even*2*pi;
             self.freqs_rad_even = ws_even;
             
            % Interpolate the FRF we have onto the evenly spaced frequencies.
            K_freq_max = find(ws_even <= freqs_rad(end), 1, 'last')
            K_freq_min = find(ws_even <= freqs_rad(1), 1, 'last')

            self.frf_even = interp1(freqs_rad, frf, ws_even(1:K_freq_max), 'spline');

            % Remove delay from he experimental FRF
            if Nd
                gdelay = (exp(j*ws_even*opts.Ts)).^Nd;
                self.frf_even = self.frf_even.*gdelay(1:K_freq_max);
            end

            % Window the FRF for frequencies past what we have, because
            % interpolation becomes inaccurate.
            self.frf_even(K_freq_max:length(ws_even)) = 0;

            % Get the impulse response
            self.impulse_experimental = frf2impulse(self.frf_even);
            
            % Do the ERA on the impulse response
            [Ad, Bd, Cd] = impulse2ss(reshape(self.impulse_experimental,1,1, []), opts.r, opts.s);
            self.ss_full.A = Ad;
            self.ss_full.B = Bd;
            self.ss_full.C = Cd;
        end
        
        function sys = realize(self, Ns)
            A = self.ss_full.A(1:Ns,1:Ns);
            B = self.ss_full.B(1:Ns);
            C = self.ss_full.C(1:Ns);
            D = 0;
            sys = ss(A, B, C, D, self.opts.Ts, 'InputDelay', round(self.Nd));
            
            
            
            
        end
        
    end
end













