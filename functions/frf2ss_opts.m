classdef frf2ss_opts
  %   opts = frf2ss_opts(name1, value1,...)
  % Builds an options structure for for frf2ss via name-value pairs. Options
  % are:
  % 'Ts', sampling_time
  % 'nd', integer number of delay. This does not appear to get used in frf2ss. 
  %       This should be fixed.
  % 'impulse_length_sec', default = 0.1
  % 'r', integer
  % 's', integer
  %
  % r and s and the width and hieght of the hankel matrix. Defaults to 500x500.
  % 
   properties
       Ts;
       nd;
       impulse_length_sec;
       r;
       s;
       
   end
    
   methods
       function self = frf2ss_opts(varargin)
          p = inputParser;
          p.addParameter('Ts', 40e-6);
          p.addParameter('nd', 0);
          p.addParameter('impulse_length_sec', 0.1);
          p.addParameter('r', 500);
          p.addParameter('s', 500);
          
          parse(p, varargin{:});
          
          self.Ts = p.Results.Ts;
          self.nd = p.Results.nd;
          self.impulse_length_sec = p.Results.impulse_length_sec;
          self.r = p.Results.r;
          self.s = p.Results.s;
          
       end
   end
    
    
end