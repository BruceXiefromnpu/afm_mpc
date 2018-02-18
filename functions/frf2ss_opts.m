classdef frf2ss_opts
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