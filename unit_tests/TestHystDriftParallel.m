classdef TestHystDriftParallel < matlab.unittest.TestCase
  
  
 properties
   r;
   w;
   d
   ws
   gdrift;
 end
 methods
   function self = TestHystDriftParallel()
     self.r  = [0    0.3870    0.7739    1.1609    1.5478    1.9348    2.3217];
     self.w = [0.5962    0.2095    0.1473    0.0857    0.1028    0.0581    0.2072];
     self.d = [-6.2078   -3.8375   -1.3246         0    1.3246    4.1386    7.2146];
     self.ws = [-0.0477    0.0790    0.0844    0.6319   -0.0760   -0.0967   -0.0822];
     self.gdrift = zpk([0.9857, 0.09982, 0.999874129919322],...
       [0.923102847666306, 0.990415930939774, 0.998549077001794, 0.999904798380564],...
       0.0134772619280941, 40e-6);
   end
 end
 methods(Test)
   function verify_hyst_sat_drift_paral_inverse(self)
     % addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'));
     % modelFit_file = fullfile(PATHS.sysid, 'FRF_data', 'x-axis_sines_infoFourierCoef_5-30-2018-01.mat');
     % load(modelFit_file, 'modelFit');
     % hyst = modelFit.models.hyst_drift_paral;
     % gdrift = ss(hyst.gdrift);
     % D = gdrift.d;
     % gdrift.d = 0;
     Ts = self.gdrift.Ts;
     
     t = (0:5000)'*Ts;
     To = 1/50;
     ud = sin((2*pi/To)*t);
     U_des = timeseries(ud, t);
     
     [rp, wp] = PIHyst.invert_hyst_PI(self.r, self.w);
     [dp, wsp] = PIHyst.invert_sat(self.d, self.ws);
     w = self.w;
     r = self.r;
     d = self.d;
     ws = self.ws;
     
     trun = t(end);
     sim('test_hyst_drift_parallel')
     
     self.verifyEqual(U_des.Data, Y.Data, 'AbsTol', 1e-9);
     % figure(1); clf
     % plot(U_des.Time, U_des.Data)
     % hold on, grid on
     % plot(U_full.Time, U_full.Data);
     % plot(Y.Time, Y.Data, '--')

   end
 end

end
