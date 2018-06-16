classdef TestHystDriftParallel < matlab.unittest.TestCase
  
  
 properties
   r;
   w;
   d
   ws
   gdrift;
   lv_data_path = ['C:\Users\arnold\Documents\matlab\afm_mpc_journal',...
     '\labview\UnitTests\fpga_harnesses\Sat_Hyst_TestCase\',...
     'hyst_drift_sat_parallel.csv'];
 end
 methods
   function self = TestHystDriftParallel()
     addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'));
     modelFit_file = fullfile(PATHS.sysid, 'FRF_data', 'x-axis_sines_infoFourierCoef_5-30-2018-01.mat');
     load(modelFit_file, 'modelFit');
     hyst = modelFit.models.hyst_drift_paral;
     self.r = hyst.r;
     self.w = hyst.w;
     self.d = hyst.d;
     self.ws = hyst.ws;
     self.gdrift = hyst.gdrift;
     
%      self.r  = [0    0.3870    0.7739    1.1609    1.5478    1.9348    2.3217]';
%      self.w = [0.5962    0.2095    0.1473    0.0857    0.1028    0.0581    0.2072]';
%      self.d = [-6.2078   -3.8375   -1.3246         0    1.3246    4.1386    7.2146]';
%      self.ws = [-0.0477    0.0790    0.0844    0.6319   -0.0760   -0.0967   -0.0822]';
%      
%      z = [0.984498204056153, 0.998186564672070, 0.999887866540897];
%      p =[0.656592343924253, 0.987479959610815, 0.998404249820942, 0.999904150810044];
%      k = 0.157983741634271;
%      self.gdrift = zpk(z, p, k, 40e-6);

   end
   
   function write_data_for_labview(self)
     
     % addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'));
     % modelFit_file = fullfile(PATHS.sysid, 'FRF_data', 'x-axis_sines_infoFourierCoef_5-30-2018-01.mat');
     % load(modelFit_file, 'modelFit');
     % hyst = modelFit.models.hyst_drift_paral;
     % gdrift = ss(hyst.gdrift);
     % D = gdrift.d;
     % gdrift.d = 0;
     Ts = self.gdrift.Ts;
     
     t = (0:500)'*Ts;
     To = 1/50;
     ud = sin((2*pi/To)*t)*4;
     U_des = timeseries(ud, t);
     
     [rp, wp] = PIHyst.invert_hyst_PI(self.r, self.w);
     rp = fi(rp, 1, 16, 11);
     wp = fi(wp, 1, 16, 11);
     [dp, wsp] = PIHyst.invert_sat(self.d, self.ws);
     dp = fi(dp, 1, 16, 11);
     wsp = fi(wsp, 1, 16, 11);
     
     w = self.w;
     r = self.r;
     d = self.d;
     ws = self.ws;
     gdrift = self.gdrift;
     [num, den] = tfdata(gdrift, 'v');
     trun = t(end);
     
     [num, den] = tfdata(gdrift, 'v');
      
     options = simset('SrcWorkspace','current');
     sim('test_hyst_drift_parallel', [], options)
     
     % Easiest to just re-load these, than convert back to double.
    [rp, wp] = PIHyst.invert_hyst_PI(self.r, self.w);
    [dp, wsp] = PIHyst.invert_sat(self.d, self.ws);
    
     fid = fopen(self.lv_data_path, 'w+');
     fprintf(fid, '%d, %d, %d\n', length(den)-1, length(rp), length(dp));
     write_vector(fid, num);
     write_vector(fid, den);
     write_vector(fid, rp);
     write_vector(fid, wp);
     write_vector(fid, dp);
     write_vector(fid, wsp);
     write_vector(fid, U_des.Data);
     write_vector(fid, U_full.Data);
     fclose(fid);
     
     figure(1); clf
     plot(U_des.Time, U_des.Data)
     hold on, grid on
     plot(U_full.Time, U_full.Data);
     plot(Y.Time, Y.Data, '--')
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
     gdrift = self.gdrift;
     
        
     trun = t(end);
     options = simset('SrcWorkspace','current');
     sim('test_hyst_drift_parallel', [], options)
     
     self.verifyEqual(U_des.Data, Y.Data, 'AbsTol', 1e-9);
     % figure(1); clf
     % plot(U_des.Time, U_des.Data)
     % hold on, grid on
     % plot(U_full.Time, U_full.Data);
     % plot(Y.Time, Y.Data, '--')

   end
 end

end

function write_vector(fid, vec)
  vec = vec(:)';
  fprintf(fid, '%.16f,', vec(1:end-1));
  fprintf(fid, '%.16f\n', vec(end));
  
end
