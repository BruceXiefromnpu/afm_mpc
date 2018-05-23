classdef CanonPlants
  
  properties
    path = fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat');
    plants;
  end
  
  
  methods (Static)
    function plants = plants_with_drift_inv(with_hyst)
    % plants = plants_with_drift_inv(with_hyst)
    % Builds the standard versions of the plants. plants is a
    % structure with the following fields
    % 
    % plants.PLANT
    % plants.SYS
    % plants.sys_recyc
    % plants.sys_nodelay
    % plants.sys_recyc_nodelay
    % plants.gdrift
    % plants.gdrift_inv
    % plants.hyst
    % 
    % The model data is loaded from 
    % fullfile(PATHS.sysid, 'hysteresis/steps_hyst_model.mat')
    % if with_hyst == true
    % Otherwise, it is loaded from 
    % fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'))
    % 
    % Both all of the models are converted to a balanced
    % realization and then scaled. This is necessary for the PLANT
    % because in the compare_max_ref_XX scripts, we do not have an
    % observer, the PLANT and control law need the same state
    % space. 
      
      if with_hyst
        load(fullfile(PATHS.sysid, ['hysteresis/steps_hyst_model.mat']));
        hyst.rp = rp;
        hyst.wp = wp;
        hyst.r = r;
        hyst.w = w;
        gdrift_inv = 1/gdrift;  
        
        % Grab the current model.
        mf = load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'));
        plants.G_delu2powI = mf.modelFit.models.G_deluz2powI;
        
      else
        load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'));
        [Gvib, gdrift_inv] = eject_gdrift(modelFit.models.G_uz2stage);
        
        plants.G_delu2powI = modelFit.models.G_deluz2powI;
        
        gdrift = 1/gdrift_inv;
        hyst.r = [];
        hyst.w = [];
        hyst.rp = [];
        hyst.wp = [];
      end
      
      
      %PLANT = ss(Gvib);
      SYS = ss(Gvib);
      
      plants.gdrift = gdrift;
      
      plants.gdrift_inv = gdrift_inv;
      
      SYS = balreal(SYS);
      Nx = SSTools.getNxNu(SYS);
      T = diag(1./Nx)/10;
      SYS = ss2ss(SYS, T);
      PLANT = SYS;
      
      Nd = 9;
      SYS.iodelay = 0;
      SYS.InputDelay = Nd;
      
      plants.sys_nodelay = SYS;
      
      SYS = absorbDelay(SYS);
      PLANT.InputDelay = Nd;
      PLANT = absorbDelay(PLANT);
      
      plants.PLANT = PLANT;
      plants.SYS = SYS;
      plants.hyst = hyst;
      
      plants.sys_recyc=SSTools.deltaUkSys(SYS);
      plants.sys_recyc_nodelay=SSTools.deltaUkSys(plants.sys_nodelay);
      plants.Nd = Nd;
      
    end

    function plants_with_drift_internal()
      error(['This function is not yet implemented. In particular, ' ...
             'the experimental labview code is not setup to handle ' ...
             'this the drift inside the main system due to number ' ...
             'of states required'])
      
    end

    
  end
  
  
  
end
