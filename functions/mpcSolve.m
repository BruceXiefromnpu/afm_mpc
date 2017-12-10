function solveMpc(block)
% Level-2 MATLAB file S-Function for solving linear constrained MPC problem
% with classes defined in condensedMPC and FGMprob.

  setup(block);
end  

function setup(block)
  
  % Register dialog parameter: mpcProblem class
  block.NumDialogPrms = 1;
  block.DialogPrmsTunable = {'Nontunable'}; %This is critical to allow class parameter
  
  % Regieste number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;

  block.SampleTime = [-1 0];
  % Setup functional port properties to dynamically inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;

  block.InputPort(1).Complexity   = 'Real'; 
  block.InputPort(1).DataTypeId   = 0;
  block.InputPort(1).SamplingMode = 'Sample';
  block.InputPort(1).Dimensions   = block.DialogPrm(1).Data.ns;
  
  block.OutputPort(1).Complexity   = 'Real';
  block.OutputPort(1).DataTypeId   = 0;
  block.OutputPort(1).SamplingMode = 'Sample';
  block.OutputPort(1).Dimensions   = block.DialogPrm(1).Data.nu;

  
  % Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  % Register methods
  block.RegBlockMethod('CheckParameters',         @CheckPrms);
  %   block.RegBlockMethod('ProcessParameters',       @ProcessPrms);
  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
  block.RegBlockMethod('Start',                   @Start);  
  block.RegBlockMethod('Outputs',                 @Outputs);
  
  % Block runs on TLC in accelerator mode.
  block.SetAccelRunOnTLC(true);
end  


function CheckPrms(block)
  mpcProb = block.DialogPrm(1).Data;
  
  if ~isa(mpcProb, 'sparseMPCprob') && ~isa(mpcProb, 'condensedMPCprob') 
    error('Need condensedMPCprob or FGMprob data type.');
  end
end  


function DoPostPropSetup(block)
  % Setup Dwork  
  mpcProb = block.DialogPrm(1).Data;
  block.NumDworks = 1;
  block.Dwork(1).Name = 'ziyi'; 
  block.Dwork(1).Dimensions      = mpcProb.n_warmstart;
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = true;
  
  % Register all tunable parameters as runtime parameters.
  block.AutoRegRuntimePrms;
end


function ProcessPrms(block)
  block.AutoUpdateRuntimePrms;
end


function Start(block)
  
  % Initialize Dwork
  mpcProb = block.DialogPrm(1).Data;
  N = mpcProb.N_mpc;
  % Length of the total decision variable.

  if N>0
      block.Dwork(1).Data = zeros(mpcProb.n_warmstart, 1);
  else
      block.Dwork(1).Data = 0;
  end

end


function Outputs(block)
  mpcProb = block.DialogPrm(1).Data;
  n_warmstart = mpcProb.N_mpc;
  
  xk_1 = block.InputPort(1).Data;
    
%   ziyi = block.Dwork(1).Data;
  if n_warmstart>0
      [U] = mpcProb.solve(xk_1);
      u_k = U(1:mpcProb.nu);
%       block.Dwork(1).Data = double(ziyi);
  else
      u_k = mpcProb.solve(xk_1, ziyi);
      block.Dwork(1).Data = 0;
  end
  
  
  block.OutputPort(1).Data = double(u_k);

end  


