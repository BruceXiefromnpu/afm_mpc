function mpcSolve(block)
% Level-2 MATLAB file S-Function for solving linear constrained MPC problem
% with classes defined in condensedMPC and FGMprob.

  setup(block);
end  

function setup(block)
  
  % Register dialog parameter: mpcProblem class
  block.NumDialogPrms = 1;
  block.DialogPrmsTunable = {'Nontunable'}; %This is critical to allow class parameter
  
  % Register number of input and output ports
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
  % block.RegBlockMethod('ProcessParameters',       @ProcessPrms);
  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
  block.RegBlockMethod('Start',                   @Start);  
  block.RegBlockMethod('Outputs',                 @Outputs);
  
  % Block runs on TLC in accelerator mode.
  block.SetAccelRunOnTLC(true);
end  


function CheckPrms(block)
% make sure we know how to handle the parameter that gets passed.
  mpcProb = block.DialogPrm(1).Data;
  ok_objs = {'sparseMPCprob', 'sparseMPCprob_OA',...
             'condensedMPCprob', 'condensedMPCprob_OA',...
             'condensedMPCprob_QP', 'FGMprob_fxp_1', 'FGMprob_1'};
  if ~ismember(class(mpcProb), ok_objs)
    str = sprintf(repmat('%s, ', 1, length(ok_objs)), ok_objs{:});
    error('Need one on %s data type.', str);
  end
end  


function DoPostPropSetup(block)
  % Register all tunable parameters as runtime parameters.
  block.AutoRegRuntimePrms;
end


function ProcessPrms(block)
  block.AutoUpdateRuntimePrms;
end


function Start(block)
  
% Used to initialize Dwork here for ziyi hotstart. Now, store hotstart data
% inside handle class property.
end


function Outputs(block)
  mpcProb = block.DialogPrm(1).Data;
  
  xk_1 = block.InputPort(1).Data;
    
  [U] = mpcProb.solve(xk_1);
  u_k = U(1:mpcProb.nu);
  
  block.OutputPort(1).Data = double(u_k);

end  


