%#codegen
function solveMpc_fxp(block)
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
%   block.InputPort(1).DataTypeId   = 0;
  block.InputPort(1).SamplingMode = 'Sample';
  block.InputPort(1).Dimensions   = block.DialogPrm(1).Data.ns;
  
  block.OutputPort(1).Complexity   = 'Real';
%   block.OutputPort(1).DataTypeId   = 0;
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
  % Inherit datatype id for non-built-in datatypes
  block.RegBlockMethod('SetInputPortDataType',    @SetInpPortDataType);
  block.RegBlockMethod('SetOutputPortDataType',   @SetOutputPortDataType);
  
  % Block runs on TLC in accelerator mode.
  block.SetAccelRunOnTLC(true);
end  
function SetOutputPortDataType(block, idx, dt)
  
  block.OutputPort(idx).DataTypeID = dt;
end


function SetInpPortDataType(block, idx, dt)

  block.InputPort(idx).DataTypeID = dt;
  % Obtain data types of input ports
  dt1 = block.InputPort(1).DataTypeID;
%   dt2 = block.InputPort(2).DataTypeID;
  
  % Calculate output data types when both input port data types are ready
  if dt1 ~= -1 && block.DataTypeIsFixedPoint(dt1)
            
%     p1 = block.FixedPointNumericType(dt1);
%     p2 = block.FixedPointNumericType(dt2);
%     
%     w1 = p1.WordLength;
%     f1 = p1.FractionLength;
%     s1 = p1.Signed;
    
    w1 = 32;
    f1 = 30;
    s1 = 1;
    
    newId = block.RegisterDataTypeFxpSlopeBias(s1, w1, 1.0*2^(-f1), 0, false);
    
    block.OutputPort(1).DataTypeID = newId;
  end
end


function CheckPrms(block)
  mpcProb = block.DialogPrm(1).Data;
  
  if ~isa(mpcProb, 'condensedMPCprob') && ~isa(mpcProb, 'FGMprob')
    error('Need condensedMPCprob or FGMprob data type.');
  end
end  


function DoPostPropSetup(block)
  % Setup Dwork  
  mpcProb = block.DialogPrm(1).Data;
  block.NumDworks = 1;
  block.Dwork(1).Name = 'ziyi'; 
  block.Dwork(1).Dimensions      = max(mpcProb.n_warmstart, 1);
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
  N = mpcProb.n_warmstart;
  if N>0
      block.Dwork(1).Data = zeros(N, 1);
  else
      block.Dwork(1).Data = 0;
  end

end


function Outputs(block)
  mpcProb = block.DialogPrm(1).Data;
  n_warmstart = mpcProb.n_warmstart;
  nu = mpcProb.nu;
  xk_1 = block.InputPort(1).Data;
    
  ziyi = block.Dwork(1).Data;
%   if mpcProb.n_warmstart>0
%       [u_k, ziyi] = mpcProb.solve(xk_1, ziyi);
    z_i = ziyi(1:n_warmstart/2);
    y_i = ziyi(n_warmstart/2+1:end);
    % shift zi and yi for warm start
    z_i = fi([z_i(nu+1:end); z_i(end-nu+1:end)], 1, mpcProb.nw, mpcProb.nf);
    y_i = fi([y_i(nu+1:end); y_i(end-nu+1:end)], 1, mpcProb.nw, mpcProb.nf);

    f = fi(mpcProb.ML*xk_1, 1, mpcProb.nw, mpcProb.nf);

    [z_i_1, y_i_1] = clqr_fgm_fxp_fict(mpcProb.I_HL, f, mpcProb.beta,...
                     mpcProb.maxIter, mpcProb.uMax, z_i, y_i);
    ziyi = [z_i_1; y_i_1];
    u_k = z_i_1(1);
      block.Dwork(1).Data = ziyi.data;
%   else
%       u_k = mpcProb.solve(xk_1, ziyi);
%       block.Dwork(1).Data = 0;
%   end
  
%   tmp = u_k
% keyboard
  block.OutputPort(1).Data = u_k;

end  


