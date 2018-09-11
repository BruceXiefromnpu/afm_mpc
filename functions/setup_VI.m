% [e, vi ] = setup_VI(vipath, verbose, varargin )
% given a path to the VI, sets a variable number of front panel controls,
% and returns a reference to the vi and activeX server.
%
%  Inputs:
%   vipath: path to the VI instrument to run.
%   varargin: name value pairs, 'name of control', val1...
%
% Outputs:
%   e: handle? to the activex server connection. Or something like that.
%   vi: virtual instrument connection object. Use, vi.Run to run the VI.
% 
function [e, vi ] = setup_VI(vipath, verbose, varargin )
    e  = actxserver('LabVIEW.Application');
    vi = invoke(e, 'GetVIReference', vipath);
    assert(isa(verbose, 'logical'));
    
    for pair = reshape(varargin, 2, [])
        try
            if iscell(pair{2}) && verbose
                printfcell(pair{1}, pair{2})
            elseif verbose
              if length(pair{2}) ==1
                fprintf('control: %s,  val: %f\n', pair{1}, pair{2});
              else
                fprintf('control: %s,  val: MATRIX DATA', pair{1});
              end
            end
            vi.SetControlValue(pair{1}, pair{2}) 
        catch ME
            fprintf('Error for name value pair %s, \n', pair{1})
            keyboard
            rethrow(ME); 
        end
    end
    fprintf('\n');
end


function printfcell(name, val_cell)

    for k=1:length(val_cell)
%         keyboard
        fprintf('control{%d}: %s, val: %f\n', k, name, val_cell{k});
    end
end

