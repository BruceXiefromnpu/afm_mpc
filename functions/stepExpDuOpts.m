classdef stepExpDuOpts
% This class holds the option set for the stepExp class. 
%
% Usage 1: Set each property with name value pairs
%   stepExpOpts('pstyle', '--k', 'name', 'afm mpc', ...)
%
% Usage 2: If the first argument is another stepExpOpts object, then all
% the properties of this object will be copied. Any other name value pairs
% will then modify those properties.
%   stepExpOpts(stepOpts1, 'pstyle', '-r'...)
    
    properties
        pstyle;
        name;
        exptype;
        params;
        controller;
        datanames;
        yunits;
        uunits;
        yscaling;
        step_ref;
        TOL;
    end
    
    methods
        function obj = stepExpDuOpts(varargin)
            % If first input is already a stepExpOpts, use that data as
            % default. 
            if isa(varargin{1}, 'stepExpDuOpts')
                oldOpts = varargin{1};
                props = properties(oldOpts); % or properties('foo') or fieldnames(myObj)
                opts = cell2struct(cellfun(@(prop) oldOpts.(prop), props, 'UniformOutput', false), props, 1);
                varargin(1) = [];
            else
            % otherwise, set standard defaults.
               step_ref_def = StepRef([0], 100);
               opts = struct('exptype', 'not specified', 'name', 'none', 'step_ref', step_ref_def,...
                       'TOL', 0.01, 'params', [], 'controller', [], 'pstyle', '-b',...
                       'datanames', {{'x', 'y'}}, 'yunits', '[v]', 'uunits', '[v]',...
                        'yscaling', 1);  % set default options            
            end
            % Update with any different options.
%             keyboard
           for pair = reshape(varargin, 2, [])
              inpName = pair{1};
              opts.(inpName) = pair{2};
           end
           
           % And copy into the class fields
            names = fieldnames(opts);
            try
              for name = names'
                obj.(name{1}) = opts.(name{1});
              end
            catch
              error(['You supplied an input name-value pair that is not',...
                'defined for this class. This is possibly because the ',...
                'stepExp class has been updated and your script is assuming ',...
                'the old function behaivior. For a quick fix, use ',...
                'stepExp_old instead\n';])
            end
        end
        
    end
    
end

