% frfopts = frfbodeoptset('name1', val1, 'name2', val2...)
%
% Sets options for the frfBode. Will return a structure of defualt options
% if called with no arguments.
%
%   Options:
%       coherenceOn: 0, (no), 1 (yes). 1 will plot coherence.
%       magOnly: 1 will only plot the magnitude response. 0 will plot mag &
%                phase.
%       singleFigure: 1 will plot all MIMO input/output pairs into a single
%                     plot. 0 will plot each MIMO pair into its own figure.
%       color: 'b', 'g', 'r' etc.
%       lineStyle: '-', ':', '--' etc.


classdef frfbodeoptset

    
    properties
        coherenceOn;
        magOnly;
        singleFigure;
        color;
        lineStyle;
        ax;
    end
    
    methods
        function obj = frfbodeoptset(varargin)
            if isa(varargin{1}, 'frfbodeoptset')
                obj = varargin{1};
                varargin(1) = [];
            else
                defaults = struct('coherenceOn', 0, 'magOnly', 0, 'singleFigure', 1,...
                      'color', 'b', 'lineStyle', '-');
              for name = reshape(fieldnames(defaults), 1, [])% for goes across columns. fieldnames comes as rows.
%                   keyboard
                    % name is a cell array with a single entry string, which is
                    % not a valid struct dynamic reference. We have to do
                    % name{1} to get just the string.
                  obj.(name{1}) = defaults.(name{1});
              end
            end
              for pair = reshape(varargin, 2, [])
                 try 
                     obj.(pair{1}) = pair{2};
                 catch
                    error('No property named %s.\n Availible option fields are: %s', pair{1}, 'none') 
                 end
                  
              end
        end
    end
    
end

