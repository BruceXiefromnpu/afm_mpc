function [ output_args ] = buildOptions( input_args )
%BUILDOPTIONS Summary of this function goes here
%   Detailed explanation goes here
           opts = struct('textOn', 1, 'yunits', '[v]', 'yscaling', 1);  % set default options
           for pair = reshape(varargin, 2, [])
              inpName = pair{1};
              opts.(inpName) = pair{2};
           end
           
            

end

