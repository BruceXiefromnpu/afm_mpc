classdef StepData
    
    properties
        params;
        file;
        fig_files;
        savedata;
        results;
    end
    
    methods
        function self = StepData(Params, varargin)
            p = inputParser;
            p.addParameter('savedata', true)
            p.addParameter('file', '')
            parse(p, varargin{:});
            
            
            self.params = Params;
            self.file = p.Results.file;
            self.savedata =  p.Results.savedata;
            self.results = [];
        end
        
    end
    
end
