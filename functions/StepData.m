classdef StepData
    
    properties
        params;
        file;
        fig_files;
        savedata;
        verbose;
        results;
    end
    
    methods
        function self = StepData(Params, varargin)
            p = inputParser;
            p.addParameter('verbose', 0)
            p.addParameter('savedata', true)
            p.addParameter('file', '')
            p.addParameter('fig_files', '')
            parse(p, varargin{:});
            
            
            self.params = Params;
            self.file = p.Results.file;
            self.fig_files = p.Results.fig_files;
            self.verbose =  p.Results.verbose;
            self.savedata =  p.Results.savedata;
            self.results = [];
        end
        
    end
    
end
