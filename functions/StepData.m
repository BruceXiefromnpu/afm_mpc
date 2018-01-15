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
        
        function plot_single_traj(self, index, ax1, ax2)
            if ~exist('ax1', 'var') || ~exist('ax2', 'var')
                fig = gcf();
                ax1 = subplot(211);
                ax2 = subplot(212);
            end
            
            self.plot_single_ytraj(index, ax1);
            self.plot_single_utraj(index, ax2);
        
        end
        
        
        function status = stepdata_struct_unchanged(self)
        % This function will return 1 if
        % (a) data_struct.file exists
        % (b) the params fields in the loaded data struct (at data_struct.file) 
        % and the provided data are the same. 
            
            if ~exist(self.file, 'file')
                status = 0;
                % keyboard
                return 
            else
                load(self.file)
                % Provides: step_data
            end
            params_other = step_data.params;
            if ~isequal(self.params, params_other)
                status = 0;
                return
            end
            
            status = 1;
        end


    end
    
end



