classdef StepData
    
    properties
        params;
        file;
        savedata;
        results;
        logger;
        ProgBar
    end
    
    methods
        function self = StepData(Params, varargin)
            p = inputParser;
            p.addParameter('savedata', true)
            p.addParameter('file', '')
            p.addParameter('logger', @fprintf);
            p.addParameter('ProgBar', @ProgressBar);
            parse(p, varargin{:});
            
            self.params = Params;
            self.file = p.Results.file;
            self.savedata =  p.Results.savedata;
            self.logger = p.Results.logger;
            self.ProgBar = p.Results.ProgBar;
            
            self.results = [];
        end
        
        function plot_single_traj(self, ref_idx, ax1, ax2, varargin)
        % plot_single_traj(self, ref_idx, ax1, ax2, varargin)    
            if ~exist('ax1', 'var') || ~exist('ax2', 'var')
                fig = gcf();
                ax1 = subplot(211);
                ax2 = subplot(212);
            end
            
            self.plot_single_ytraj(ref_idx, ax1, varargin{:});
            self.plot_single_utraj(ref_idx, ax2, varargin{:});
        end

        function status = stepdata_struct_unchanged(self)
        % This function will return 1 if
        % (a) data_struct.file exists
        % (b) the params fields in the loaded data struct (at data_struct.file) 
        % and the provided data are the same. 
            warning('off', 'MATLAB:structOnObject');
            if ~exist(self.file, 'file')
                status = 0;
                self.logger('self.file does not exist!')
                return 
            else
                load(self.file);
                % Provides: step_data
            end
            params_other = struct(step_data.params); 
            tmp_self_params = struct(self.params); % Do this in case I ever make it handle class.
            
            if isfield(tmp_self_params, 'Q')
                if ~isfield(params_other, 'Q')
                    status = 0;
                    self.logger('Matrix Q does not exist!')
                    return
                end
                Qstatus = isMatrixAlmostEqual(self.params.Q, params_other.Q);
                tmp_self_params = rmfield(tmp_self_params, 'Q');            
                params_other = rmfield(params_other, 'Q');
            else
                Qstatus = 1;
            end
            
            if ~isequal(tmp_self_params, params_other)
                status = 0;
                return
            end
            
            status = 1 & Qstatus;
        end


    end
    
end



function status = isMatrixAlmostEqual(A, B, varargin)
%  status = isMatrixAlmostEqual(A, B, varargin)
% quick check to see if matrices A and B are almost equal. Computed is 
% max(max(abs(A-B)))) >? tol, tol={0.001 | varargin{1} }

    if length(varargin) == 1
        tol = varargin{1}
    else
        tol =0.001;
    end
    if (size(A,1) ~= size(B,1) ) || ( size(A,2) ~= size(B,2) )
      status = 0;
      return
    end
    q = A-B;
    maxq = max(max(abs(q)));
    fprintf('A, B differ by %f\n', maxq)
    if maxq > tol
        status = 0;
    else
        status = 1;
    end

end

