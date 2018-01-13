% upd = ProgressBar(max_iter, varargin)
% 
% progressbar(..., 'start_str', start_str) a label for the start of the
% progress bar
%
% progressbar(..., 'bar_len', 35) the number of characters wide that the 
% bar is
%
% prgressbar(..., 'fid', fid) A file id to write the progress to a
% file. Default is fid=1, which writes to standard terminal.
% 
% Returns
% -------
%   upd : a function handle that will update the progress bar.
% 
% Example:
% -------
% N = 1000
% upd = progressbar(N, 'start_str', 'my label')
%
% for i=1:N
%     upd(i)
% end
classdef ProgressBar < handle
    properties
        max_iter;
        start_str;
        bar_len;
        fid;
        prior_length;
        prior_msg;
        has_initialized;
        prior_perc_done;
    end
    
    methods
        function self = ProgressBar(max_iter, varargin)
            % self = ProgressBar(max_iter, varargin)
            p = inputParser;
            default_bar_len = 35;
            default_start_str = '';
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            addParameter(p, 'bar_len', default_bar_len, validScalarPosNum);
            addParameter(p, 'start_str', default_start_str);
            addParameter(p, 'fid', 1);
            
            parse(p, varargin{:})

            self.bar_len = p.Results.bar_len;
            self.start_str = p.Results.start_str;
            self.fid = p.Results.fid;
            self.max_iter = max_iter;
            
            bar_init = sprintf('[%s]', repmat(' ', 1, self.bar_len));
            whole_bar = sprintf('%s: %s%s', self.start_str, bar_init);
            self.prior_msg = sprintf('%s %3.1f', whole_bar, 0); 
            self.prior_length = length(self.prior_msg);
            self.has_initialized = 0;
        end
        function upd(self, iter)
            % Update the progress bar.
            fracDone = iter / self.max_iter;
            percentDone = round(100 * fracDone, 0);
            
            
            done_bar = repmat('+', 1, floor(self.bar_len * fracDone));
            undone_bar = repmat(' ', 1, self.bar_len - floor(self.bar_len * fracDone));
            whole_bar = sprintf('%s: [%s%s]', self.start_str, done_bar, undone_bar);
            msg = sprintf('%s %3.0f', whole_bar, percentDone); 

            if ~self.has_initialized || self.fid > 1
                fprintf(self.fid, '%s\n', msg);
                self.prior_msg = msg;
                self.prior_length = length(msg);
                self.has_initialized = 1;
            elseif percentDone > self.prior_perc_done;
                reverseStr = repmat(sprintf('\b'), 1, 1+length(self.prior_msg));
                fprintf(self.fid, '%s%s\n', reverseStr, msg);
                self.prior_msg = msg;
                self.prior_length = length(msg);
            end

            self.prior_perc_done = percentDone;
            
        end
    end
end