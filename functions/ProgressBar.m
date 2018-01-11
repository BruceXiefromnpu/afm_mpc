% upd = ProgressBar(max_iter, varargin)
% 
% progressbar(..., 'start_str', start_str) a label for the start of the
% progress bar
%
% progressbar(..., 'bar_len', 35) the number of characters wide that the 
% bar is
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
        prior_length;
        prior_msg;
        has_initialized;
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

            parse(p, varargin{:})

            self.bar_len = p.Results.bar_len;
            self.start_str = p.Results.start_str;
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
            percentDone = 100 * fracDone;

            done_bar = repmat('+', 1, floor(self.bar_len * fracDone));
            undone_bar = repmat(' ', 1, self.bar_len - floor(self.bar_len * fracDone));
            whole_bar = sprintf('%s: [%s%s]', self.start_str, done_bar, undone_bar);
            msg = sprintf('%s %3.1f', whole_bar, percentDone); 

            if ~self.has_initialized
                fprintf('%s\n', msg);
                self.prior_msg = msg;
                self.prior_length = length(msg);
                self.has_initialized = 1;
            else
                reverseStr = repmat(sprintf('\b'), 1, 1+length(self.prior_msg));
                fprintf('%s%s\n', reverseStr, msg);
                self.prior_msg = msg;
                self.prior_length = length(msg);
            end

            %         if iter >= max_iter
            %             fprintf('\n')
            %         end
        end
    end
end