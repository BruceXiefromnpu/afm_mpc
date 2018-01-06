% upd = progressbar(max_iter, varargin)
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
% upd = progressbar(N, 'start_str', 'my label'
%
% for i=1:N
%     upd(i)
% end

function upd = progressbar(max_iter, varargin)

    p = inputParser;
    default_bar_len = 35;
    default_start_str = '';
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    addParameter(p, 'bar_len', default_bar_len, validScalarPosNum);
    addParameter(p, 'start_str', default_start_str);
    
    
    parse(p, varargin{:})
    
    bar_len = p.Results.bar_len;
    start_str = p.Results.start_str;

    if ~exist('bar_len', 'var')
        bar_len = 35;
    end

    function update(iter)
       fracDone = iter / max_iter;
       percentDone = 100 * fracDone;
        
        done_bar = repmat('+', 1, floor(bar_len * fracDone));
        undone_bar = repmat(' ', 1, 35 - floor(bar_len * fracDone));
        whole_bar = sprintf('%s: [%s%s]', start_str, done_bar, undone_bar);
        msg = sprintf('%s %3.1f', whole_bar, percentDone); %Don't forget this semicolon
        
        if iter == 0
            reverseStr = sprintf('');
        else
            reverseStr = repmat(sprintf('\b'), 1, length(msg)+2);
        end
        fprintf('%s%s%%\n', reverseStr, msg);
        if iter >= max_iter
            fprintf('\n\n')
        end

    end

upd = @update;
end