classdef EchoFile< handle
    properties
        fname;
    end
    methods 
        function self = EchoFile(fname, varargin)
            p=inputParser;
            p.addParameter('overwrite', true);
            parse(p, varargin{:});
            self.fname = fname;
            
            if p.Results.overwrite
                cmd = sprintf('echo > %s', self.fname);
                system(cmd);
            end
        end
        
        function [ status] = echo_file(self, varargin)
            if isempty(varargin)
                error('echo_file requires at least one argument')
            end

            if strcmp(varargin{1}, '-set')
                self.fname = varargin{2};
                status = 0;
                return
            elseif strcmp(varargin{1}, '-get')
                status = self.fname;
                return
            end

            str = sprintf(varargin{:});
            if self.fname == 0
                fprintf(1, '%s', str);
            else
                status = sys_echo(str, self.fname);
            end
        end
    end
end

function status = sys_echo(str, fname)
% echo adds a newline anyway, so strip any trailing newlines
    ind = regexp(str, '\n$');
    str = str(1:ind-1);
    cmd = sprintf('echo "%s" >> %s', str, fname);
    status = 1;
    while status
       status = system(cmd);
    end

end