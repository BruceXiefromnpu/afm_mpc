function [ status ] = ssendmail(subject, message_file, varargin)
% securely (meaning my password is not stored in clear text in matlabs 
% config file) sends email.
% with my msmpt recipient, subject, message_file
%
% [ status ] = ssendmail(subject, message, varargin)
%
%              ssendmail('to', 'abraker@fastmail.com'); Default is
%              robr9299@colorado.edu. For multiple reciepients, this should
%              be a cell array of email addresses. 
%          
%              ssendmail('attachment', 'somefile.ext'). Currently not
%              implemented.
p = inputParser;
p.addParameter('to', 'robr9299@colorado.edu');
p.addParameter('attachments', '');
parse(p, varargin{:});
to = p.Results.to;
attachments = p.Results.attachments;

% msg = sprintf(['From: arnold@rabraker.com\n'...
%                'Subject: %s: %s\n',...
%                '%s\n'], date(), subject, message);
subject = sprintf('%s: %s', date, subject);
sndmail = '/usr/bin/msmtp -a arnold-rabraker';
if ~iscell(to)
    to = {to};
end
for to_ = to
    % cmd = sprintf('echo "%s" |%s %s', msg, sndmail, to_{1});
    cmd = sprintf('cat %s|./sendattach.sh -t %s -s "%s" ', message_file,...
                  to_{1}, subject);
                  
    cmd = add_attachments(cmd,  attachments);
    status = system(cmd);
end


end

function cmd = add_attachments(cmd, attachments)
    if ~isempty(attachments) % need this?
        for attach = attachments
            % k=1:length(attachments)
            % attach = attachments{k}
            % keyboard
            cmd = sprintf('%s -a %s', cmd, attach{1});
        end
    end

end
