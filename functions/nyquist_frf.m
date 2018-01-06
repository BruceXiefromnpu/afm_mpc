function [ output_args ] = nyquist_frf(G_frf, w_s)
%NYQUIST_FRF Summary of this function goes here


plot(real(G_frf), imag(G_frf), '-')
hold on
plot(real(G_frf), -imag(G_frf), '-')
plot(-1, 0, 'r+')

dcm = datacursormode(gcf);
% datacursormode on;

set(dcm, 'UpdateFcn', {@local_callback_function, w_s})

end

function output_txt = local_callback_function(~,event_obj, w_s)
    % Display the position of the data cursor
    % obj          Currently not used (empty)
    % event_obj    Handle to event object
    % w_s          list of frequencies the nyquiest is evaluated at.
    
    pos = get(event_obj,'Position');
    
    % Why doesn't this show up if we do event_obj 'enter' in command
    % window??
    ind = get(event_obj, 'DataIndex'); 
    freq = w_s(ind)/2/pi;
    
    
    output_txt = {['Re: ',num2str(pos(1),4)],...
        ['Im: ',num2str(pos(2),4)],...
        ['Frequency: ', num2str(freq),' Hz']};


    % If there is a Z-coordinate in the position, display it as well
    if length(pos) > 2
        ind_str = output_txt{end}
        output_txt{end} = ['Z: ',num2str(pos(3),4)];
        output_txt{end+1} = ind_str;
    end

end