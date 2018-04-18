classdef frdExt_local < frd
    %FRFEXT Summary of this class goes here
    %   Detailed explanation goes here
    % 
    %    frdExt(respData, freqs, Ts).
    % Calling sequence is exactly the same as frd
    % 
    % You can also include coherence data:
    %   frdExt(respData, freqs, CoherenceData, Ts)
    % Methods:
    %   frfBode:  overloads standard frfBode function. Always plots into
    %   units of Hz.
    %   [hMag, hPhase] = frfBode(obj, F1, plotstyle)
    %
    % See Also: frd, frfBode
    properties
        Coherence;
                
    end
    
    methods
        function obj = frdExt_local(varargin)
             if length(varargin{3}) == length(varargin{1})
                % The only reason the third argument should be "long", ie
                % same length as the response data is if the user included
                % coherence infor
                Coherence = varargin{3};
                
                varargin(3) = [];
            else
                Coherence = [];
             end
            
            [Gz, w_s, ind_s] = monotonicFRF(varargin{1}, varargin{2});
            Coherence = monotonicCoherence(Coherence, ind_s);
            
            obj@frd(Gz, w_s, varargin{3:end});

            % In the single-axis case, Coherence may come in as a
            % vector, rather than 3d matrix. Test for and fix that.

            if size(Coherence, 3) ~= size(obj.ResponseData, 3) && ~isempty(Coherence)
               Coherence = reshape(Coherence, [], size(obj.ResponseData,2),...
                           size(obj.ResponseData, 3));
            end                                
         obj.Coherence = Coherence;

        end
        function H = removeDelay2(obj, Nd) % Removes delay from the frf. Assumes delay is constant for all inputs to outputs.
            sz = size(obj.ResponseData);
            no = sz(1);
            ni = sz(2);
            Ts = obj.Ts;
            gdelay = (exp(1i*obj.Frequency*Ts)).^Nd;
            max(max(max(obj.Coherence(:,:,:))))
            for numIn =1:ni
               for numOut=1:no
                   obj.ResponseData(numOut, numIn,:) = squeeze(obj.ResponseData(numOut, numIn,:)).*gdelay;
               end
            end
            H = obj;
                        max(max(max(H.Coherence(:,:,:))))
        end
        % ==================================================================
        function PlotHandles = frfBode(obj, F_s, plotstyle)
            sz = size(obj.ResponseData);
            no = sz(1);
            ni = sz(2);
            if ~exist('F1', 'var') || isempty('F1')
                for i =1:no+ni
                    F_s(i) = figure(100+i);
                end
            end
            if ~exist('plotstyle', 'var') || isempty('plotstyle')
                plotstyle = 'b';
            end
            
            tmpSys = chgFreqUnit(obj, 'Hz');
            if isempty(obj.Coherence)
                cohFlag = 0;
            else
                cohFlag = 1;
            end

%             subsBase = 100*ni+10*no;
            cnter = 1;
           
            for inNum=1:ni
                for outNum=1:no
                    Subs = [211;212];
%                     Subs = subsBase + cnter;
                    PlotHandles(cnter) = frfBode_internal(tmpSys, inNum,outNum, cohFlag, F_s(cnter),Subs, plotstyle);
                    
                    cnter = cnter+1;
                end
            end
            
        end
        % ==================================================================
        function PlotHandles = frfBodeMag(obj, F1, plotstyle)
            if ~exist('F1', 'var') || isempty('F1')
                F1 = figure(100);
            end
            if ~exist('plotstyle', 'var') || isempty('plotstyle')
                plotstyle = 'b';
            end
            
            tmpSys = chgFreqUnit(obj, 'Hz');
            if isempty(obj.Coherence)
                cohFlag = 0;
            else
                cohFlag = 1;
            end
            sz = size(obj.ResponseData);
            no = sz(1);
            ni = sz(2);
            subsBase = 100*ni+10*no;
            cnter = 1;
           
            for inNum=1:ni
                for outNum=1:no
%                     Subs = [211;212];
                    Subs = subsBase + cnter;
                    PlotHandles = frfMag_internal(tmpSys, inNum,outNum, cohFlag, F1,Subs, plotstyle);
                                  
                    cnter = cnter+1;
                end
            end
            
        end
        % ==================================================================        
        function H = interpFRF(obj, w_s, forceZeroPhase)
            if ~exist('forceZeroPhase')
                forceZeroPhase = 0
            end
            sz = size(obj.ResponseData);
            no = sz(1);
            ni = sz(2);
            for numIn=1:ni
                for numOut=1:no
                    H_frf(numOut,numIn,:) = interp1(obj.Frequency, squeeze(obj.ResponseData(numOut,numIn,:)), w_s, 'spline');
                    if forceZeroPhase
                        H_frf(numOut,numIn,1) = abs(H_frf(numOut,numIn,1));  %At dc, imaginary part should be zero.
                    end
                    Coh_i= interp1(obj.Frequency, squeeze(obj.Coherence(numOut,numIn,:)), w_s, 'spline');
                    Coh_i = max(0, Coh_i);
                    Coh_i = min(1, Coh_i);
                    Coh(numOut,numIn,:) = Coh_i;
                end
            end
             H = frdExt_local(H_frf, w_s,Coh, obj.Ts);
            H.InputName = obj.InputName;
            H.OutputName = obj.OutputName;

        end
        
    end
    
end

%                             (tmpSys, inNum,outNum, cohFlag, F_s(cnter),Subs, plotstyle);
function [PlotHandles] = frfBode_internal(SYS,inNum,outNum, cohFlag, F1,subs, plotStyle)

Gz_frf = squeeze(SYS.ResponseData(outNum,inNum,:));
w_s    = SYS.Frequency;
coherence_i = squeeze(SYS.Coherence(outNum, inNum, :));
if ~exist('F1', 'var') || isempty(F1)
    F1 = figure;
end

if ~exist('plotStyle', 'var') || isempty(plotStyle)
    plotStyle = 'b';
end




magsDB = 20*log10(abs(Gz_frf));
phase  = unwrap(angle(Gz_frf))*180/pi;

figure(F1);
subplot(subs(1,:));
stit = sprintf('G_{%s,%s}', SYS.InputName{inNum}, SYS.OutputName{outNum});

    if ~cohFlag
        [h_mag] = semilogx(w_s, magsDB, plotStyle, 'LineWidth', 1.5);
        PlotHandles.h_mag = h_mag;
        grid on; hold on;
        xlim([w_s(1), w_s(end)]);
        title(stit)
    else
        
        [AX, h_mag, h_coh] = plotyy(w_s, magsDB, w_s, coherence_i, 'semilogx', 'semilogx');

        grid on; hold on;
        AX(1).XLim = [w_s(1), w_s(end)];
        AX(2).XLim = [w_s(1), w_s(end)];
%         AX(2).YLim = [0,2];
        AX(1).YLim = [AX(1).YLim(1), max(magsDB)*1.1];
        [style, color] = splitLineStyle(plotStyle);
        if ~isempty(style)
        h_mag.LineStyle = style;
        end
        if ~isempty(color)
           h_mag.Color = color;
           AX(1).YColor = color;
        end
        title(AX(1), stit)
        ylabel(AX(1), 'Mag [dB]');
        ylabel(AX(2), 'Coherence');
        PlotHandles.h_mag = h_mag;
        PlotHandles.h_coh = h_coh;
        PlotHandles.AX    = AX;
        
        
    end
subplot(subs(2,:));
    h_phase = semilogx(w_s, phase, plotStyle, 'LineWidth', 1.5);
    grid on; hold on;
    xlim([w_s(1), w_s(end)]);
% if exist('freqUnit', 'var')
    xlabel(SYS.FrequencyUnit)
% end
    PlotHandles.h_phase = h_phase;
    
    ylabel('Phase [deg]')
    PlotHandles.F = F1;
end
% --------------------------------------------------------------------------
%                        frfMag_internal(tmpSys, j,i, F1,Subs, plotstyle);
function [PlotHandles] = frfMag_internal(SYS, inNum, outNum, cohFlag, F1,subs, plotStyle)

Gz_frf = squeeze(SYS.ResponseData(outNum,inNum,:));
w_s    = SYS.Frequency;
coherence_i = squeeze(SYS.Coherence(outNum,inNum,:));

if ~exist('F1', 'var') || isempty(F1)
    F1 = figure;
end

if ~exist('plotStyle', 'var') || isempty(plotStyle)
    plotStyle = 'b';
end




magsDB = 20*log10(abs(Gz_frf));
% phase  = unwrap(angle(Gz_frf))*180/pi;

figure(F1);
subplot(subs);
stit = sprintf('G_{%s,%s}', SYS.InputName{inNum}, SYS.OutputName{outNum});

    if ~cohFlag
        [h_mag] = semilogx(w_s, magsDB, plotStyle, 'LineWidth', 1.5);
        PlotHandles.h_mag = h_mag;
        grid on; hold on;
        xlim([w_s(1), w_s(end)]);
        'll';
    else
        
        [AX, h_mag, h_coh] = plotyy(w_s, magsDB, w_s, coherence_i, 'semilogx', 'semilogx');

        grid on; hold on;
        AX(1).XLim = [w_s(1), w_s(end)];
        AX(2).XLim = [w_s(1), w_s(end)];
%         AX(2).YLim = [0,2];
        AX(1).YLim = [AX(1).YLim(1), max(magsDB)*1.1];
        [style, color] = splitLineStyle(plotStyle);
        if ~isempty(style)
        h_mag.LineStyle = style;
        end
        if ~isempty(color)
           h_mag.Color = color;
           AX(1).YColor = color;
        end

        title(AX(1), stit)
        ylabel(AX(1), 'Mag [dB]');
        ylabel(AX(2), 'Coherence');
        PlotHandles.h_mag = h_mag;
        PlotHandles.h_coh = h_coh;
        PlotHandles.AX    = AX;
        
        
    end

    xlabel(SYS.FrequencyUnit)


end

function Coherence = monotonicCoherence(Coherence, ind_s)


    if ~isempty(Coherence)
        if length(size(Coherence)) > 2 % Check if Gz_frf contains mimo data. We have to 
                                        % handle that case differently.
            Coherence(:,:,ind_s) = [];
        else
            Coherence(ind_s) = []; 
        end
    end;

end


function [style, color] = splitLineStyle(st)
    inds = isletter(st);
    if inds(end) == 1
        % contains color
        color = st(end);
        style = st(1:end-1);
    else
        % no color specified
       style = st;
       color = [];
    end
    


end




