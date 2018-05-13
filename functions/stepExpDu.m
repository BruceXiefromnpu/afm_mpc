% This is a class to house (hopefully) all relevant data and parameters
% for a step experiment. 
%   properties:
%       y, u, du:    output, control input, diff(u). Should be y.Time, y.Data
%       pstyle:  e.g., '--g'
%       name:    e.g. 'linear feedback (sim)'. This field is used to
%           construct the legend, ie, h.DisplayName = obj.name
%       exptype: e.g., 'linear, sim', 'CLQR, exp', 'MPC'
%       datanames: names for the outputs, y. 'x', or {{'x', 'y'}}. These
%       fields are used to construct the y-axis labels
%       yunits:  string of units for output, '[v]' 
%       uunits:  string of units for input,  '[v]'
%       y_ref:   vector of references: [1.5; 1]
%       TOL:    Settling tolerance, 0.01
%
%       controller: eg, a linear fdbk struct, ptos structure, etc.
%                  this is optional for backward compatibility. If not
%                  included, obj.controller = [];
%
%   Construction: y & u are the only required fields. The rest are including
%   with name value pairs, e.g., stepExp(y, u, 'pstyle', '--k',
%   'controller', 'K_fdbk'). All negelected fields are set with defaults:
%       pstyle: '-b'
%       'exptype': 'not specified'
%       'name':    'none'
%       'controller': []
%       'paraps': []
%     Defaults are set in the class stepExpOpts. You can construct all
%     these options from that class. You can use the same option set, and
%     pass in only a modification:
%       stepExp(y2, u2, stepOpts1, 'pstyle, '--k'), will apply the same options that
%       are contained in stepOpts1 but modify the 'pstyle' field to '--k'
% 
%
%   Methods:
%   [H_line, F1, H_string]= plot(sepExp_1, varargin). Plots the output y, and control u. 
%       example: plot(stepExp_1, F10, 'yscaling', 5, 'yunits_scaled', '[$\m um]$)
%       example: plot(stepExp_1, 'yscaling', 5, 'yunits_scaled', '[$\m um]$)
%
%            You can supply a figure handle as the first optional input
%            argument. The remaining arguments are supplied with name-value
%            pairs, given here with defaults: 
%            'textOn', 0: If this is set to 1, the settle times are printed in the bottom subplot. 
%            'yscaling', 1:  1 is the scale factor, the function plots
%                         y*yscaling. Eg, yscaling = 50micrometers/10v   
%           'yunits_scaled', '[v]':  New Unit string for the scaled units.
%
%       Outputs:
%           H_line: an array of line handles. The DisplayName property of
%           each handle is preset to obj.name. 
%           F1: Figure handle to whatever figure the function plotted to.
%           H_string: a handle to the strings printed. If 'textOn' = 0,
%           this should be empty.
%
%
%   [H, F1, H_s] = plotZoom(obj, xlm, ylm, 'textOn', 0,...
%                  'yunits_scaled', '[v]', 'showRefInfo', 0) 
%       The same inputs and outputs as before. The only change is in
%       addition to u(t) and y(t) a third subplot of y(t) is plotted around
%       around the setpoint to examine the steady state properties. Define
%       the zoom region with xlm, and ylm.
%
%   settle_time(stepExp_1) yields the settle time. This is calculated every
%   time the method is called, based on the current value held in
%   params.TOL.
%
%
%
%
classdef stepExpDu
  
  properties
    y;
    u;
    du;
    pstyle;
    name;
    exptype;
    params;
    controller;
    datanames;
    yunits;
    uunits;
    yscaling;
    y_ref;
    TOL;
    
  end
  
  methods
    function obj = stepExpDu(y, u, du, varargin)
      opts = stepExpOpts(varargin{:});
      
      obj.y = y;
      obj.u = u;
      obj.du = du;
      names = fieldnames(opts);
      for name = names'
        obj.(name{1}) = opts.(name{1});
      end
    end % end constructor
    % =================================================================
    function [H_line, F1, H_text] = plot(obj,  varargin)
      p = inputParser();
      p.addParameter('textOn', false);
      p.addParameter('yunits_scaled', obj.yunits);
      p.addParameter('lineWidth', 1);
      p.addParameter('umode', 'both');
      if length(varargin)>= 1 && isa(varargin{1}, 'matlab.ui.Figure')
        tmp = varargin{1};
        varargin(1) = [];
        if  isvalid(tmp)
          F1 = tmp;
        else
          F1 = figure;
        end
      else
        F1 = figure;
      end
      p.parse(varargin{:});
      preopts = p.Results;
      
      nChan = size(obj.y, 2);
      
      opts = buildStepExpPlotOpts(obj, preopts);
      if opts.textOn
        opts.TS_s = settle_time(obj);
      end
      if strcmp(opts.umode, 'u')
        [H_line, H_text]= plot_internal(obj.y, obj.u, opts, F1);
      elseif strcmp(opts.umode, 'du')
        [H_line, H_text]= plot_internal(obj.y, obj.u, opts, F1);
      elseif strcmp(opts.umode, 'both')
        [H_line, H_text]= plot_both_internal(obj.y, obj.u, obj.du, opts, F1);
      end
      
      
      
    end % end standard plot function
    
    % =================================================================
    function [H_line, F1, H_text] = plotZoom(obj, xlm, ylm, varargin)
      if length(varargin)>= 1
        if isgraphics(varargin{1},'figure')
          F1 = varargin{1};
          varargin(1) = [];
        else
          F1 = figure;
        end
      else
        F1 = figure;
      end
      preopts = struct('textOn', 0, 'yunits_scaled', obj.yunits,...
        'showRefInfo', 0);  % set default options
      for pair = reshape(varargin, 2, [])
        inpName = pair{1};
        preopts.(inpName) = pair{2};
      end
      opts = buildStepExpPlotOpts(obj, preopts);
      opts.xlm = xlm;
      opts.ylm = ylm;
      if opts.textOn
        opts.TS_s = settle_time(obj);
      end
      %             keyboard
      [H_line, H_text, H_zoomSubplots] = plotZoom_internal(obj.y, obj.u, opts, F1);
      
      if opts.showRefInfo
        plotSettleBoundary(obj, opts,F1, H_zoomSubplots)
      end
    end % end plotZoom function
    % =================================================================
    function TS_s = settle_time(obj)
      for kk = 1:size(obj.y.Data, 2)
        TS_s(kk) = settle_time(obj.y.Time, obj.y.Data(:,kk)*obj.yscaling,...
          obj.y_ref(kk)*obj.yscaling, obj.TOL*obj.y_ref(kk)*obj.yscaling);
        
      end
    end
    
  end %end methods
  
  
end

% =========================================================================
% =========================================================================
function [H_line, H_text] = plot_both_internal(y, u, du, opts, F1)
  no = size(y.Data, 2);
  figure(F1);
  H_line = [];
  subs_base = 300+no*10; % 2,no_col,
  for kk = 1:no % across columns of figure, no inputs, square system only.
    
    subplot(subs_base+kk); hold on; grid on;
    h1 = plot(y.Time, y.Data(:,kk)*opts.yscaling, opts.pstyle, 'LineWidth', opts.lineWidth);
    s_ylab = sprintf('%s-axis %s', opts.datanames{kk}, opts.yunits);
    ylabel(s_ylab, 'interpreter', 'latex');
    xlabel('t [s]', 'interpreter', 'latex');
    
    subplot(subs_base+no+kk); hold on; grid on;
    h2 = plot(u.Time, u.Data(:,kk), opts.pstyle, 'LineWidth', opts.lineWidth);
    s_ulab = sprintf('input: $u_{%s}$ %s', opts.datanames{kk}, opts.uunits);
    ylabel(s_ulab, 'interpreter', 'latex');
    xlabel('t [s]', 'interpreter', 'latex');
    
    subplot(subs_base+2*no+kk); hold on; grid on;
    h3 = plot(du.Time, du.Data(:,kk), opts.pstyle, 'LineWidth', opts.lineWidth);
    s_ulab = sprintf('$\\Delta u_{%s}$ %s', opts.datanames{kk}, opts.uunits);
    ylabel(s_ulab, 'interpreter', 'latex');
    xlabel('t [s]', 'interpreter', 'latex');
    
    h3.DisplayName = opts.name;
    h2.DisplayName = opts.name;
    h1.DisplayName = opts.name;
    
    H_line = [H_line,[h1;h2; h3]];
    H_text = [];
    if opts.textOn
      if ~isnan(opts.TS_s(kk))
        s1 = sprintf('%s settle time: %f', opts.name, opts.TS_s(kk));
      else
        s1 = sprintf('%s settle time: n/a', opts.name);
      end
      ss1 = catTexts(s1, F1);
      H_text = text(.3, .1, ss1, 'Units', 'normalized', 'VerticalAlignment', 'bottom');
    end
  end
end

% =========================================================================
% =========================================================================
function [H_line, H_text] = plot_internal(y, u, opts, F1)
  no = size(y.Data, 2);
  figure(F1);
  H_line = [];
  subs_base = 200+no*10; % 2,no_col,
  for kk = 1:no % across columns of figure, no inputs, square system only.
    
    subplot(subs_base+kk); hold on; grid on;
    h1 = plot(y.Time, y.Data(:,kk)*opts.yscaling, opts.pstyle, 'LineWidth', opts.lineWidth);
    s_ylab = sprintf('%s-axis %s', opts.datanames{kk}, opts.yunits);
    ylabel(s_ylab, 'interpreter', 'latex');
    xlabel('t [s]', 'interpreter', 'latex');
    
    subplot(subs_base+no+kk); hold on; grid on;
    h2 = plot(u.Time, u.Data(:,kk), opts.pstyle, 'LineWidth', opts.lineWidth);
    s_ulab = sprintf('input: $u_{%s}$ %s', opts.datanames{kk}, opts.uunits);
    ylabel(s_ulab, 'interpreter', 'latex');
    xlabel('t [s]', 'interpreter', 'latex');
    
    h2.DisplayName = opts.name;
    h1.DisplayName = opts.name;
    
    H_line = [H_line,[h1;h2]];
    H_text = [];
    if opts.textOn
      if ~isnan(opts.TS_s(kk))
        s1 = sprintf('%s settle time: %f', opts.name, opts.TS_s(kk));
      else
        s1 = sprintf('%s settle time: n/a', opts.name);
      end
      ss1 = catTexts(s1, F1);
      H_text = text(.3, .1, ss1, 'Units', 'normalized', 'VerticalAlignment', 'bottom');
    end
  end
end

% =========================================================================
% =========================================================================

function [H_line, H_text,  H_zoomSubplots] = plotZoom_internal(y, u, opts, F1)
    no = size(y.Data, 2);
    figure(F1);
    H_line = gobjects(3, no);
    subs_base = 300+no*10; % 2,no_col, 
    for kk = 1:no % across columns of figure, no inputs, square system only. 
        
        subplot(subs_base+kk); hold on; grid on;
            h1 = plot(y.Time, y.Data(:,kk)*opts.yscaling, opts.pstyle);
            s_ylab = sprintf('%s-axis %s', opts.datanames{kk}, opts.yunits);
            ylabel(s_ylab, 'interpreter', 'latex');
            xlabel('t [s]', 'interpreter', 'latex');

         H_zoomSubplots(kk) = subplot(subs_base+no+kk); hold on; grid on;
            h2 = plot(y.Time, y.Data(:,kk)*opts.yscaling, opts.pstyle);
            s_zoomylab = sprintf('(zoomed) %s-axis %s', opts.datanames{kk}, opts.yunits);
            ylabel(s_zoomylab, 'interpreter', 'latex');
            xlabel('t [s]', 'interpreter', 'latex');
            if opts.ylm(kk,1) ==0 && opts.ylm(kk,2) ==0
                ylim([-0.01, 0.01])
            else
                ylim(opts.ylm(kk,:));
            end
            xlim(opts.xlm(kk,:));
        subplot(subs_base+2*no+kk); hold on; grid on;
            h3 = plot(u.Time, u.Data(:,kk), opts.pstyle);
            s_ulab = sprintf('input: $u_{%s}(k)$ %s', opts.datanames{kk}, opts.uunits);
            ylabel(s_ulab, 'interpreter', 'latex');
            xlabel('t [s]', 'interpreter', 'latex');

            h2.DisplayName = opts.name;
            h1.DisplayName = opts.name;
            h3.DisplayName = opts.name;

            H_line(:,kk) = [h1;h2; h3];
            H_text = [];
            if opts.textOn
                if ~isnan(opts.TS_s(kk))
                    s1 = sprintf('%s settle time: %f', opts.name, opts.TS_s(kk));
                else
                    s1 = sprintf('%s settle time: n/a', opts.name);
                end
                ss1 = catTexts(s1, F1);
                H_text = text(.3, .1, ss1, 'Units', 'normalized', 'VerticalAlignment', 'bottom');
            end
    end
end


% =========================================================================
% =========================================================================

function plotSettleBoundary(obj, opts, F1, H_zoomSubplots)
    for kk = 1:length(H_zoomSubplots)
       
        y_upper = [obj.y_ref(kk) + obj.y_ref(kk)*obj.TOL,...
                   obj.y_ref(kk) + obj.y_ref(kk)*obj.TOL]*opts.yscaling;
        y_lower = [obj.y_ref(kk) - obj.y_ref(kk)*obj.TOL,...
                   obj.y_ref(kk) - obj.y_ref(kk)*obj.TOL]*opts.yscaling;
        subplot(H_zoomSubplots(kk))
            x_s = xlim;
            plot(x_s,y_upper, ':k')
            plot(x_s, y_lower, ':k')
            plot(x_s, [obj.y_ref(kk), obj.y_ref(kk)]*opts.yscaling, '--k')
    end
end
% =========================================================================
% =========================================================================
function ss1 = catTexts(s1, F1)

    c = findobj(F1, 'Type', 'Text');
    if ~isempty(c)
        [rowsz] = size(c(1).String);
        if rowsz >1
        % Sometimes, matlab makes this thing have rows. deal with that case
        % here.
            ss1 = [s1];
            for i = 1:length(c(1).String)
                if ~strcmp(s1, c(1).String(i))
                    keyboard
                    ss1 = [ss1, '\newline', c(1).String(i)];
                end
            end
        else % mostly, its a single row. That case here;
            kk = strfind(c(1).String, '\newline');
            try
                kk = kk(1);
            catch
                kk = length(c(1).String);
            end
            if isempty(strfind(c(1).String(1:kk), s1))
                ss1 = [s1, '\newline', c(1).String];
            else
                ss1 = c(1).String;
            end
        end
        delete(c)
    else
        ss1 = s1;
    end


end


function opts = buildStepExpPlotOpts(obj, pre_opts)
    % concatenates the opts from obj, (ie all field names except for 'y',
    % and 'u', with the extra options derived or defined  inside the plotting function'

    fnames_obj = fieldnames(obj);
    fnames_pre = fieldnames(pre_opts);
    fnames = [fnames_obj; fnames_pre];
    
    
    data_pre = struct2cell(pre_opts);
    props_obj = properties(obj);
    data_obj = cellfun(@(props_obj) obj.(props_obj), props_obj, 'UniformOutput', false);
    data = [data_obj; data_pre];

    opts = cell2struct(data, fnames);
    opts = rmfield(opts, {'u', 'y'});

    opts.yunits = opts.yunits_scaled;
        
%         y;
%         u;
%         pstyle;
%         name;
%         exptype;
%         params;
%         controller;
%         datanames;
%         yunits;
%         uunits;
%         yscaling;
%         y_ref;
%         TOL;
    
% Ensure we have sensible defaults, in case needed fields are empty.
nChan = size(obj.y,2);

defs = {'datanames', repmat({' '}, nChan, 1);
        'yunits', ' ';
        'uunits', ' ';
        'yscaling', 1;
        'TOL', 0;};

for pair=defs'
    field = pair{1};
    if isempty(opts.(field))
        opts.(field) = pair{2};
    end
end 

end



