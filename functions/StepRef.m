classdef StepRef
  %STEPREF Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    n_space;        % 
    step_idx;       % index of each step
    step_diff_amps; % differential amplitudes
    step_amps;      % absolute amplitudes
    yref;
    yscaling;
    
  end
  
  methods
    function self = StepRef(step_amps, n_space, varargin)
      % self = StepRef(step_amps, n_space);
      % step_amps is a vector of amplitueds.
      % n_space is an integer, number of samples for each ref. Should be
      % scalar. 
      p = inputParser;
      p.addParameter('yscaling', 1);
      p.parse(varargin{:});
      
      self.yscaling = p.Results.yscaling;
      step_amps = step_amps(:);
      Ts = StageParams.Ts;
      
      self.n_space = n_space;
      
      step_amps = [0;step_amps];
      self.step_diff_amps = [0; diff(step_amps)];
      self.step_amps = step_amps;
      
      
      N_imp = length(self.step_amps)-1;
      % 1 is for default 0-ref of 1 sample.
      self.step_idx= [1; (2:n_space:N_imp*n_space)']; 
      yref_vec = zeros((N_imp+1)*(n_space)+1, 1);

      yref_vec(self.step_idx) = self.step_diff_amps;
      yref_vec = cumsum(yref_vec);
      
      t_vec = (0:length(yref_vec)-1)'*Ts;
      
      self.yref = timeseries(yref_vec, t_vec);
      
    end
    
    function [Fig, h] =  plot(self, fig_ax, varargin)
    % plot(self, fig_ax, varargin)
%       if ~exist('Fig', 'var') || ~isvalid(Fig)
%         Fig = figure();
%       end
      if ~exist('fig_ax', 'var') || ~isvalid(fig_ax)
        Fig = figure;
        ax = gca();
      elseif isa(fig_ax, 'matlab.ui.Figure')
        Fig = figure(fig_ax);
        ax = gca();
      elseif isa(fig_ax, 'matlab.graphics.axis.Axes')
        Fig = fig_ax.Parent();
        ax = fig_ax;
      else
        Fig = figure();
        ax = gca();
        varargin = {fig_ax, varargin{:}};
      end
      
      figure(Fig); 
      h = plot(ax, self.yref.Time, self.yref.Data*self.yscaling, varargin{:});
      hold(ax, 'on'); 
      grid(ax, 'on');
      
    end
    
    function [Fig] = plot_settle_boundary(self, fig_ax, TOL, tol_mode)
    % [Fig, h_line] = plot_settle_boundary(self, fig_ax, TOL, tol_mode)
      if ~strcmp(tol_mode, 'abs') && ~strcmp(tol_mode, 'rel')
        error('tol_mode must be one of "abs" or "rel"\n');
      end
      
      

      ref_s = self.step_amps;
      
      if ~exist('fig_ax', 'var') || ~isvalid(fig_ax)
        Fig = figure;
        ax = gca();
      elseif isa(fig_ax, 'matlab.ui.Figure')
        Fig = figure(fig_ax);
        ax = gca();
      elseif isa(fig_ax, 'matlab.graphics.axis.Axes')
        Fig = fig_ax.Parent();
        ax = fig_ax;
      else
        Fig = figure();
        ax = gca();
        varargin = {fig_ax, varargin{:}};
      end
      %figure(Fig);
      
      
      r = self.yref;
      for k=2:length(ref_s)
        idx_start = self.step_idx(k);
        if k == length(ref_s)
          idx_end = length(r.Time);
        else
          idx_end = self.step_idx(k+1)-1;
        end
        
        ref = ref_s(k);
        
        t0 = r.Time(idx_start);
        t1 = r.Time(idx_end);
        
        if strcmp(tol_mode, 'abs')
          TOL_k = TOL;
        else
          TOL_k = TOL*(ref - ref_s(k-1));
        end
        
        plot(ax, [t0, t1], [ref+TOL_k, ref+TOL_k]*self.yscaling, ':k');
        plot(ax, [t0, t1], [ref-TOL_k, ref-TOL_k]*self.yscaling, ':k');
        
      end

    end
    
  end % methods
    

  
end

