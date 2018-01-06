% [traj_s, settle_times] = opt_traj_gen(Q, R, N_traj, sys, ref_s, ref_0_s, varargin )
% The goal of this function is to generate a set of optimal trajectories
% over a LONG horizon (for a sequence of setpoints). The immediate goal
% here is to look at what is the minimum acheivable settling time for a set
% of setpoints, given a (Q,R) pair, laying aside all issues about the MPC
% control horizon.
% 
% Required Inputs
% ------
%  Q, R: state and control weighting matrices
%  
%  N_traj : the TOTAL length of the optimal trajectory. Eg., 400 is
%  usuall reasonable for the AFM
% 
%  sys : discrete time dynamical system. The optimal control problem will
%  be generated based on sys.
%
%  ref_s : a list of reference inputs to generate trajectories for.
%  
%  ref_0_s : a list of starting points. We will assume the system is at
%  steady for y = ref_0(iter). Can be either a single point, in which case
%  the same ref_0 is used for all trajectories, or can be a list the same
%  length as ref_s.
%
%  Optional Inputs
%  ---------------
%    opt_traj_gen(..., 'S', S) : cross weighting matrix for the optimal
%    control problem
%
%    opt_traj_gen(..., verbose, 1) (default 0) flag indiciating what plots to make
%        verbose = 1, plot output
%        verbose = 2, plot output and control
%        verbose = 3, plot output, control and accumulated control (useful
%        if doing optimization over incremental form)
%
%   opt_traj_gen(..., 'uMax', umax) specifies that a box constraint will be
%   added to the quadratic program of [-umax, umax]. If optimzing with the
%   incremental form, use (..., 'uMax', du_max)
%
%   opt_traj_gen(..., 'slewMax', slewmax) specifies that a slew-rate input
%   constraint will be added to the problem.



function [traj_s, settle_times] = opt_traj_gen(Q, R, N_traj, sys, ref_s, ref_0_s, varargin )



ns = size(sys.b, 1);
nu = size(sys.b, 2);
no = size(sys.c, 1);

if length(ref_0_s) == 1
    ref_0_s = 0*ref_s + ref_0_s;
elseif length(ref_0_s) ~= length(ref_s)
    error('length(ref_0_s) must be either one or the same as length(ref_s)')
end

p = inputParser;
defaultS = zeros(size(sys.b, 1), size(sys.b, 2));
defaultUMax = 0;
defaultSlewMax = 0;
addParameter(p, 'S', defaultS);
addParameter(p, 'verbose', 0);
addParameter(p, 'uMax', defaultUMax);
addParameter(p, 'slewMax',defaultSlewMax);

parse(p, varargin{:})
S = p.Results.S;
uMax = p.Results.uMax;
slewMax = p.Results.slewMax;
verbose = p.Results.verbose;

Qp = dare(sys.a, sys.b, Q, R);
NLQR_prob = sparseMPCprob(sys, N_traj, Q, Qp, R, S);

if uMax ~= 0
    NLQR_prob.add_U_constraint('box', [-uMax, uMax]);
end
if slewMax ~= 0 
    NLQR_prob.add_U_constraint('slew', slewMax);
end

if verbose > 0
%     mkfig(verbose)
    Fig = figure();
else
    Fig = 0;
end

[Nx, Nu] = SSTools.getNxNu(sys);

N_refs = length(ref_s);
Y_vec_s = cell(1, N_refs);
U_vec_s = cell(1, N_refs);
X_vec_s = cell(1, N_refs);

tvec = [0:1:N_traj-1]*sys.Ts;
settle_times = zeros(1, length(ref_s));
start_str = sprintf('gamma: %.0f', R);
upd = progressbar(length(ref_s), 'start_str', start_str);
upd(0);
for iter = 1:length(ref_s)
    ref_f = ref_s(iter);
    ref_0 = ref_0_s(iter);
    
    x0_err = ref_0*Nx - ref_f*Nx;
    
    [U, Xerr] = NLQR_prob.solve(x0_err); 
    
    % Convert the solution trajectory from error coords to regular coords.
    X = Xerr + ref_f*Nx;
    Y = sys.c*X;
    t_settle = settle_time(tvec, Y, ref_f, 0.01*ref_f,...
                                  [], [], 30);
    settle_times(iter) = t_settle;
    % We have to transpose to columns, otherwise matlab will create a 3d
    % matrix out of this stuff, for some bizare reason.
    Uvec = timeseries(U', tvec);
    
    % We have one extra element N+1 states, but only N controls. So drop
    % the last one. 
    Xvec = timeseries(X(:,1:end-1)', tvec);
    Yvec = timeseries(Y(:,1:end-1)', tvec);
    
    U_vec_s{iter} = Uvec;
    X_vec_s{iter} = Xvec;
    Y_vec_s{iter} = Yvec;
    
    if verbose > 0
        plot_local(Xvec, Yvec, Uvec, verbose, Fig)
    end
    
    upd(iter);
    
end
% fprintf('PUT\n')
traj_s.X_vec_s = X_vec_s;
traj_s.Y_vec_s = Y_vec_s;
traj_s.U_vec_s = U_vec_s;

end

function plot_local(Xvec, Yvec, Uvec, verbose, Fig)

if verbose == 1
   subplot(111); hold on
   plot(Yvec.time, Yvec.data)
   ylabel('y(k)')
   xlabel('time')
   drawnow()
   
elseif verbose ==2
    subplot(211); hold on
    plot(Yvec.time, Yvec.data)
    ylabel('y(k)')
    
    subplot(212); hold on
    plot(Uvec.time, Uvec.data)
    ylabel('u(k)')
    xlabel('time')
    drawnow()
elseif verbose == 3
    subplot(311); hold on
    plot(Yvec.time, Yvec.data)
    ylabel('y(k)')
    
    subplot(312); hold on
    plot(Uvec.time, Uvec.data)
    ylabel('u(k)')
    
    subplot(313); hold on
    plot(Uvec.time, cumsum(Uvec.data))
    ylabel('accum(u(k))')
    xlabel('time')
    drawnow()    
end



end




