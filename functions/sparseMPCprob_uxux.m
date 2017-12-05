classdef sparseMPCprob
%         H;
%         M;
%         Ainq;
%         binq;
%         N_mpc;
%         kappa;
%
%
% Construction: condensedMPCprob(sys,N, Q,Qp, R)
%               condensedMPCprob(sys,N, Q,Qp, R, S)
%
% This sets H, M, N_mpc, and kappa fields. We leave it to you to fill in
% Ainq and binq.


    properties
        H;
        Aeq;
        beq;
        Ainq;
        binq;
        N_mpc;
        kappa;
        ns;
        nu;
        
    end
    
    methods
        function obj = sparseMPCprob(sys,N, Q,Qp, R, S)
            if ~exist('S', 'var')
                ns = size(sys.b,1);
                nu = size(sys.b, 2);
                S = zeros(ns, nu);
            end
            [H, Aeq, beq] = clqrProblem_local(sys,N, Q, R, Qp, S);
            obj.H     = H;
            obj.Aeq   = Aeq;
            obj.beq   = beq;
            obj.N_mpc = N;
            obj.kappa = cond(H);
            obj.ns = size(sys.b,1);
            obj.nu = size(sys.b,2);
        end
        
        function UX = solve(obj, xk_1)
            
            beq_xk = obj.beq;
            beq_xk(1:obj.ns) = xk_1;
            UX = mpcSolve_local(obj.H, obj.Aeq, beq_xk);
            
            upad = zeros(obj.nu,1)+NaN;
            UX = reshape([UX;upad], obj.ns+obj.nu,[]);
            
        end
        
        
    end
    
end

function UX = mpcSolve_local(H, Aeq, beq)

    f  = zeros(size(H,2),1);
    UX = quadprog(H, f, [], [], Aeq, beq);
end


function [H, Aeq, beq] = clqrProblem_local(sys, N, Q, R, Qp, S)
    if sys.Ts == 0
        error('system "sys" should be discrete time dynamical system')
    end
    
    
    Ns = size(sys.b, 1);
    Nu = size(sys.b, 2);

    [a, b] = ssdata(sys);

    % -------------------------------------------------------------------------
    % Sparse cLQR problem

    % Build the cost function. THis should look like:
    %      [Q  S        ]
    %      [S' R        ]
    %  H = [     Q  S   ]
    %      [     S' R   ]
    %      [          Qp]
    
    QR = [Q, S; 
          S', R];
    bigI = eye(N-1);
    H = kron(bigI, QR);
    H = blkdiag(H, Qp);
    
    % Now build the equality constraint. This should look like:
    %  [I 0 0  0 0 0 ][x0]   [x0]
    %  [A B -I 0 0 0 ][u0] = [0 ]
    %  [0 0  A B -I  ][x1]   [0 ]
    %                 [u1]
    %                 [x_N] 
    AB_I = [a, b, -eye(Ns)];
    
    Aeq = zeros(Ns*N, (Ns+Nu)*N-Nu);
    beq = zeros(Ns*N, 1);

    colStart = 1;
    for iter = 1:Ns:Ns*N
        if iter==1
            Aeq(1:Ns,1:Ns) = eye(Ns);
        else
            Aeq(iter:Ns+iter-1, colStart:(colStart+2*Ns+Nu-1)) = AB_I;
            colStart = colStart+Ns+Nu;
        end
    end
end

