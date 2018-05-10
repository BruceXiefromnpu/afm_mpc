% allmatrix = packMatrix(sys_obs, L, K)
%
% Packs observer system matrices, control gain, and observer gain into a
% single column vector. To be read by labview for experiement.
% This is performs the same functionality but in addition to packing 
% [B;L;K';Arows;]
% it appends xss to the end s.t
% [B;L;K';Arows; xss]
%
% ASSUMES A IS PROVIDED AS A~ = A-LC
function AllMatrix = packMatrixDistEstMPC_singleAxis(sys_obs, L, Nx)
    A = sys_obs.a;
    B = sys_obs.b;
    % C = sys_obs.c;

    %Ni = size(B,2);
    %No = size(C,1);
    Ns = size(B, 1);

    AllMatrix = [B; 
                 L;
                 B*0]; % filler for K
    
    for i=1:Ns
       AllMatrix = [AllMatrix; A(i,:)']; 
    end
    
    AllMatrix = [AllMatrix;
                 Nx];

end