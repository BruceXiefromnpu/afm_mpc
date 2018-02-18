% [Ad, Bd, Cd] = impulse2ss(yy_exp, r, s)
%
% Performs the eigensystem realization algorithm from the paper,
%   "Multivariable Model Identification from Frequency Response Data", by
%   Robert N. Jacques and David W. Miller, CDC, 1993.
%
% Inputs:
%   YY_imp: the impulse response data. This should be a 3d array with No
%   rows, Ni columns and length(t) pages.
%
%   r, s are integers, of rows and columns of hankel. Or something. Go look
% this up.




function [Ad, Bd, Cd] = impulse2ssMIMO(yy_IMP, r, s, Ns)

no = size(yy_IMP, 1);
ni = size(yy_IMP, 2);
fprintf('Assuming number if outputs = %d, number inputs = %d\n', no, ni) 
% Form Hankel matrices and start the Eigen system Realization
k = 1;
H0 = impulse2hankel(yy_IMP, k, r, s,ni, no);
k = 2;
Hmag1 = impulse2hankel(yy_IMP, k, r, s, ni,no);

fprintf('size of Hankel is %d x %d \n', size(Hmag1,1), size(Hmag1,2));
[P, D, Q] = svd(H0, 'econ');
% % H = P*D*Q'
% keyboard
I_no = eye(no);
I_ni = eye(ni);
Z_no = zeros(no,no);
Z_ni = zeros(ni, ni);


Ep = repmat(Z_no, r,1);
Ep(1:no,1:no) = I_no;

Em = repmat(Z_ni, s, 1);
Em(1:ni,1:ni) = I_ni;

% SISO case:
% Ep    = zeros(r,1);
% Ep(1) = 1;
% Em    = zeros(r,1);
% Em(1) = 1;

if exist('Ns', 'var')
   Dp = zeros(size(D));
   D_neg_sqrt = Dp;
   D_sqrt     = Dp;
   
   D_neg_sqrt(1:Ns, 1:Ns) = (D(1:Ns, 1:Ns))^-0.5;
   D_sqrt(1:Ns,1:Ns)      = sqrt(D(1:Ns, 1:Ns));
else
    D_neg_sqrt = D^-0.5;
    D_sqrt     = D^0.5;
end


try
    Ad = D_neg_sqrt*P'*Hmag1*Q*D_neg_sqrt;
catch
    keyboard
end
try
    Bd = D_sqrt*Q'*Em;
catch
    keyboard
end
Cd = Ep'*P*D_sqrt;


if exist('Ns', 'var')
    Ad = Ad(1:Ns,1:Ns);
    Bd = Bd(1:Ns,:);
    Cd = Cd(:,1:Ns);
end
% keyboard


end

