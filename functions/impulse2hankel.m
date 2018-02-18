% Hk = getHankelMIMO(y, k, Ni, No)
%  This function will generate an No*(r-1) X Ni*(s-1) Hankel matrix from
%  impulse response data. 
%
%  Inputs:
%  y = Impulse response data. We assume mimo data. Y has Ni rows, No cols, 
%  and length(t) pages.
%  k = hankel index: H(k)
%  r,s are the number of (block) rows and columns respectively.
%  
%   Outputs:
%
%  Hrs(k-1) = [Y(k), ...., Y(k+s-1)
%            | :,             :
%            | Y(k+r-1), .... Y(k+r+s-1)
%
%  
%
%   This definition is taken from the paper:
%   "Multivariable Model Identification from Frequency Response Data", by
%   Robert N. Jacques and David W. Miller, CDC, 1993.
% 
% Arnold Braker
% 5-7-2016

function Hk = getHankelMIMO(Y, k, r, s, ni, no)

    % SISO case (Old)
    % for i =0:r-1
    %    for j = 0:s-1
    %       Hk(i+1,j+1) = y(k+i+j); 
    %    end    
    % end
    if length(size(Y))<3% Ensure data is stored as *pages*. Catch the old 
                             % Where it's stored as rows. 
        error('Impulse data does not appear to be store as pages. Please re-arrange it.')
    end

                             
    % In MIMO, data for different times stored in different pages, row col is
    % different inputs/outputs.
    Hk = zeros(r*no, s*ni);


    for i=0:r-1
       for j =0:s-1 

            row_i = i*no+1; % i=0, get 1, i=1, get 3, for no= 2
            col_j = j*ni+1; % j=0, get 1, j=1, get 3, for ni= 2
            Hk(row_i:row_i+no-1, col_j:col_j+ni-1) = Y(:,:,k+i+j);
       end
    end

end