% [ Chat, Dhat, Q, S, Rmin ] = fict_zeros(sys, zdes)
%
% Returns Fictitious output matrices Chat and Dhat such that the zeros of
% the system
% |sys.A  sys.B|
% |Chat   Dhat | 
% 
% are at the locations specified in zdes. 
%
% This is useful for pole placment in LQR. To effect such a design, this
% function also returns
%
% Q = Chat^T * Chat
% S = Chat^T
% Rmin = Dhat^2
%
% In general, we can choose Dhat := 1 and then solve for Chat.
%
% In LQR, this corresponds to 
%
% min_{U} sum (1/2) ( yhat^That + u*R*u )
%         = (1/2) ( (Chat*x + Dhat*u)^T(Chat*x + Dhat*u) + u*R*u )
%         = (1/2)x^T*Chat^T*Chat)*x + x^T*Chat^T*Dhat*u + u^T(Dhat^T*Dhat + R)u
%
% Thus, we set Q = Chat^T*Chat  
%              S = Chat^T*Dhat = Chat^T
%              Rmin = 1


function [ Chat, Dhat, Q, S, Rmin ] = fict_zeros(sys, zdes)


Chat = place(sys.A, sys.B, zdes);
Dhat = 1;
Q = Chat'*Chat;
S = Chat';

Rmin = Dhat^2;


end

