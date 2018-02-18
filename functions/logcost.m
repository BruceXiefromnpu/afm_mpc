
function J = logcost(theta, P, w)

    w_weight = [w(1);w;w(end)];
    
   logW = log(w_weight(3:end)) - log(w_weight(1:end-2));
%     logW = log(W);

   
    H_frf = H(w, theta);
    J_vec = abs(log(H_frf./P)).^2;
%     J_vec = J_vec.*logW;
%     keyboard
%     J_vec = J_vec
    J = sum(J_vec);



end


function H_frf = H(w, theta)
Ts = 40e-6;
nz = 9;
nnum = 10;

g = tf(theta(1:nnum)', [1 theta(nnum+1:end)'], Ts);
H_frf = squeeze(freqresp(g, w));

end