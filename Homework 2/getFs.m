% Inputs: 
% q - DOF vector (4*nv - 1 x 1)
%
% Outputs:
% Fs - stretching force vector (4*nv - 1 x 1)
% Js - stretching Jacobian vector (4*nv - 1 x 4*nv - 1)

function [Fs, Js] = getFs(q)
    global EA refLen

    %define dimensions
    nv = (length(q)+1) / 4;
    
    %initialize Fs and Js
    Fs = zeros(size(q)); %stretching force vector 
    Js = zeros(length(q), length(q)); %jacobian of stretching force
    
    for i=1:nv-1 % loop over nodes
        n1 = transpose(q(4*i-3:4*i-1)); % current node
        n2 = transpose(q(4*i+1:4*i+3)); % next node

        %compute stretching energies using appendix functions
        [dF, dJ] = gradEs_hessEs(n1, n2, refLen(i), EA);

        ind = [4*i-3, 4*i-2, 4*i-1, 4*i+1, 4*i+2,4*i+3]; % 6 numbers, exclude the edge DOF

        %append dF and dJ to force vector and jacobian
        Fs(ind) = Fs(ind) - dF;
        Js(ind, ind) = Js(ind, ind) - dJ;
    end
end
