% Inputs: 
% q - DOF vector (4*nv - 1 x 1)
% refTwist - twisting  (ne x 3)
%
% Outputs:
% Ft - twisting force vector (4*nv - 1 x 1)
% Jt - twisting Jacobian vector (4*nv - 1 x 4*nv - 1)

function [Ft, Jt] = getFt(q, refTwist)
    global GJ voronoiLength

    %define dimensions
    nv = (length(q)+1) / 4;

    %initialize Ft and Jt
    Ft = zeros(size(q)); %twisting force vector
    Jt = zeros(length(q), length(q)); %jacobian of twisting energy
    
    for i=2:nv-1 % Compute twisting force at each node except first and last nodes
        n0 = transpose(q(4*i-7:4*i-5)); % previous node
        n1 = transpose(q(4*i-3:4*i-1)); % current node
        n2 = transpose(q(4*i+1:4*i+3)); % next node
        te = q(4*i-4); %twisting angles of previous edge
        tf = q(4*i); %twisting angles of current edge

        %compute twisting energies using appendix functions
        [dF, dJ] = gradEt_hessEt(n0, n1, n2, te, tf, refTwist(i), voronoiLength(i), GJ);
        
        ind = 4*i-7:4*i+3; % 11 numbers

        %append dF and dJ to force vector and jacobian
        Ft(ind) = Ft(ind) - dF;
        Jt(ind, ind) = Jt(ind, ind) - dJ;
    end
end
