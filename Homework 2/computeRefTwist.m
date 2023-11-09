% Inputs:
% a1 - 1st time-parallel reference director on each edge (ne x 3)
% tangent - tangent vectors on each edge (ne x 3)
% refTwist - twist for each node (nv x 1)

% Outputs:
% Computes vector with reference twist of each node (nv x 1)

function refTwist = computeRefTwist(a1, tangent, refTwist)
    
    %define dimensions
    ne = size(a1); 
    nv = ne + 1;

    for i=2:ne %loop over edges
        u0 = a1(i-1,:); %a1 of prior edge
        u1 = a1(i,:); %a1 of current edge
        t0 = tangent(i-1,:); %tangent of prior edge
        t1 = tangent(i,:); %tangent of current edge
        ut = parallel_transport(u0, t0, t1); %space parallel vector 
    
        %compute reference twist using signedAngle()
        refTwist(i) = signedAngle(ut,u1,t1);
    end
end