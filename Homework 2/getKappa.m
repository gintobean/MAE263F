% Inputs: 
% q  - DOF Vector (4*nv - 1 x 1)
% m1 - 1st material director (ne x 3)
% m2 (2nd matrerial director (ne x 3)

% Outputs:
% The value of kappa, which is the discrete curvature of the rod (nv x 2)

function kappa = getKappa(q,m1,m2)

    %define dimensions
    nv = (length(q) + 1) / 4;
    ne = nv - 1;

    %initialize kappa vector
    kappa = zeros(nv,2);

    for i = 2:ne %loop over all nodes except first and last node
        n0 = q(4*i-7:4*i-5); %node prior to turning node
        n1 = q(4*i-3:4*i-1); %turning node
        n2 = q(4*i+1:4*i+3); %node after turning node

        m1e = m1(i-1,:); %material director 1 of prior edge
        m2e = m2(i-1,:); %material director 2 of prior edge
        m1f = m1(i,:); %material director 1 of following edge
        m2f = m2(i,:); %material director 2 of following edge

        %invoke computekappa() to get kappa values
        currentKappa = computekappa(n0,n1,n2,m1e,m2e,m1f,m2f);
        
        %transfer kappa values to curvature array
        kappa(i,1) = currentKappa(1);
        kappa(i,2) = currentKappa(2);
    end
end