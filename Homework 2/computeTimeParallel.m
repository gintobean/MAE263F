% Inputs:
% a1_old - matrix containing time-reference directors at each node, from previous time step (ne x 3)
% q0 - DOF vector at previous time step (4*nv - 1 x 1)
% q - DOF vector at new time step (4*nv - 1 x 1)

% Outputs:
% a1 - first time-parallel reference director (a1) at new time step (ne x 3)
% a2 - second time-parallel reference director (a2) at new time step (ne x 3)

function [a1, a2] = computeTimeParallel(a1_old, q0, q)

    %define dimensions
    nv = (length(q)+1)/4;
    ne = nv - 1;
    
    %compute tangent vector at each edge
    tangent0 = computeTangent(q0); % Tangents at old step
    tangent = computeTangent(q); % Tangents at new step
    
    %initialize the time-parallel reference directors
    a1 = zeros(ne, 3);
    a2 = zeros(ne, 3);
    
    for i=1:ne % loop over edges
        t0 = tangent0(i,:); % tangent of ith edge at prior step
        t = tangent(i,:); % tangent of ith edge at current step
        tempOldA1 = a1_old(i,:); % old a1 director on ith edge
        tempA1 = parallel_transport( tempOldA1, t0, t); %new a1 director
  
        % Ensure that a1 and t are orthogonal
        tempA1 = tempA1 - dot(tempA1, t) * t; %graham schmidt
        tempA1 = tempA1 / norm(tempA1); %normalize
        
        %store results in a1 and a2 vectors
        a1(i,:) = tempA1; 
        a2(i,:) = cross(t, tempA1);
    end
end
