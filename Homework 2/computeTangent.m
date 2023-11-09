% Inputs:
% q - DOF vector (4*nv - 1 x 1)

% Outputs: 
% tangent vector for each edge (ne x 3)

function tangent = computeTangent(q)
    %define dimensions
    nv = (length(q) + 1) / 4;
    ne  = nv - 1;
    
    %initialize tangent vector
    tangent = zeros(ne,3);
    
    for i = 1:ne %loop over edges
        x0 = q(4*i-3:4*i-1); %position of current node
        x1 = q(4*i+1:4*i+3); %position of next node
        edge = x1 - x0; %vector representing the edge
        tangent(i,:) = edge / norm(edge); %insert normalized edge into edge vector
    end
end