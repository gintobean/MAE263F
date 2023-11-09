% Inputs:
% a1 - first time parallel reference director (ne x 3)
% a2 - second time parallel reference director (ne x 3)
% theta - twist angle vector (ne x 1)

% Outputs:
% m1 - first material directors on each edge (ne x 3)
% m2 - second material directors on each edge (ne x 3)

function [m1,m2] = computeMaterialDirectors(a1,a2,theta)

    %set dimensions
    ne = length(theta);
    
    %set material directors
    m1 = zeros(ne, 3);
    m2 = zeros(ne, 3);
    
    for i = 1:ne %loop over edges
        %create orthogonal m1, m2 by rotating each by 90 degrees
        m1(i,:) = cos(theta(i)) * a1(i,:) + sin(theta(i)) * a2(i,:); 
        m2(i,:) = -sin(theta(i)) * a1(i,:) + cos(theta(i)) * a2(i,:);
    end
end

