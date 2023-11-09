clear all; close all; clc

%% Global variables
global Fg mMat dt
global kappaBar EI voronoiLength % Bending
global GJ % Twisting
global EA refLen % Stretching

%% Inputs
nv = 50; % number of nodes
ne = nv - 1; % number of edges
ndof = 3*nv + ne; % number of DOF = 4*nv - 1
dt = 0.01; % Time step size

% Rod Geometry
RodLength = 0.2;
natR = 0.02; % Natural radius (m)
r0 = 0.001; % Cross-sectional radius
rho = 1000; % Density 

% Material parameters
Y = 10e6; % Young's Modulus
nu = 0.5; % Poisson's ratio
G = Y / (2 * (1 + nu)); % Shear modulus

% Gravity
g = [0; 0; -9.81];
totalTime = 5; % Simulation time

%% Stiffness variables
EI = Y * pi * r0^4/4; % Bending stiffness
GJ = G * pi * r0^4/2; % Shearing stiffness
EA = Y * pi * r0^2; % Stretching stiffness

%% Tolerance
tol = EI / RodLength^2 * 1e-6; %characteristic bending force estimate

%% Mass
totalM = pi * r0^2 * RodLength * rho; %total mass of rod
dm = totalM / ne; % mass per edge

%create mass vector
massVector = zeros(ndof, 1);

for c = 1:nv % Loop over nodes
    ind = [4*c-3; 4*c-2; 4*c-1]; %cth node

    if c == 1 
        massVector(ind) = dm/2; %end node
    elseif c == nv
        massVector(ind) = dm/2; %end node
    else
        massVector(ind) = dm; %all internal nodes
    end
end

for c=1:ne % Loop over edges
    massVector(4*c) = 1/2 * dm * r0^2;
end
mMat = diag(massVector); % ndof x ndof sized mass matrix

%% Initial Geometry
nodes = zeros(nv, 3);
dTheta = (RodLength / natR) * (1/ne); %angle from center = arc length / radius

%Generate node positions
for c = 1 :nv
    nodes(c,1) = natR * cos( (c-1) * dTheta);
    nodes(c,2) = natR * sin( (c-1) * dTheta);
    nodes(c,3) = 0; %Pitch = 0
end

%% Initial DOF vector
q0 = zeros(ndof, 1); % ndof = 4N-1
for c=1:nv
    ind = [4*c-3; 4*c-2; 4*c-1]; %cth node
    q0(ind) = nodes(c,:);
end

%% Gravity
%Create gravity matrix before looping
Fg = zeros(ndof, 1);

for c=1:nv % loop over the nodes
    ind = [4*c-3; 4*c-2; 4*c-1];
    Fg(ind) = massVector(ind) .* g;
end

%% Reference frame (Space parallel transport at t=0)
a1 = zeros(ne, 3); % First reference director for all edges
a2 = zeros(ne, 3); % Second reference director for all edges
tangent = computeTangent( q0 ); % Tangent

t0 = tangent(1,:); % Tangent on first edge
t1 = [0;0;-1]; % "arbitrary" vector
a1Tmp = cross(t0, t1); % is perpendicular to both t0 and t1

if abs(a1Tmp) < 1e-6 %just in case a1Tmp is parallel to [0; 1; 0]
    t1 = [0;1;0]; % arbitrary
    a1Tmp = cross(t0, t1);
end

a1(1,:) = a1Tmp / norm( a1Tmp ); % Plug into my big a1 matrix
a2(1,:) = cross(tangent(1,:), a1(1,:)); %a2 is orthogonal to a1

% Done with the first edge
% Space parallel transport to construct the reference frame
for c=2:ne %loop over remaining edges
    t0 = tangent(c-1,:); % tanget on c-1 edge
    t1 = tangent(c,:); % tanget on c-th edge
    a1_0 = a1(c-1, :); %previous a1 director on c-1 edge
    a1_l = parallel_transport( a1_0, t0, t1); %parallel transport to subsequent edge
    a1(c,:) = a1_l / norm( a1_l ); %normalize a1 and a2 again
    a2(c,:) = cross( t1, a1(c,:) );
end

%% Material frame
theta = q0(4:4:end);
[m1, m2] = computeMaterialDirectors(a1, a2, theta);

%% Reference twist
refTwist = zeros(nv, 1); % Get starting reference twist

%% Natural curvature
kappaBar = getKappa(q0, m1, m2); %Get natural curvature

%% Reference length or edge length, used in stretching force
refLen = zeros(ne, 1); %create edge length vector

for c = 1:ne % loop over the edges
    dx = nodes(c+1,:) - nodes(c,:); %dx = x_(c+1) - x_c
    refLen(c) = norm( dx ); %get length
end

%% Voronoi length (Length associated with each node), used in benging and twisting
voronoiLength = zeros(nv, 1); %create voronoi length vector to get node "length"

for c=1:nv % loop over the nodes
    if c==1
        voronoiLength(c) = 0.5 * refLen(c);
    elseif c==nv
        voronoiLength(c) = 0.5 * refLen(c-1);
    else
        voronoiLength(c) = 0.5 * (refLen(c-1) + refLen(c));
    end
end

%% Fixed and free DOFs
fixedIndex = 1:7; %first two nodes are fixed 
freeIndex = 8:ndof;

%% Time stepping scheme
Nsteps = round(totalTime/dt); %define number of timesteps
ctime = 0; % current time (utility variable)
endZ = zeros(Nsteps, 1); % z-coordinate of last node

% Initialize position and velocity
q = q0; % Position
u = zeros(size(q)); % Velocity

for timeStep = 1:Nsteps
    fprintf('Current time=%f\n', ctime);
    % Guess q
    q = q0;
    iter = 1; % Number of iterations
    error = 10 * tol;

    while error > tol
        % Compute reference frame
        [a1Iterate, a2Iterate] = computeTimeParallel(a1, q0, q);
        
        % Compute reference twist
        tangent = computeTangent(q);
        refTwist_iterate = computeRefTwist(a1, tangent, refTwist);
        
        % Compute material frame
        theta = q(4:4:end);
        [m1, m2] = computeMaterialDirectors(a1Iterate, a2Iterate, theta);
        
        % Elastic Force and Elastic Jacobian calculation
        [Fb, Jb] = getFb(q, m1, m2);
        [Ft, Jt] = getFt(q, refTwist_iterate);
        [Fs, Js] = getFs(q);
        Forces = Fb + Ft + Fs + Fg; %Sum all forces and Jacobians
        JForces = Jb + Jt + Js;
        
        % Equations of motion
        f = mMat / dt * ( (q-q0) / dt - u) - Forces;
        
        % Jacobian
        J = mMat / dt^2 - JForces;
        f_free = f(freeIndex); %only iterate the free indices!
        J_free = J(freeIndex, freeIndex);
        
        % Newton's update
        dq_free = J_free \ f_free; 
        q(freeIndex) = q(freeIndex) - dq_free;
    
        % Error
        error = sum( abs( f_free ) );
        fprintf('Iter=%d, error=%f\n', iter, error);
        iter = iter + 1; %iterate time
    end
    
    %update velocity
    u = (q - q0) / dt;
    a1 = a1Iterate; %iterate the reference frame
    a2 = a2Iterate;
    ctime = ctime + dt; %iterate current time
    
    % Update q
    q0 = q;

    % Store
    endZ(timeStep) = q(end); %store the z-position of the last node
    if mod(timeStep, 100) == 0 %plot the rod every second
        theta = q(4:4:end);
        [m1, m2] = computeMaterialDirectors(a1, a2, theta); %update the material directors
        plotrod(q, a1, a2, m1, m2, ctime);
    end
end

%% Visualization
figure(2);
timearray = (1:1:Nsteps) * dt;
plot( timearray, endZ, 'ro-');
box on
xlabel('Time, t [sec]');
ylabel('z-coord of last node, \delta_z [m]');
