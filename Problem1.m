%Filename: Problem1.m
%Simulation of Spheres and Beam in Viscous Flow: 
%Explicit and Implicit Methods (3 Node Case)

%Simulation Size
N = 3; %number of nodes
ndof = N * 2; %2-D problem

%Physical Values
rhoMetal = 7000; %metal density kg/m^3
rhoFluid = 1000; %fluid density kg/m^3
viscosity = 1000; %viscosity Pa-s
g = 9.8; %gravity acceleration m/s^2

%Radii of Spheres m
R1 = 0.005; 
R2 = 0.025;
R3 = 0.005;

%Beam Parameters
l = 0.1; %length of beam m
deltaL = l / (N-1); %distance between each beam m
r0 = 0.001; %radius of beam cross-section m
E = 1e9; %elastic modulus of beam Pa
A = pi() * r0^2; %cross sectional area of beam m^2
I = pi() * r0^4 / 4; %moment of inertia beam m^4

%Simulation Parameters
dt = 1e-2; %time step
eps = E*I/l^2 * 1e-3; %tolerance based on small force
maxTime = 10;

%Mass Matrix
M = zeros(ndof, ndof);
M(1,1) = 4/3*pi*R1^3*rhoMetal; %mass of spheres kg
M(2,2) = 4/3*pi*R1^3*rhoMetal;
M(3,3) = 4/3*pi*R2^3*rhoMetal;
M(4,4) = 4/3*pi*R2^3*rhoMetal;
M(5,5) = 4/3*pi*R3^3*rhoMetal;
M(6,6) = 4/3*pi*R3^3*rhoMetal;

%Weight matrix
W = zeros(ndof,1);
W(2) = 4/3*pi*R1^3*(rhoMetal-rhoFluid)*g; %weight of spheres kg
W(4) = 4/3*pi*R2^3*(rhoMetal-rhoFluid)*g;
W(6) = 4/3*pi*R3^3*(rhoMetal-rhoFluid)*g;


%Damping Matrix
C = zeros(6,6);
C(1,1) = 6*pi*viscosity*R1; %damping coefficient
C(2,2) = 6*pi*viscosity*R1;
C(3,3) = 6*pi*viscosity*R2;
C(4,4) = 6*pi*viscosity*R2;
C(5,5) = 6*pi*viscosity*R3;
C(6,6) = 6*pi*viscosity*R3;

%%Implicit Method -------------------------------------------

%Initial DOF Vector 
q0 = zeros(ndof,1); 
for c=1:N
    q0(2*c-1) = (c-1) * deltaL; %initial x configuration
    q0(2*c) = 0; %initial y configuration
end

%Initial Velocity Vector
u0 = zeros(ndof,1); 

%Number of Steps
steps = round(maxTime / dt); 

%Middle-node Position and Velocity Arrays
midNodePosition = zeros(steps,1);
midNodeVelocity = zeros(steps,1); 

%Time-Stepping Scheme
for k = 1:steps-1
    q = q0; %update old DOFS
    u = u0; %update old velocities
    err = 10 * eps; %starting error

    %Newton-Rhapson Iteration
    while err > eps
        %Initialize Function and Jacobian
        f = M / dt * ((q-q0)/dt - u) + W + C * (q-q0)/dt; 
        J = M / dt^2 + C / dt;

        %Linear Spring Between Nodes 1,2
        dF = gradEs(q(1),q(2),q(3),q(4),deltaL,E*A);
        dJ = hessEs(q(1),q(2),q(3),q(4),deltaL,E*A);
        f(1:4) = f(1:4) + dF; 
        J(1:4,1:4) = J(1:4,1:4) + dJ;

        %Linear Spring Between Nodes 2,3
        dF = gradEs(q(3),q(4),q(5),q(6),deltaL,E*A);
        dJ = hessEs(q(3),q(4),q(5),q(6),deltaL,E*A);
        f(3:6) = f(3:6) + dF;
        J(3:6,3:6) = J(3:6,3:6) + dJ;

        %Bending Spring, Centered at Node 2
        curvature = 0;
        dF = gradEb(q(1),q(2),q(3),q(4),q(5),q(6),curvature,deltaL,E*I);
        dJ = hessEb(q(1),q(2),q(3),q(4),q(5),q(6),curvature,deltaL,E*I);
        f(1:6) = f(1:6) + dF;
        J(1:6,1:6) = J(1:6,1:6) + dJ;

        %Update DOF Vector and Error
        q = q - J\f;
        err = sum(abs(f));
    end
    
    %Store new positions and velocities
    u = (q-q0)/dt;
    q0 = q;
    u0 = u;

    %Store Mid-Node Position and Velocity (Node 2)
    midNodePosition(k+1) = q(4);
    midNodeVelocity(k+1) = u(4);


    %Plot beam positions (real time)
    figure(1)
    if (k == 1|| k == 2 || k == 6 || k == 11 || k == 101 || k == 999) 
    plot(q(1:2:end), q(2:2:end), 'ro-');
    axis equal;
    xlabel('x (meter)');
    ylabel('y (meter)');
    drawnow
    hold on
    end
end

%Plot mid-node y-position
t = (0:steps-1) * dt; %time array
figure(2)
plot(t,midNodePosition, 'k-'); 
xlabel('t [s]'); 
ylabel('y [m]');


%Plot mid-node y-velocity
t = (0:steps-1) * dt; %time array
figure(3)
plot(t,midNodeVelocity, 'k-'); 
xlabel('t [s]'); 
ylabel('u_y [m/s]');
ylim([-0.01 0]);

%%Explicit Method ------------------------------------------------

dt = 1e-5; %decrease dt for explicit method to converge

%Initial DOF Vector 
q0 = zeros(ndof,1); %make the q vector
for c=1:N
    q0(2*c-1) = (c-1) * deltaL; %initial x configuration
    q0(2*c) = 0; %initial y configuration
end

%Initial Velocity Vector
u0 = zeros(ndof,1); %velocity vector

%Time-Stepping Parameters
steps = round(maxTime / dt); %number of steps

%Middle-node Position and Velocity Arrays
midNodePosition = zeros(steps,1);
midNodeVelocity = zeros(steps,1); %middle-node velocity array

%Time-Stepping Scheme
for k = 1:steps-1
    q = q0; %update old DOFS
    u = u0; %update old velocities
    
    %Gradient of Elastic Energy
    dF1 = gradEs(q(1),q(2),q(3),q(4),deltaL,E*A);
    dF2 = gradEs(q(3),q(4),q(5),q(6),deltaL,E*A);
    dB = gradEb(q(1),q(2),q(3),q(4),q(5),q(6),0,deltaL,E*I);
    
    dEdq = dB;
    dEdq(1:4) = dEdq(1:4) + dF1;
    dEdq(3:6) = dEdq(3:6) + dF2;
    
    %Update DOF Vector
    q = q + dt * (u - dt * (M \ (W + C*u + dEdq)));

    %store new positions and velocities
    u = (q-q0)/dt;
    q0 = q;
    u0 = u;

    %Mid-Node Velocity and Position (Node 2)
    midNodePosition(k+1) = q(4);
    midNodeVelocity(k+1) = u(4);

end

%plot beam position (not in real time, since explicit is slow)
figure(4)
clf;
plot(q(1:2:end), q(2:2:end), 'ro-');
axis equal;
xlabel('x (meter)');
ylabel('y (meter)');

%plot mid-node y-velocity
t = (0:steps-1) * dt; %time array
figure(5)
plot(t,midNodePosition, 'k-'); 
xlabel('t [s]'); 
ylabel('u_y [m/s]');

%plot mid-node y-velocity
t = (0:steps-1) * dt; %time array
figure(6)
plot(t,midNodeVelocity, 'k-'); 
xlabel('t [s]'); 
ylabel('u_y [m/s]');
ylim([-0.01 0]);
