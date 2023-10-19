%Filename: Problem2.m
%Simulation of Spheres and Beam in Viscous Flow
%Explicit and Implicit Methods (N Node Case)

%----------------------- N loop -------------------------------
%Use this for loop to see relationship between N and terminal velocity
%Comment out if already doing the dt loop.
maxN = 33; %maximum N

%Terminal Velocity Array
terminalVelocity = zeros(maxN-3+1,1); 

%N Array
NArray = zeros(maxN-3 + 1,1); %Array for N
for N = 3:maxN
    disp(N); %status check if you're impatient

%----------------------- dt loop --------------------------------
% %Use this for loop to see relationship between dt and terminal velocity
% %Comment out if already doing the N loop.
% maxDt = -3; %maxmimum power for dt
% 
% %Terminal Velocity Array
% terminalVelocity = zeros(1 - maxDt,1); 
% 
% %dt Array
% dtArray = zeros(1 - maxDt, 1);
% 
% for a = 1:1-maxDt
%     dt = 10^(1-a); %dt's power is 1 - a
%     disp(1-a); %status check if you're impatient
%----------------------------------------------------------------

%Simulation Size
%N = 21; %number of nodes (comment out if looping N)
ndof = N * 2; %2-D problem

%Physical Parameters
rhoMetal = 7000; %metal density kg/m3
rhoFluid = 1000; %fluid density kg/m3
viscosity = 1000; %viscosity Pa-s
g = 9.8; %gravity acceleration m/s2

%Beam Parameters
l = 0.1; %length of beam m
deltaL = l / (N-1);
r0 = 0.001; %radius of beam cross-section m
E = 1e9; %elastic modulus of beam Pa
A = pi() * r0^2; %cross sectional area of beam m2
I = pi() * r0^4 / 4; %moment of inertia beam m4

%Simulation Parameters
dt = 1e-2; %time step (comment out if looping dt)
eps = E*I/l^2 * 1e-3; %tolerance based on small force
maxTime = 10;
midNode = round((N + 1) / 2);

%Radius Vector
R = zeros(N,1);
R(:) = deltaL / 10; %radius of spheres in m
R(midNode) = 0.025; 

%Geometry and initial configuration
nodes = zeros(N,2);
for c=1:N % %loop over nodes
    nodes(c,1) = (c-1) * deltaL; %initial x config
    nodes(c,2) = 0; %initial y config
end

%Mass Matrix
M = zeros(ndof, ndof);
for i = 1:N
    M(2*i-1,2*i-1) = 4/3*pi*R(i)^3*rhoMetal; %mass of spheres kg
    M(2*i,2*i) = 4/3*pi*R(i)^3*rhoMetal;
end

%Weight matrix
W = zeros(ndof,1);
for i = 1:N
    W(2*i-1) = 0;%no x-component in weight
    W(2*i) = 4/3*pi*R(i)^3*(rhoMetal-rhoFluid)*g; %weight of spheres kg
end

%Damping Matrix
C = zeros(ndof,ndof);
for i = 1:N
    C(2*i-1,2*i-1) = 6*pi*viscosity*R(i);
    C(2*i,2*i) = 6*pi*viscosity*R(i);
end

%%Implicit Method--------------------------------------------------------

%Initial DOF vector
q0 = zeros(ndof,1);
for c=1:N
    q0(2*c-1) = nodes(c,1); %initial x configuration
    q0(2*c) = nodes(c,2); %initial y configuration
end

%Initial Velocity Vector
u0 = zeros(ndof,1); %velocity vector

%Number of Steps
steps = round(maxTime / dt);

%Middle-node Velocity and Position Array
midNodePosition = zeros(steps,1);
midNodeVelocity = zeros(steps,1); 

%Time Stepping Scheme
for k = 1:steps-1
    q = q0; %update old DOFS
    u = u0; %update old velocities
    err = 10 * eps; %starting error

    %Newton-Rhapson Iteration
    while err > eps
        %Initialize function and Jacobian
        f = M / dt * ((q-q0)/dt - u) + W + C * (q-q0)/dt; 
        J = M / dt^2 + C / dt;
        
        %Elastic Energy of Linear Springs
        for i = 1:N-1
            dF = gradEs(q(2*i-1),q(2*i),q(2*i+1),q(2*i+2),deltaL,E*A);
            dJ = hessEs(q(2*i-1),q(2*i),q(2*i+1),q(2*i+2),deltaL,E*A);
            f(2*i-1:2*i+2) = f(2*i-1:2*i+2) + dF;
            J(2*i-1:2*i+2,2*i-1:2*i+2) = J(2*i-1:2*i+2,2*i-1:2*i+2) + dJ;
        end

        %Elastic Energy of Bending Springs
        for i = 2:N-1
            curvature = 0;
            dF = gradEb(q(2*i-3),q(2*i-2),q(2*i-1),q(2*i),q(2*i+1),q(2*i+2),curvature,deltaL,E*I);
            dJ = hessEb(q(2*i-3),q(2*i-2),q(2*i-1),q(2*i),q(2*i+1),q(2*i+2),curvature,deltaL,E*I);
            f(2*i-3:2*i+2) = f(2*i-3:2*i+2) + dF;
            J(2*i-3:2*i+2,2*i-3:2*i+2) = J(2*i-3:2*i+2,2*i-3:2*i+2) + dJ;
        end 

        %Update DOF Vector and Error
        q = q - J\f;
        err = sum(abs(f));
    end
    
    %store new positions and velocities
    u = (q-q0)/dt;
    q0 = q;
    u0 = u;

    %Store Mid-node Velocity
    midNodePosition(k+1) = q(2*midNode);
    midNodeVelocity(k+1) = u(2*midNode);
     
end

%---------N update (comment out this section if not looping N) ---------
    %Record Terminal Velocity for N array 
    terminalVelocity(N-2) = u(end); 

    %Update N array
    NArray(N-2) = N;
end 

%---------dt update (comment out this section if not looping dt) ---------
%     %Record Terminal Velocity for dt array
%     terminalVelocity(a) = u(end);
% 
%     %Update dt array 
%     dtArray(a) = 1-a;
% end 
% -----------------------------------------------------------------------
% Plot the final beam position
figure(1)
clf;
plot(q(1:2:end), q(2:2:end), 'ro-');
axis equal;
xlabel('x (meter)');
ylabel('y (meter)');
drawnow

%plot mid-node position
t = (0:steps-1) * dt; %time array
figure(2)
plot(t,midNodePosition, 'k-'); 
xlabel('t [s]'); 
ylabel('y [m]');

%plot mid-node y-velocity
t = (0:steps-1) * dt; %time array
figure(3)
plot(t,midNodeVelocity, 'k-'); 
xlabel('t [s]'); 
ylabel('u_y [m/s]');
ylim([-0.01 0]);

%plot N vs. terminal y-velocity (comment out if not doing N loop)
figure(4)
plot(NArray,terminalVelocity, 'ko-'); 
xlabel('N'); 
ylabel('y_{terminal} [m/s]');

% %plot dt vs. terminal y-velocity (comment out if not doing dt loop)
% figure(5)
% plot(dtArray,terminalVelocity, 'ko-'); 
% xlabel('log(dt)'); 
% ylabel('y_{terminal} [m/s]');

% %%Explicit Method-----------------------------------------------
% 
% dt = 1e-6; %decrease dt for explicit method to converge
% 
% %Initial DOF vector
% q0 = zeros(ndof,1);
% for c=1:N
%     q0(2*c-1) = nodes(c,1); %initial x configuration
%     q0(2*c) = nodes(c,2); %initial y configuration
% end
% 
% %Initial Velocity Vector
% u0 = zeros(ndof,1); %velocity vector
% 
% %Number of Steps
% steps = round(maxTime / dt);
% 
% %Middle-node Velocity and Position Array
% midNodePosition = zeros(steps,1);
% midNodeVelocity = zeros(steps,1); 
% 
% 
% %Time stepping Scheme
% for k = 1:steps-1
%     q = q0; %update old DOFS
%     u = u0; %update old velocities
%     
%     %Initialize Gradient of Elastic Energy Vector
%     dEdq = zeros(ndof,1);
% 
%     %Linear Springs
%     for i = 1:N-1
%         dF = gradEs(q(2*i-1),q(2*i),q(2*i+1),q(2*i+2),deltaL,E*A);
%         dEdq(2*i-1:2*i+2) = dEdq(2*i-1:2*i+2) + dF;
%     end
% 
%     %Bending Springs
%     for i = 2:N-1
%         dB = gradEb(q(2*i-3),q(2*i-2),q(2*i-1),q(2*i),q(2*i+1),q(2*i+2),curvature,deltaL,E*I);
%         dEdq(2*i-3:2*i+2) = dEdq(2*i-3:2*i+2) + dB;
%     end
%     
%     %Update DOF Vector
%     q = q + dt * (u - dt * (M \ (W + C*u + dEdq)));
% 
%     %store new positions and velocities
%     u = (q-q0)/dt;
%     q0 = q;
%     u0 = u;
% 
%     %Store Mid-node Velocity
%     midNodePosition(k+1) = q(2*midNode);
%     midNodeVelocity(k+1) = u(2*midNode);
% end
% 
% %plot beam position (not real time, since explicit is slow)
% figure(6)
% % clf;
% plot(q(1:2:end), q(2:2:end), 'ro-');
% axis equal;
% xlabel('x (meter)');
% ylabel('y (meter)');
% 
% %plot mid-node position
% t = (0:steps-1) * dt; %time array
% figure(7)
% plot(t,midNodePosition, 'k-'); 
% xlabel('t [s]'); 
% ylabel('y [m]');
% 
% %plot mid-node y-velocity
% t = (0:steps-1) * dt; %time array
% figure(8)
% plot(t,midNodeVelocity, 'k-'); 
% xlabel('t [s]'); 
% ylabel('u_y [m/s]');
% ylim([-0.01 0]);
