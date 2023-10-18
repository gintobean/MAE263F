%Filename: Problem3.m
%Simulation of Elastic Beam Bending

N = 50; %number of nodes
ndof = N * 2; %2-D problem

%Physical Values
rhoMetal = 2700; %aluminum density kg/m3
g = 9.8; %gravity acceleration m/s2
l = 1; %length of beam m
c = 0.75; %location of force from origin mclc
deltaL = l / (N-1);
R = 0.013; %outer radius m
r = 0.011; %inner radius m
Py = 0; %external vertical force in N
E = 70e9; %elastic modulus of beam Pa
A = pi() * (R^2 - r^2); %cross sectional area of beam m2
I = pi() * (R^4 - r^4) / 4; %moment of inertia beam m4
dt = 1e-2; %time step
eps = E*I/l^2 * 1e-3; %tolerance based on small force
maxTime = 1;

%Geometry and initial configuration
nodes = zeros(N,2);
for i=1:N % %loop over nodes
    nodes(i,1) = (i-1) * deltaL; %initial x config
    nodes(i,2) = 0; %initial y config
end

%Mass Matrix
M = zeros(ndof, ndof);
for i = 1:N
    M(i,i) = pi*(R^2-r^2)*l*rhoMetal/(N-1); %mass of spheres kg
    M(2*i,2*i) = M(i,i);
end

%Weight matrix
W = zeros(ndof,1);
for i = 1:N
    W(2*i) = M(2*i,2*i)*g; %weight of spheres
end

%Boundary Conditions, Fixed Indices
fixed = [1;2;2*N];
free = 3:2*N-1;

%Initial DOF vector

q0 = zeros(ndof,1);
for i=1:N
    q0(2*i-1) = nodes(i,1); %make the q vector
    q0(2*i) = nodes(i,2);
end

u0 = zeros(ndof,1); %velocity vector
steps = round(maxTime / dt); %number of steps
ymax = zeros(steps,1);%ymax array
t = zeros(steps,1); %time array

%Containers used for varying P with ymax
maxP = 10000;
increment = 100;
Psteps = maxP / increment;
yTrue = zeros(Psteps+1,1); %true ymax array
yEuler = zeros(Psteps+1,1); %theoretical (Euler-Bernoulli) ymax array
PyArray = zeros(Psteps+1,1); %external force array

%Iterating over different P values
for j = 1:Psteps+1

    %External Force Matrix
    P = zeros(ndof,1);
    nodeOfAction = round(c / l * (N-1)); %Location that Py acts: location/length * number of edges
    P(2*nodeOfAction) = Py;
    
    %Time Stepping
    for k = 1:steps
        q = q0; %update old DOFS
        qFree = q(free);
        u = u0; %update old velocities
        err = 10 * eps; %starting error
        while err > eps
            f = M / dt * ((q-q0)/dt - u) + W + P; 
            J = M / dt^2;
        
            for i = 1:N-1
                %Linear spring
                dF = gradEs(q(2*i-1),q(2*i),q(2*i+1),q(2*i+2),deltaL,E*A);
                dJ = hessEs(q(2*i-1),q(2*i),q(2*i+1),q(2*i+2),deltaL,E*A);
                f(2*i-1:2*i+2) = f(2*i-1:2*i+2) + dF;
                J(2*i-1:2*i+2,2*i-1:2*i+2) = J(2*i-1:2*i+2,2*i-1:2*i+2) + dJ;
            end
            for i = 2:N-1
                %Bending spring
                curvature = 0;
                dF = gradEb(q(2*i-3),q(2*i-2),q(2*i-1),q(2*i),q(2*i+1),q(2*i+2),curvature,deltaL,E*I);
                dJ = hessEb(q(2*i-3),q(2*i-2),q(2*i-1),q(2*i),q(2*i+1),q(2*i+2),curvature,deltaL,E*I);
                f(2*i-3:2*i+2) = f(2*i-3:2*i+2) + dF;
                J(2*i-3:2*i+2,2*i-3:2*i+2) = J(2*i-3:2*i+2,2*i-3:2*i+2) + dJ;
            end 

            %Update
            fFree = f(free);
            JFree = J(free,free);
            fFixed = f(fixed);
            qFree = qFree - JFree\fFree;
            q(free) = qFree;
            err = sum(abs(fFree));
        end
    
        %store new positions and velocities
        u = (q-q0)/dt;
        q0 = q;
        u0 = u;

        %At Py=2000, find relationship between ymax and time
        if Py == 2000
            disp(abs(min(q)));
            ymax(k) = abs(min(q)); 
            t(k) = (k-1) * dt; %t(k=1) = 0
        end
    end

    %Store values of Py and corresponding ymax
    PyArray(j) = Py;
    yTrue(j) = abs(min(q));
    yEuler(j) = Py*(l-c)*(l^2-(l-c)^2)^1.5 / (9*sqrt(3)*E*I*l);
    Py = Py + increment;
end

% Plot ymax as a function of time.
figure(1)
clf;
plot(t, ymax, 'r*-');
axis equal;
xlabel('t (s)');
ylabel('y (meter)');
drawnow

%Plot theoretical vs. true ymax
figure(2)
plot(PyArray,yEuler, 'ko-');
hold on
plot(PyArray,yTrue, 'bo-');
hold off
xlabel('P_y (N)');
ylabel('y_(max) (m)');
