clear variables; 
clc;

%Section 1: Setup Variables

%For time,
tFinal = 1; 
t = 0;

%For Boundary Conditions
a = -1; 
b = 1; 
c = -1; 
d = 1;

%Number of points
D = 0.0005;
nX = 160; 
nY = 160;

%For Velocity Equations
vx = @(t,x,y)  cos((pi*x)/2)*cos((pi*x)/2)*sin(pi*y);
vy = @(t,x,y) -cos((pi*y)/2)*(-cos((pi*y)/2))*sin(pi*x);

%For exact differential equation, S, fc and fd
exactU = @(t, x, y)   exp(-t)*(sin(x) + sin(y));
S      = @(t, x, y) -exp(-t)*(sin(x)+sin(y)) ... 
    + (cos((pi*x)/2)*cos((pi*x)/2)*sin(pi*y))*exp(-t)*cos(x) + ...
    (-cos((pi*y)/2)*(-cos((pi*y)/2))*sin(pi*x))*...
     exp(-t)*cos(y) + D*exp(-t)*sin(x) + D*exp(-t)*sin(y);                     
fc     = @(t, x, y)   D*exp(-t)*cos(y);
fd     = @(t, x, y)  -D*exp(-t)*cos(y);

%declaring x and y variables
x = linspace(a, b, nX);
dx = x(2)-x(1);

y = linspace(c, d, nX);
 dy = y(2)-y(1);


%Other declarations for velocity, etc.
un     = zeros(nX,     nY);
uExact = zeros(nX,     nY);
unp1   = zeros(nX   ,    nY);
rhs    = zeros(nX*nY,     1);
u      = zeros(nX*nY,     1);
A      = sparse(nX*nY, nX*nY);
dt     = 0.5*min(dx,dy);

%Section 2: Loops and Plots

%for looping statement involving exact differential and un
for i = 1:nX
    for j = 1:nY
        un(i,j)     = exactU(0, x(i), y(j));
        uExact(i,j) = exactU(0, x(i), y(j));
    end
end

%Plot for graphing results
cla;
surf(x, y, un');
s = sprintf('Solution at t=%2.2f', t);
title(s); xlabel('x'); ylabel('y'); zlabel('u');
axis([a b c d 0 90]); colormap('hot'); shading interp;
pause(.1)

%Declarations for derivatives and C,L,B,R,T
dtOverDx = dt/dx/dx;
dtOverDy = dt/dy/dy;

C = (1 + 2*D*dtOverDx + 2*D*dtOverDy); 
L = -D*dtOverDx; 
R = L;
B = -D*dtOverDy; 
T = B;

%main while loop to run the program
while t < tFinal
    if t+dt > tFinal
        dt = tFinal - t;
        dtOverDx = dt/dx/dx;
        dtOverDy = dt/dy/dy;
        C = (1 + 2*D*dtOverDx + 2*D*dtOverDy);
        L = -D*dtOverDx; 
        R = L;
        B = -D*dtOverDy; 
        T = B;
    end

    % Boundary conditions:
    % On the Left wall:
    for j=1:nY
        i = 1;
        p = (j-1)*nX + i;
        A(p,p) = 1; rhs(p) = exactU(t+dt, a, y(j));
    end

    % On the Right wall:
    for j=1:nY
        i = nX;
        p = (j-1)*nX + i;
        A(p,p) = 1; rhs(p) =exactU(t+dt, x(i), y(j));
    end
    % On the Bottom wall:
    for i=2:nX-1
        j = 1;
        p = (j-1)*nX + i;
            A(p,p) = C;
            A(p,p-1 ) = L;
            A(p,p+1 ) = R;
            A(p,p+nX) = B+T;
            if vx(t, x(i), y(j)) > 0
                dutimesdx = (un(i,j)-un(i-1,j))/dx;
            else
                dutimesdx = (un(i+1,j)-un(i,j))/dx;
            end
            rhs(i,j) = un(i,j) + dt*S(t, x(i), y(j))- vx(t, x(i), y(j))*dutimesdx + 2*B*dy/D*fc(t, x(i), y(j));
            rhs(p) = rhs(i,j);
    end
    % On the Top wall:
    for i=2:nX-1
        j = nY;
        p = (j-1)*nX + i;
            A(p,p) = C;
            A(p,p-1 ) = L;
            A(p,p+1 ) = R;
            A(p,p-nX) = B+T;
            if vx(t, x(i), y(j)) > 0
                dutimesdx = (un(i,j)-un(i-1,j))/dx;
            else
                dutimesdx = (un(i+1,j)-un(i,j))/dx;
            end
            rhs(i,j) = un(i,j) + dt*S(t, x(i), y(j))- vx(t, x(i), y(j))*dutimesdx + 2*T*dy/D*fd(t, x(i), y(j));
            rhs(p) = rhs(i,j);
    end
    
    % Interior points:
    for i=2:nX-1
        for j=2:nY-1
            % Index p of the row associated with the grid point
            % (i,j):
            p = (j-1)*nX + i;
            A(p,p) = C;
            A(p,p-1 ) = L;
            A(p,p+1 ) = R;
            A(p,p-nX) = B;
            A(p,p+nX) = T;
 
             if vx(t, x(i), y(j)) > 0 
                 dutimesdx = (un(i,j)-un(i-1,j))/dx;
             else
                 dutimesdx = (un(i+1,j)-un(i,j))/dx;
             end 
             
             if vy(t, x(i), y(j)) > 0 
                 dudy = (un(i,j)-un(i,j-1))/dy;
             else
                 dudy = (un(i,j+1)-un(i,j))/dy;
             end
             rhs(i,j) = un(i,j)-vx(t, x(i), y(j))*dt*dutimesdx-vy(t, x(i), y(j))*dt*dudy + dt*S(t,x(i),y(j));
             rhs(p) = rhs(i,j);
        end
    end
     


    % Solve for the linear system:
    u = A \ rhs;
    
    for i=1:nX
        for j=1:nY
            p = (j-1)*nX + i;
            unp1(i,j) = u(p);
        end
    end
    
    
%Creates updates for time variables with respect to time step
t = t+dt;
un = unp1;
    
% Graphs the results
cla;
surf(x, y, un');
s = sprintf('Solution at t=%2.2f', t);
title(s); xlabel('x'); ylabel('y'); zlabel('u');
axis([a b c d 0 90]); colormap('hot'); shading interp;
pause(.1)
    
    for i = 1:nX
        for j = 1:nY
            uExact(i,j) = exactU(t, x(i), y(j));
        end
    end
end

maxErrorrate = max(max(abs(un - uExact)));
fprintf('The maximum error is %2.8f \n', maxErrorrate);

