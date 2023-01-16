clear variables; 
clc;

%Section 1: Setup Variables

%For time,
tFinal = 2; 
t = 0;

%For Boundary Conditions
a = -1; 
b = 1; 
c = -1; 
d = 1;

%Number of points
D = 0.0005;
nX = 100; 
nY = 100;

%For Velocity Equations
Vx = @(t,x,y)  cos((pi*x)/2)*cos((pi*x)/2)*sin(pi*y);
Vy = @(t,x,y) -cos((pi*y)/2)*(-cos((pi*y)/2))*sin(pi*x);

exactU = @(t, x, y)   exp(-t)*(sin(x) + sin(y));
S      = @(t, x, y) 0                       
fc     = @(t, x, y)   0.1%D*exp(-t)*cos(y);
fd     = @(t, x, y)   0.2%-D*exp(-t)*cos(y);

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
dt     = 0.5*dx;

%Section 2: Loops and Plots

%for looping statement involving exact differential and un
for i=1:nX
   for j=1:nY
	un(i,j) = 15*x(i) + 75;
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
    
    % Build the linear system A*unp1 = rhs
    % Boundary conditions:
    % Left wall:
    for j=1:nY
        i = 1;
        p = (j-1)*nX + i;
        A(p,p) = 1; rhs(p) = 60 ;%exactU(t+dt, a, y(j));
    end
    % Right wall:
    for j=1:nY
        i = nX;
        p = (j-1)*nX + i;
        A(p,p) = 1; rhs(p) = 90;%exactU(t+dt, x(i), y(j));
    end
    % Bottom wall:
    for i=2:nX-1
        j = 1;
        p = (j-1)*nX + i;
            A(p,p   ) = C;
            A(p,p-1 ) = L;
            A(p,p+1 ) = R;
            A(p,p+nX) = B+T;
            if Vx(t, x(i), y(j)) > 0
                dutimesdx = (un(i,j)-un(i-1,j))/dx;
            else
                dutimesdx = (un(i+1,j)-un(i,j))/dx;
            end
            rhs(i,j) = un(i,j) + dt*S(t, x(i), y(j))- Vx(t, x(i), y(j))*dutimesdx + 2*B*dy/D*fc(t, x(i), y(j));
            rhs(p) = rhs(i,j);
    end
    % Top wall:
    for i=2:nX-1
        j = nY;
        p = (j-1)*nX + i;
            A(p,p   ) = C;
            A(p,p-1 ) = L;
            A(p,p+1 ) = R;
            A(p,p-nX) = B+T;
            if Vx(t, x(i), y(j)) > 0
                dutimesdx = (un(i,j)-un(i-1,j))/dx;
            else
                dutimesdx = (un(i+1,j)-un(i,j))/dx;
            end
            rhs(i,j) = un(i,j) + dt*S(t, x(i), y(j))- Vx(t, x(i), y(j))*dutimesdx + 2*T*dy/D*fd(t, x(i), y(j));
            rhs(p) = rhs(i,j);
    end
    
    % Interior points:
    for i=2:nX-1
        for j=2:nY-1
            % Index p of the row associated with the grid point
            % (i,j):
            p = (j-1)*nX + i;
            A(p,p   ) = C;
            A(p,p-1 ) = L;
            A(p,p+1 ) = R;
            A(p,p-nX) = B;
            A(p,p+nX) = T;
 
             if Vx(t, x(i), y(j)) > 0 
                 dutimesdx = (un(i,j)-un(i-1,j))/dx;
             else
                 dutimesdx = (un(i+1,j)-un(i,j))/dx;
             end 
             
             if Vy(t, x(i), y(j)) > 0 
                 dutimesdy = (un(i,j)-un(i,j-1))/dy;
             else
                 dutimesdy = (un(i,j+1)-un(i,j))/dy;
             end
             rhs(i,j) = un(i,j)-Vx(t, x(i), y(j))*dt*dutimesdx-Vy(t, x(i), y(j))*dt*dutimesdy + dt*S(t,x(i),y(j));
             rhs(p) = rhs(i,j);
        end
    end
     
    % Solve the linear system:
    u = A \ rhs;
    
    for i=1:nX
        for j=1:nY
            p = (j-1)*nX + i;
            unp1(i,j) = u(p);
        end
    end
    
%Creates updates for time variables with respect to time step
t = t+dt;
un = unp1
    
% Graphs the results
cla;
surf(x, y, un');
s = sprintf('Solution at t=%2.2f', t);
title(s); xlabel('x'); ylabel('y'); zlabel('u');
axis([a b c d 0 90]); colormap('hot'); shading interp;
pause(.1)
    
end

maxErrorrate = max(max(abs(un - uExact)));
fprintf('The maximum error is %2.8f \n', maxErrorrate);
