clear,figure(1),clf 

% Mesh parameters
Nx = 64;                 % number of cells
dx = 1/Nx;               % mesh step size
xc = dx/2 : dx : 1-dx/2; % cell centers (x_{i+1/2})
xn = 0    : dx : 1;      % nodes (x_i)

% Time parameters
tbeg = 0;
tend = 1000000;
dt   = 0.25 * dx;             % time step size for stability
Nt   = (tend-tbeg) / dt;

% Boundary conditions for u
a = 1; % u(0) 
b = 0; % u(1)

% Initial conditions
%P  = rand(Nx,1);
P  = ones(Nx,1) * 0.5;   % zero initial pressure
V  = zeros(Nx+1,1);      % velocity
%dv = y*0;

% Convergence parameters
eps = 1e-3; 

for it = 1:Nt
    
    % find velocity at inner nodes
    V(2:end-1) = (1-5/Nx)*V(2:end-1) - dt*diff(P)/dx;
    % find boundary velocity
    % dV/dt = -grad(P) - c*V
    % V^{n+1} - V^n = 
    % V(1)
    V(1)       =  - dt*(P(1) - a)   / (dx/2);
    V(end)     =  - dt*(b - P(end)) / (dx/2);
    
    %dv         = dv*(1-5/Nx) - diff(q)/dx * dt;
    
    % dP/dt + diff(V)/dx = 0
    P     = P - diff(V)/dx*dt;
  
    if mod(it,10)==0
       plot(xc,P,'-d'),title(it*dt),drawnow
        
       % check error
       err = max(abs(diff(V)/dx))
       if err < eps
         printf("Converged in %d steps", it);
         break;
       end 
    end
end

%plot(x,C,'-g',x,C_it,'-d'),title(time),drawnow
%plot(xc,y,'-d'),title(time),drawnow