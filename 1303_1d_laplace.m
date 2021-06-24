clear,figure(1),clf 

% Mesh parameters
Nx = 16;                 % number of cells
dx = 1/Nx;               % mesh step size
xc = dx/2 : dx : 1-dx/2; % cell centers (x_{i+1/2})
xn = 0    : dx : 1;      % nodes (x_i)

% Time parameters
tbeg = 0;
tend = 1000000;
dt   = 0.25 * dx^2;             % time step size for stability
Nt   = (tend-tbeg) / dt;

% Boundary conditions for u
a = 1; % u(0) 
b = 0; % u(1)

% Initial conditions
y = ones(Nx,1) * 0.5;   % zero initial temperature
q = zeros(Nx+1,1);

% Convergence parameters
eps = 1e-3; 
err0 = 0;

for it = 1:Nt
    
    % find fluxes at inner nodes
    q(2:end-1) = -diff(y)/dx;
    
    % find boundary fluxes 
    q(1)       = -(y(1) - a)   / (dx/2);
    q(end)     = -(b - y(end)) / (dx/2);
    
    
    % du/dt + diff(q)/dx = 0
    y     = y - diff(q)/dx * dt;
  
    if mod(it,1)==0
        plot(xc,y,'-d'),title(it*dt),drawnow
        
        % check error
        err = max(abs(diff(q)/dx));
        if it == 1
          err0 = err;
        end
        printf("iter %d, r = %e (%e)\n", it, err, err/err0);
        if err < eps
            printf("Converged in %d steps", it);
            break;
        end 
    end
end