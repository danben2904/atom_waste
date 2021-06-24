clear,figure(1),clf 

% Mesh parameters
Nx = 2^9;                % number of cells
dx = 1/Nx;               % mesh step size
xc = dx/2 : dx : 1-dx/2; % cell centers (x_{i+1/2})
xn = 0    : dx : 1;      % nodes (x_i)

% Transport parameters
%q = 1;
q = rand(Nx-1,1);
%q = (xn(2:end-1))';


% Time parameters
tbeg = 0;
tend = 1000000;
dt   = dx / max(abs(q));             % time step size for stability
Nt   = (tend-tbeg) / dt;

% Initial conditions
C   = zeros(Nx,1);   % zero initial concentration
qC  = zeros(Nx+1,1); % approximation of (qC) on nodes
for i = 1:Nx
  x = (i+0.5)*dx;
  if x > 0.4 && x < 0.6
    C(i) = 1;
  endif
endfor


plot(xc,C,'-d'),title('Initial'),drawnow


% Action
for it = 1:Nt
    Cn          = (q>0) .* C(1:end-1) + (q<0) .* C(2:end); % conc at nodes: upwind
    qC(2:end-1) = q .* Cn;
    qC(1)       = 0; % Zero conc
    qC(end)     = 0; % Zero conc
    adv         = diff(qC)/dx; 
    C           = C - dt * adv;
    if mod(it,10) == 1
        plot(xc,C,'-d'),title(it*dt),drawnow
    endif
endfor
