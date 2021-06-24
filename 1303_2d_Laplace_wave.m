clear,figure(1),clf,colormap(jet)
% Physics
Lx  = 1;
Ly  = 1;
Dx  = 1;%1e-1; % diffusion coefficient
Dy  = 1e-1*Dx;

% Numerics
nx  = 2^8; % number of cells
ny  = nx;

dx  = Lx/nx
dy  = Ly/ny

%preprocessing
[x y]     = ndgrid(dx/2:dx:Lx-dx/2, dy/2:dy:Ly-dy/2);
[x_v y_v] = ndgrid(0   :dx:Lx,      dy/2:dy:Ly-dy/2);
[x_h y_h] = ndgrid(dx/2:dx:Lx-dx/2, 0   :dy:Ly     );
niter     = 3000000;

% Initial conditions
P      = zeros(nx, ny); % Pressure
Vx     = zeros(nx+1,ny);  % x-flux
Vy     = zeros(nx,ny+1);  % y-flux


% Boundary conditions
Pya = rand(ny, 1); Pyb = rand(ny, 1);
Pxa = 0; Pxb = 1;

err0  = 1e20;
 
% Action
tic
tau = min(dx,dy)/4.1;

for nit = 1:niter
  Vx(2:end-1, :) = Vx(2:end-1,:) - tau*diff(P, 1, 1)/dx;
  Vy(:, 2:end-1) = Vy(:,2:end-1) - tau*diff(P, 1, 2)/dy;

     
  % BCs for flux
  Vy(:,1)     = Vy(:,1) - tau*(P(:, 1) - Pya)   / (dy/2);
  Vy(:,end)     = Vy(:,end) - tau*(Pyb - P(:,end)) / (dy/2);

  %Vx(1, :)     = Vx(1, :) - tau*(P(1, :) - Pxa)   / (dx/2);
  %Vx(end, :)     = Vx(end, :) - tau*(Pxb - P(end, :)) / (dx/2);
  
  %Vx(1,:)    = -Kx(1,:).*((P(1,:)-Pa)/(dx/2));
     
  div_q       = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
  r           = div_q;
  err         = max(max(abs(r))); 
  if nit == 1
    err0 = err;
  endif
  if(mod(nit,10)==1  )
    printf("iter %d, r = %e (%e)\n", nit, err, err/err0);
    subplot(221)
    pcolor(x,y,P);shading flat;colorbar;title("Pressure");drawnow;
    subplot(222)
    quiver(Vx(1:end-1, :),Vy(:, 1:end-1));title("Flux");drawnow
  endif
  if(err < 1e-8 || err < 1e-6*err0)
    nit
    err
    break;
  end
  
  %dV  = dV*(1-2/nx) - div_q*tau;
  P   = P - r*tau;
end



toc