clear,figure(1),clf,colormap(jet)
% Physics
Lx  = 1;
Ly  = 1;

% Numerics
nx  = 2^7; % number of cells
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
qx     = zeros(nx+1,ny);  % x-flux
qy     = zeros(nx,ny+1);  % y-flux


% Boundary conditions
Pya = rand(ny, 1); Pyb = rand(ny, 1);
Pxa = rand(1, nx); Pxb = rand(1, nx);

err0  = 1e20;
 
% Action
tic
dt = min(dx,dy)/4.1;

Kf = (1-5/nx);

for nit = 1:niter
  qx(2:end-1, :) = Kf*qx(2:end-1,:) - dt*diff(P, 1, 1)/dx;
  qy(:, 2:end-1) = Kf*qy(:,2:end-1) - dt*diff(P, 1, 2)/dy;

     
  % BCs for flux
  qy(:,1)     = Kf*qy(:,1)   - dt*(P(:, 1) - Pya)  / (dy/2);
  qy(:,end)   = Kf*qy(:,end) - dt*(Pyb - P(:,end)) / (dy/2);

  qx(1, :)     = qx(1, :) - dt*(P(1, :) - Pxa)   / (dx/2);
  qx(end, :)     = qx(end, :) - dt*(Pxb - P(end, :)) / (dx/2);
  
  %qx(1,:)    = -Kx(1,:).*((P(1,:)-Pa)/(dx/2));
     
  div_q       = diff(qx,1,1)/dx + diff(qy,1,2)/dy;
  r           = div_q;
  err         = max(max(abs(r))); 
  if nit == 1
    err0 = err;
  endif
  %if(mod(nit,10)==1  )
  %  printf("iter %d, r = %e (%e)\n", nit, err, err/err0);
  %  subplot(221)
  %  pcolor(x,y,P);shading flat;colorbar;title("Pressure");drawnow;
  %  subplot(222)
  %  quiver(qx(1:end-1, :),qy(:, 1:end-1));title("Flux");drawnow
  %endif
  if(err < 1e-8 || err < 1e-4*err0)
    nit
    err
    break;
  endif
  
  %dV  = dV*(1-2/nx) - div_q*dt;
  P   = P - r*dt;
endfor

subplot(221)
pcolor(x,y,P);shading flat;colorbar;title("Pressure");drawnow;


dt   = dx / max(max(abs(qy(:))), max(abs(qx(:))));
Nt   = niter / dt;

C    = zeros(nx,ny);
qCx  = zeros(nx+1,ny);
qCy  = zeros(nx,ny+1);

for i = 1:nx
  xx = (i+0.5)*dy;
  if xx > 0.7 && xx < 0.9
    for j = 1:ny
      yy = (j+0.5)*dy;
      if yy > 0.6 && yy < 0.8
        C(i, j) = 1;
      endif
    endfor
  endif
endfor

for it = 1:Nt
  Cx              = (qx(2:end-1, :)>0) .* C(1:end-1, :) + (qx(2:end-1, :)<0) .* C(2:end, :);
  qCx(2:end-1, :) = qx(2:end-1, :) .* Cx;
  qCx(1,   :)     = 0;
  qCx(end, :)     = 0;
  Cy              = (qy(:, 2:end-1)>0) .* C(:, 1:end-1) + (qy(:, 2:end-1)<0) .* C(:, 2:end);
  qCy(:, 2:end-1) = qy(:, 2:end-1) .* Cy;
  qCy(:,   1)     = 0;
  qCy(:, end)     = 0;
  adv             = diff(qCx, 1, 1)/dx + diff(qCy, 1, 2)/dy; 
  C               = C - dt * adv;
  if mod(it,200) == 0
    subplot(223)
    pcolor(x,y,C);shading flat;colorbar;title("Water");drawnow;
  endif
endfor


toc