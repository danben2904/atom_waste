clear,figure(1),clf,colormap(jet)
% Physics
Lx  = 1;
Ly  = 1;
Dx  = 1;          %1e-1; % diffusion coefficient
Dy  = 1e-1*Dx;

% Numerics
nx  = 2^4; % number of cells
ny  = nx;

dx  = Lx/nx
dy  = Ly/ny

K  = ones(nx,ny);

Kx = Dx*ones(nx+1,ny);
Ky = Dy*ones(nx,ny+1);





%preprocessing
[x y]     = ndgrid(dx/2: dx : Lx-dx/2, dy/2 : dy : Ly-dy/2);
[x_v y_v] = ndgrid(0   :dx:Lx,      dy/2:dy:Ly-dy/2);
[x_h y_h] = ndgrid(dx/2:dx:Lx-dx/2, 0   :dy:Ly     );
niter     = 3000000;

% Initial conditions
P      = ones(nx,ny) * 0; % Pressure
qx     = zeros(nx+1,ny);  % x-flux
qy     = zeros(nx,ny+1);  % y-flux


% Boundary conditions
Pa = 0; Pb = 1;

err0  = 1e20;
 
% Action
tic
tau = min(dx^2,dy^2)/4.1;
% tau = tau / sqrt(K)
tau = tau / sqrt(max(max(Kx(:)),max(Ky(:))));

for nit = 1:niter
  qx(2:end-1,:)      = -Kx(2:end-1,:) .* diff(P,1,1)/dx; % -dP/dx
  qy(:,2:end-1)      = -Ky(:,2:end-1) .* diff(P,1,2)/dy; % -dP/dy
     
  % BCs for flux
  qy(:,1)     = -Ky(:,1)   .* (P(:,1) -   Pa) / (dy/2);... % P = Pa
  qy(:,end)   = -Ky(:,end) .* (Pb - P(:,end)) / (dy/2); % P = Pb
  
  %qx(1,:)     = -Kx(1,:).*((P(1,:)-Pa)/(dx/2));
     
  div_q       = diff(qx,1,1)/dx + diff(qy,1,2)/dy;
  r           = div_q;
  err         = max(max(abs(r))); 
  if nit == 1
    err0 = err;
  endif
  if(mod(nit,100)==1)
    printf("iter %d, r = %e (%e)\n", nit, err, err/err0);
    subplot(221)
    pcolor(x,y,P);shading flat;colorbar;title("Pressure");drawnow;
    subplot(222)
    %quiver(qx,qy);title("Flux");drawnow
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