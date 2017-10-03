clear global; clear functions;
clear all;
global a a_x a_x_half a_y a_y_half abscissa_in_pml abscissa_normalized alpha_max_pml alpha_prime_x alpha_prime_x_half alpha_prime_y alpha_prime_y_half angle_force b_x b_x_half b_y b_y_half courant_number cp cs d0_x d0_y d_x d_x_half d_y d_y_half degrees_to_radians deltat deltax deltay density dist distval epsilon_xx epsilon_xy epsilon_yy f0 factor force_x force_y hugeval i irec isource it it_display ix_rec iy_rec j jsource k_max_pml k_x k_x_half k_y k_y_half lambda lambda_half_x lambda_plus_two_mu_half_x memory_dsigmaxx_dx memory_dsigmaxy_dx memory_dsigmaxy_dy memory_dsigmayy_dy memory_dvx_dx memory_dvx_dy memory_dvy_dx memory_dvy_dy mu mu_half_x mu_half_y npoints_pml npower nrec nstep nx ny pi rcoef rho rho_half_x_half_y sigmaxx sigmaxy sigmayy sisvx sisvy source_term stability_threshold t t0 thickness_pml_x thickness_pml_y total_energy_kinetic total_energy_potential use_pml_xmax use_pml_xmin use_pml_ymax use_pml_ymin value_dsigmaxx_dx value_dsigmaxy_dx value_dsigmaxy_dy value_dsigmayy_dy value_dvx_dx value_dvx_dy value_dvy_dx value_dvy_dy velocnorm vx vy xdeb xfin xoriginleft xoriginright xrec xsource xspacerec xval ydeb yfin yoriginbottom yorigintop yrec ysource yspacerec yval zero ; 


sourcetype=input(' [1]  Gaussian [2] First Derivative of Gaussian [3] Second Derivative of Gaussian ' ,'s')

%%
if isempty(nx), nx = 100; end;
if isempty(ny), ny = nx; end;
% size of a grid cell
if isempty(deltax), deltax = 1.0d0; end;
if isempty(deltay), deltay = deltax; end;
% flags to add PML layers to the edges of the grid
if isempty(use_pml_xmin), use_pml_xmin = true; end;
if isempty(use_pml_xmax), use_pml_xmax = true; end;
if isempty(use_pml_ymin), use_pml_ymin = true; end;
if isempty(use_pml_ymax), use_pml_ymax = true; end;
% thickness of the PML layer in grid points
if isempty(npoints_pml), npoints_pml = nx/10; end;
% P-velocity, S-velocity and density
if isempty(cp), cp = 5956000d0; end;
if isempty(cs), cs = 3227070d0; end;
if isempty(density), density = 0.007874d0; end;

cp2= 500000
cs2= 1000000
density2= 0.002000

% display information on the screen from time to time
if isempty(it_display), it_display = 100; end;
temp_step=1

% total number of time steps
if isempty(nstep), nstep = 20001; end;
% time step in seconds
if isempty(deltat), deltat = 2.0d-8; end;
% parameters for the source
if isempty(f0), f0 = 80000.0d0; end;
if isempty(t0), t0 = 1.20d0 ./ f0; end;
if isempty(factor), factor = 1.0d2; end;
% source
if isempty(isource), isource = fix(nx./2)+1; end;
if isempty(jsource), jsource = fix(ny./2)+1; end;
if isempty(xsource), xsource =(isource - 1) .* deltax; end;
if isempty(ysource), ysource =(jsource - 1) .* deltay; end;
% angle of source force clockwise with respect to vertical (Y) axis
if isempty(angle_force), angle_force = 90.0d0; end;
% receivers
if isempty(nrec), nrec = 2; end;


%radius= fix(nx/3)+1
radius=30

% first receiver x in meters
if isempty(xdeb), xdeb = xsource; end;
fprintf(['%s %0.15g \n'], 'x first receiver = ',xdeb);

% first receiver y in meters
if isempty(ydeb), ydeb = ysource; end;
fprintf(['%s %0.15g \n'], 'y first receiver = ',ydeb);

% last receiver x in meters
if isempty(xfin), xfin = xsource; end;
fprintf(['%s %0.15g \n'], 'x second receiver  = ',xfin);

% last receiver y in meters
if isempty(yfin), yfin =  (jsource - radius) .* deltay; end;
fprintf(['%s %0.15g \n'], 'y second receiver = ',yfin);



%%%%%%%%%%%%%%%%%%%%5
% value of PI
if isempty(pi), pi = 3.141592653589793238462643d0; end;
% conversion from degrees to radians
if isempty(degrees_to_radians), degrees_to_radians = pi ./ 180.0d0; end;
% zero
if isempty(zero), zero = 0.0d0; end;
% large value for maximum
if isempty(hugeval), hugeval = 1.0d+30; end;
% velocity threshold above which we consider that the code became unstable
if isempty(stability_threshold), stability_threshold = 1.0d+25; end;
% main arrays
if isempty(vx), vx=zeros(nx,ny); end;
if isempty(vy), vy=zeros(nx,ny); end;
if isempty(sigmaxx), sigmaxx=zeros(nx,ny); end;
if isempty(sigmayy), sigmayy=zeros(nx,ny); end;
if isempty(sigmaxy), sigmaxy=zeros(nx,ny); end;
if isempty(lambda), lambda=zeros(nx,ny); end;
if isempty(mu), mu=zeros(nx,ny); end;
if isempty(rho), rho=zeros(nx,ny); end;
% to interpolate material parameters at the right location in the staggered grid cell
if isempty(lambda_half_x), lambda_half_x=0; end;
if isempty(mu_half_x), mu_half_x=0; end;
if isempty(lambda_plus_two_mu_half_x), lambda_plus_two_mu_half_x=0; end;
if isempty(mu_half_y), mu_half_y=0; end;
if isempty(rho_half_x_half_y), rho_half_x_half_y=0; end;
% for evolution of total energy in the medium
if isempty(epsilon_xx), epsilon_xx=0; end;
if isempty(epsilon_yy), epsilon_yy=0; end;
if isempty(epsilon_xy), epsilon_xy=0; end;
if isempty(total_energy_kinetic), total_energy_kinetic=zeros(1,nstep/1000); end;
if isempty(total_energy_potential), total_energy_potential=zeros(1,nstep/1000); end;
% power to compute d0 profile
if isempty(npower), npower = 2.0d0; end;
% from Gedney page 8.11
if isempty(k_max_pml), k_max_pml = 1.0d0; end;
% from festa and Vilotte
if isempty(alpha_max_pml), alpha_max_pml = 2.0d0.*pi.*(f0./2.0d0); end;
% arrays for the memory variables
% could declare these arrays in PML only to savemlv a lot of memory, but proof of concept only here
if isempty(memory_dvx_dx), memory_dvx_dx=zeros(nx,ny); end;
if isempty(memory_dvx_dy), memory_dvx_dy=zeros(nx,ny); end;
if isempty(memory_dvy_dx), memory_dvy_dx=zeros(nx,ny); end;
if isempty(memory_dvy_dy), memory_dvy_dy=zeros(nx,ny); end;
if isempty(memory_dsigmaxx_dx), memory_dsigmaxx_dx=zeros(nx,ny); end;
if isempty(memory_dsigmayy_dy), memory_dsigmayy_dy=zeros(nx,ny); end;
if isempty(memory_dsigmaxy_dx), memory_dsigmaxy_dx=zeros(nx,ny); end;
if isempty(memory_dsigmaxy_dy), memory_dsigmaxy_dy=zeros(nx,ny); end;
if isempty(value_dvx_dx), value_dvx_dx=0; end;
if isempty(value_dvx_dy), value_dvx_dy=0; end;
if isempty(value_dvy_dx), value_dvy_dx=0; end;
if isempty(value_dvy_dy), value_dvy_dy=0; end;
if isempty(value_dsigmaxx_dx), value_dsigmaxx_dx=0; end;
if isempty(value_dsigmayy_dy), value_dsigmayy_dy=0; end;
if isempty(value_dsigmaxy_dx), value_dsigmaxy_dx=0; end;
if isempty(value_dsigmaxy_dy), value_dsigmaxy_dy=0; end;
% 1D arrays for the damping profiles
if isempty(d_x), d_x=zeros(1,nx); end;
if isempty(k_x), k_x=zeros(1,nx); end;
if isempty(alpha_prime_x), alpha_prime_x=zeros(1,nx); end;
if isempty(a_x), a_x=zeros(1,nx); end;
if isempty(b_x), b_x=zeros(1,nx); end;
if isempty(d_x_half), d_x_half=zeros(1,nx); end;
if isempty(k_x_half), k_x_half=zeros(1,nx); end;
if isempty(alpha_prime_x_half), alpha_prime_x_half=zeros(1,nx); end;
if isempty(a_x_half), a_x_half=zeros(1,nx); end;
if isempty(b_x_half), b_x_half=zeros(1,nx); end;
if isempty(d_y), d_y=zeros(1,ny); end;
if isempty(k_y), k_y=zeros(1,ny); end;
if isempty(alpha_prime_y), alpha_prime_y=zeros(1,ny); end;
if isempty(a_y), a_y=zeros(1,ny); end;
if isempty(b_y), b_y=zeros(1,ny); end;
if isempty(d_y_half), d_y_half=zeros(1,ny); end;
if isempty(k_y_half), k_y_half=zeros(1,ny); end;
if isempty(alpha_prime_y_half), alpha_prime_y_half=zeros(1,ny); end;
if isempty(a_y_half), a_y_half=zeros(1,ny); end;
if isempty(b_y_half), b_y_half=zeros(1,ny); end;
if isempty(thickness_pml_x), thickness_pml_x=0; end;
if isempty(thickness_pml_y), thickness_pml_y=0; end;
if isempty(xoriginleft), xoriginleft=0; end;
if isempty(xoriginright), xoriginright=0; end;
if isempty(yoriginbottom), yoriginbottom=0; end;
if isempty(yorigintop), yorigintop=0; end;
if isempty(rcoef), rcoef=0; end;
if isempty(d0_x), d0_x=0; end;
if isempty(d0_y), d0_y=0; end;
if isempty(xval), xval=0; end;
if isempty(yval), yval=0; end;
if isempty(abscissa_in_pml), abscissa_in_pml=0; end;
if isempty(abscissa_normalized), abscissa_normalized=0; end;
% for the source
if isempty(a), a=0; end;
if isempty(t), t=0; end;
if isempty(force_x), force_x=0; end;
if isempty(force_y), force_y=0; end;
if isempty(source_term), source_term=0; end;
% for receivers
if isempty(xspacerec), xspacerec=0; end;
if isempty(yspacerec), yspacerec=0; end;
if isempty(distval), distval=0; end;
if isempty(dist), dist=0; end;
if isempty(ix_rec), ix_rec=zeros(1,nrec); end;
if isempty(iy_rec), iy_rec=zeros(1,nrec); end;
if isempty(xrec), xrec=zeros(1,nrec); end;
if isempty(yrec), yrec=zeros(1,nrec); end;
% for seismograms
if isempty(sisvx), sisvx=zeros(nstep/1000,nrec); end;
if isempty(sisvy), sisvy=zeros(nstep/1000,nrec); end;
if isempty(i), i=0; end;
if isempty(j), j=0; end;
if isempty(it), it=0; end;
if isempty(irec), irec=0; end;
if isempty(courant_number), courant_number=0; end;
if isempty(velocnorm), velocnorm=0; end;

%%
%---
%--- program starts here
%---
fprintf(['%0.15g \n']);
fprintf(['%s \n'], '2D elastic finite-difference code in velocity and stress formulation with C-PML');
fprintf(['%0.15g \n']);
% display size of the model
fprintf(['%0.15g \n']);
fprintf(['%s %0.15g \n'], 'NX = ',nx);
fprintf(['%s %0.15g \n'], 'NY = ',ny);
fprintf(['%0.15g \n']);
fprintf(['%s %0.15g \n'], 'size of the model along X = ',(nx - 1) .* deltax);
fprintf(['%s %0.15g \n'], 'size of the model along Y = ',(ny - 1) .* deltay);
fprintf(['%0.15g \n']);
fprintf(['%s %0.15g \n'], 'Total number of grid points = ',nx .* ny);
fprintf(['%0.15g \n']);
%--- define profile of absorption in PML region
% thickness of the PML layer in meters
thickness_pml_x = npoints_pml .* deltax;
thickness_pml_y = npoints_pml .* deltay;
% reflection coefficient (INRIA report section 6.1)
rcoef = 0.001d0;
% check that NPOWER is okay
if(npower < 1)
error(['stop encountered in original matlab code  ',char(10),' NPOWER must be greater than 1;']);
end;
% compute d0 from INRIA report section 6.1
d0_x = -(npower + 1) .* cp .* log(rcoef) ./(2.0d0 .* thickness_pml_x);
d0_y = -(npower + 1) .* cp .* log(rcoef) ./(2.0d0 .* thickness_pml_y);
fprintf(['%s %0.15g \n'], 'd0_x = ',d0_x);
fprintf(['%s %0.15g \n'], 'd0_y = ',d0_y);
fprintf(['%0.15g \n']);
d_x(:) = zero;
d_x_half(:) = zero;
k_x(:) = 1.0d0;
k_x_half(:) = 1.0d0;
alpha_prime_x(:) = zero;
alpha_prime_x_half(:) = zero;
a_x(:) = zero;
a_x_half(:) = zero;
d_y(:) = zero;
d_y_half(:) = zero;
k_y(:) = 1.0d0;
k_y_half(:) = 1.0d0;
alpha_prime_y(:) = zero;
alpha_prime_y_half(:) = zero;
a_y(:) = zero;
a_y_half(:) = zero;
% damping in the X direction
% origin of the PML layer (position of right edge minus thickness, in meters)
xoriginleft = thickness_pml_x;
xoriginright =(nx-1).*deltax - thickness_pml_x;
for i = 1:nx;
% abscissa of current grid point along the damping profile
xval = deltax .* (i-1);
%---------- left edge
if(use_pml_xmin)
% define damping profile at the grid points
abscissa_in_pml = xoriginleft - xval;
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_x;
d_x(i) = d0_x .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_x(i) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_x(i) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
% define damping profile at half the grid points
abscissa_in_pml = xoriginleft -(xval + deltax./2.0d0);
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_x;
d_x_half(i) = d0_x .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_x_half(i) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_x_half(i) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
end;
%---------- right edge
if(use_pml_xmax)
% define damping profile at the grid points
abscissa_in_pml = xval - xoriginright;
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_x;
d_x(i) = d0_x .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_x(i) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_x(i) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
% define damping profile at half the grid points
abscissa_in_pml = xval + deltax./2.0d0 - xoriginright;
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_x;
d_x_half(i) = d0_x .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_x_half(i) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_x_half(i) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
end;
% just in case, for -5 at the end
if(alpha_prime_x(i) < zero)
alpha_prime_x(i) = zero;
end;
if(alpha_prime_x_half(i) < zero)
alpha_prime_x_half(i) = zero;
end;
b_x(i) = exp(-(d_x(i) ./ k_x(i) + alpha_prime_x(i)) .* deltat);
b_x_half(i) = exp(-(d_x_half(i) ./ k_x_half(i) + alpha_prime_x_half(i)) .* deltat);
% this to avoid division by zero outside the PML
if(abs(d_x(i)) > 1.0d-6)
a_x(i) = d_x(i) .*(b_x(i) - 1.0d0) ./(k_x(i) .*(d_x(i) + k_x(i) .* alpha_prime_x(i)));
end;
if(abs(d_x_half(i)) > 1.0d-6)
a_x_half(i) = d_x_half(i) .*(b_x_half(i) - 1.0d0) ./(k_x_half(i) .*(d_x_half(i) + k_x_half(i) .* alpha_prime_x_half(i)));
end;
end; i =fix(nx+1);
% damping in the Y direction
% origin of the PML layer (position of right edge minus thickness, in meters)
yoriginbottom = thickness_pml_y;
yorigintop =(ny-1).*deltay - thickness_pml_y;
for j = 1:ny;
% abscissa of current grid point along the damping profile
yval = deltay .* (j-1);
%---------- bottom edge
if(use_pml_ymin)
% define damping profile at the grid points
abscissa_in_pml = yoriginbottom - yval;
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_y;
d_y(j) = d0_y .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_y(j) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_y(j) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
% define damping profile at half the grid points
abscissa_in_pml = yoriginbottom -(yval + deltay./2.0d0);
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_y;
d_y_half(j) = d0_y .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_y_half(j) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_y_half(j) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
end;
%---------- top edge
if(use_pml_ymax)
% define damping profile at the grid points
abscissa_in_pml = yval - yorigintop;
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_y;
d_y(j) = d0_y .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_y(j) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_y(j) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
% define damping profile at half the grid points
abscissa_in_pml = yval + deltay./2.0d0 - yorigintop;
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_y;
d_y_half(j) = d0_y .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_y_half(j) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_y_half(j) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
end;
b_y(j) = exp(-(d_y(j) ./ k_y(j) + alpha_prime_y(j)) .* deltat);
b_y_half(j) = exp(-(d_y_half(j) ./ k_y_half(j) + alpha_prime_y_half(j)) .* deltat);
% this to avoid division by zero outside the PML
if(abs(d_y(j)) > 1.0d-6)
a_y(j) = d_y(j) .*(b_y(j) - 1.0d0) ./(k_y(j) .*(d_y(j) + k_y(j) .* alpha_prime_y(j)));
end;
if(abs(d_y_half(j)) > 1.0d-6)
a_y_half(j) = d_y_half(j) .*(b_y_half(j) - 1.0d0) ./(k_y_half(j) .*(d_y_half(j) + k_y_half(j) .* alpha_prime_y_half(j)));
end;
end; j =fix(ny+1);


% compute the Lame parameters and density
for j = 1:ny;
for i = 1:nx;

rho(i,j) = density;
mu(i,j) = density.*cs.*cs;
lambda(i,j) = density.*(cp.*cp - 2.0d0.*cs.*cs);


if(((i-isource).^2+(j-jsource).^2)>=(radius.^2))

rho(i,j) = density2;
mu(i,j) = density2.*cs2.*cs2;
lambda(i,j) = density2.*(cp2.*cp2 - 2.0d0.*cs2.*cs2);
       
end;



%if(((i-isource).^2+(j-jsource).^2)<=(20.^2))

%rho(i,j) = 1000;
%mu(i,j) = 1000.*1481.*1481;
%lambda(i,j) = 2.*(1400.*1400 - 2.0d0.*1481.*1481);
%end;




end; i =fix(nx+1);
end; j =fix(ny+1);




% print position of the source
fprintf(['%s \n'], 'Position of the source:');
fprintf(['%0.15g \n']);
fprintf(['%s %0.15g \n'], 'x = ',xsource);
fprintf(['%s %0.15g \n'], 'y = ',ysource);
fprintf(['%0.15g \n']);
% define location of receivers
fprintf(['%s %0.15g %s \n'], 'There are ',nrec,' receivers');
fprintf(['%0.15g \n']);
xspacerec =(xfin-xdeb) ./ (nrec-1);
yspacerec =(yfin-ydeb) ./ (nrec-1);
for irec=1:nrec;
xrec(irec) = xdeb + (irec-1).*xspacerec;
yrec(irec) = ydeb + (irec-1).*yspacerec;
end; irec=fix(nrec+1);
% find closest grid point for each receiver
for irec=1:nrec;
dist = hugeval;
for j = 1:ny;
for i = 1:nx;
distval = sqrt((deltax.*(i-1) - xrec(irec)).^2 +(deltay.*(j-1) - yrec(irec)).^2);
if(distval < dist)
dist = distval;
ix_rec(irec) = fix(i);
iy_rec(irec) = fix(j);
end;
end; i =fix(nx+1);
end; j =fix(ny+1);
fprintf(['%s %0.15g %s %0.15g %0.15g \n'], 'receiver ',irec,' x_target,y_target = ',xrec(irec),yrec(irec));
fprintf(['%s %0.15g %s %0.15g %0.15g \n'], 'closest grid point found at distance ',dist,' in i,j = ',ix_rec(irec),iy_rec(irec));
fprintf(['%0.15g \n']);
end; irec=fix(nrec+1);
% check the Courant stability condition for the explicit time scheme
% R. Courant et K. O. Friedrichs et H. Lewy (1928)
courant_number = cp .* deltat .* sqrt(1.0d0./deltax.^2 + 1.0d0./deltay.^2);
fprintf(['%s %0.15g \n'], 'Courant number is ',courant_number);
fprintf(['%0.15g \n']);
if(courant_number > 1.0d0)
error(['stop encountered in original matlab code  ',char(10),' time step is too large, simulation will be unstable;']);
end;
% suppress old files (can be commented out if 'call system' is missing in your compiler)
% call system('rm -f Vx_*.dat Vy_*.dat image*.pnm image*.gif')
% initialize arrays
vx(:,:) = zero;
vy(:,:) = zero;
sigmaxx(:,:) = zero;
sigmayy(:,:) = zero;
sigmaxy(:,:) = zero;
% PML
memory_dvx_dx(:,:) = zero;
memory_dvx_dy(:,:) = zero;
memory_dvy_dx(:,:) = zero;
memory_dvy_dy(:,:) = zero;
memory_dsigmaxx_dx(:,:) = zero;
memory_dsigmayy_dy(:,:) = zero;
memory_dsigmaxy_dx(:,:) = zero;
memory_dsigmaxy_dy(:,:) = zero;
% initialize seismograms
sisvx(:,:) = zero;
sisvy(:,:) = zero;
% initialize total energy
total_energy_kinetic(:) = zero;
total_energy_potential(:) = zero;

%%
%---
%---  beginning of time loop
%---
it=0;
for temp=1:temp_step:nstep
it=it+1;
%------------------------------------------------------------
% compute stress sigma and update memory variables for C-PML
%------------------------------------------------------------
for j = 2:ny;
for i = 1:nx-1;
% interpolate material parameters at the right location in the staggered grid cell
lambda_half_x = 0.5d0 .*(lambda(i+1,j) + lambda(i,j));
mu_half_x = 0.5d0 .*(mu(i+1,j) + mu(i,j));
lambda_plus_two_mu_half_x = lambda_half_x + 2.0d0 .* mu_half_x;
value_dvx_dx =(vx(i+1,j) - vx(i,j)) ./ deltax;
value_dvy_dy =(vy(i,j) - vy(i,j-1)) ./ deltay;
memory_dvx_dx(i,j) = b_x_half(i) .* memory_dvx_dx(i,j) + a_x_half(i) .* value_dvx_dx;
memory_dvy_dy(i,j) = b_y(j) .* memory_dvy_dy(i,j) + a_y(j) .* value_dvy_dy;
value_dvx_dx = value_dvx_dx ./ k_x_half(i) + memory_dvx_dx(i,j);
value_dvy_dy = value_dvy_dy ./ k_y(j) + memory_dvy_dy(i,j);
sigmaxx(i,j) = sigmaxx(i,j) +(lambda_plus_two_mu_half_x .* value_dvx_dx + lambda_half_x .* value_dvy_dy) .* deltat;
sigmayy(i,j) = sigmayy(i,j) +(lambda_half_x .* value_dvx_dx + lambda_plus_two_mu_half_x .* value_dvy_dy) .* deltat;
end; i =fix(nx-1+1);
end; j =fix(ny+1);
for j = 1:ny-1;
for i = 2:nx;
% interpolate material parameters at the right location in the staggered grid cell
mu_half_y = 0.5d0 .*(mu(i,j+1) + mu(i,j));
value_dvy_dx =(vy(i,j) - vy(i-1,j)) ./ deltax;
value_dvx_dy =(vx(i,j+1) - vx(i,j)) ./ deltay;
memory_dvy_dx(i,j) = b_x(i) .* memory_dvy_dx(i,j) + a_x(i) .* value_dvy_dx;
memory_dvx_dy(i,j) = b_y_half(j) .* memory_dvx_dy(i,j) + a_y_half(j) .* value_dvx_dy;
value_dvy_dx = value_dvy_dx ./ k_x(i) + memory_dvy_dx(i,j);
value_dvx_dy = value_dvx_dy ./ k_y_half(j) + memory_dvx_dy(i,j);
sigmaxy(i,j) = sigmaxy(i,j) + mu_half_y .*(value_dvy_dx + value_dvx_dy) .* deltat;
end; i =fix(nx+1);
end; j =fix(ny-1+1);
%--------------------------------------------------------
% compute velocity and update memory variables for C-PML
%--------------------------------------------------------
for j = 2:ny;
for i = 2:nx;
value_dsigmaxx_dx =(sigmaxx(i,j) - sigmaxx(i-1,j)) ./ deltax;
value_dsigmaxy_dy =(sigmaxy(i,j) - sigmaxy(i,j-1)) ./ deltay;
memory_dsigmaxx_dx(i,j) = b_x(i) .* memory_dsigmaxx_dx(i,j) + a_x(i) .* value_dsigmaxx_dx;
memory_dsigmaxy_dy(i,j) = b_y(j) .* memory_dsigmaxy_dy(i,j) + a_y(j) .* value_dsigmaxy_dy;
value_dsigmaxx_dx = value_dsigmaxx_dx ./ k_x(i) + memory_dsigmaxx_dx(i,j);
value_dsigmaxy_dy = value_dsigmaxy_dy ./ k_y(j) + memory_dsigmaxy_dy(i,j);
vx(i,j) = vx(i,j) +(value_dsigmaxx_dx + value_dsigmaxy_dy) .* deltat ./ rho(i,j);
end; i =fix(nx+1);
end; j =fix(ny+1);
for j = 1:ny-1;
for i = 1:nx-1;
% interpolate density at the right location in the staggered grid cell
rho_half_x_half_y = 0.25d0 .*(rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1));
value_dsigmaxy_dx =(sigmaxy(i+1,j) - sigmaxy(i,j)) ./ deltax;
value_dsigmayy_dy =(sigmayy(i,j+1) - sigmayy(i,j)) ./ deltay;
memory_dsigmaxy_dx(i,j) = b_x_half(i) .* memory_dsigmaxy_dx(i,j) + a_x_half(i) .* value_dsigmaxy_dx;
memory_dsigmayy_dy(i,j) = b_y_half(j) .* memory_dsigmayy_dy(i,j) + a_y_half(j) .* value_dsigmayy_dy;
value_dsigmaxy_dx = value_dsigmaxy_dx ./ k_x_half(i) + memory_dsigmaxy_dx(i,j);
value_dsigmayy_dy = value_dsigmayy_dy ./ k_y_half(j) + memory_dsigmayy_dy(i,j);
vy(i,j) = vy(i,j) +(value_dsigmaxy_dx + value_dsigmayy_dy) .* deltat ./ rho_half_x_half_y;
end; i =fix(nx-1+1);
end; j =fix(ny-1+1);
% add the source (force vector located at a given grid point)
a = pi.*pi.*f0.*f0;
t = (it-1).*deltat;

% Gaussian
if sourcetype=='1'
source_term = factor * exp(-a*(t-t0).^2);
end
% first derivative of a Gaussian
if sourcetype=='2'
source_term = - factor .* 2.0d0.*a.*(t-t0).*exp(-a.*(t-t0).^2);
end
% Ricker source time function (second derivative of a Gaussian)
if sourcetype=='3'
source_term = factor * (1.0d0 - 2.0d0*a*(t-t0).^2)*exp(-a*(t-t0).^2);
end

force_x = sin(angle_force .* degrees_to_radians) .* source_term;
force_y = cos(angle_force .* degrees_to_radians) .* source_term;
% define location of the source
i = fix(isource);
j = fix(jsource);
% interpolate density at the right location in the staggered grid cell
rho_half_x_half_y = 0.25d0 .*(rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1));
vx(i,j) = vx(i,j) + force_x .* deltat ./ rho(i,j);
vy(i,j) = vy(i,j) + force_y .* deltat ./ rho_half_x_half_y;
% Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
vx(1,:) = zero;
vx(nx,:) = zero;
vx(:,1) = zero;
vx(:,ny) = zero;
vy(1,:) = zero;
vy(nx,:) = zero;
vy(:,1) = zero;
vy(:,ny) = zero;
% store seismograms
for irec = 1:nrec;
sisvx(it+1,irec) = vx(ix_rec(irec),iy_rec(irec));
sisvy(it+1,irec) = vy(ix_rec(irec),iy_rec(irec));
% solid stress xx
sissigmaxx(it+1,irec) = sigmaxx(ix_rec(irec),iy_rec(irec));
% solid stress xy
sissigmaxy(it+1,irec) = sigmaxy(ix_rec(irec),iy_rec(irec));
% solid stress yy
sissigmayy(it+1,irec) = sigmayy(ix_rec(irec),iy_rec(irec));

end; irec =fix(nrec+1);
% compute total energy in the medium (without the PML layers)
% compute kinetic energy first, defined as 1/2 rho ||v||^2
% in principle we should use rho_half_x_half_y instead of rho for vy
% in order to interpolate density at the right location in the staggered grid cell
% but in a homogeneous medium we can safely ignore it
total_energy_kinetic(it+1) = 0.5d0 .* sum(sum(rho([npoints_pml+1:nx-npoints_pml],[npoints_pml+1:ny-npoints_pml]).*(vx([npoints_pml+1:nx-npoints_pml],[npoints_pml+1:ny-npoints_pml]).^2 +vy([npoints_pml+1:nx-npoints_pml],[npoints_pml+1:ny-npoints_pml]).^2)));
% add potential energy, defined as 1/2 epsilon_ij sigma_ij
% in principle we should interpolate the medium parameters at the right location
% in the staggered grid cell but in a homogeneous medium we can safely ignore it
total_energy_potential(it+1) = zero;
for j = npoints_pml+1: ny-npoints_pml;
for i = npoints_pml+1: nx-npoints_pml;
epsilon_xx =((lambda(i,j) + 2.0d0.*mu(i,j)) .* sigmaxx(i,j) - lambda(i,j) .*sigmayy(i,j)) ./(4.0d0 .* mu(i,j) .*(lambda(i,j) + mu(i,j)));
epsilon_yy =((lambda(i,j) + 2.0d0.*mu(i,j)) .* sigmayy(i,j) - lambda(i,j) .*sigmaxx(i,j)) ./(4.0d0 .* mu(i,j) .*(lambda(i,j) + mu(i,j)));
epsilon_xy = sigmaxy(i,j) ./(2.0d0 .* mu(i,j));
total_energy_potential(it+1) = total_energy_potential(it+1) +0.5d0 .*(epsilon_xx .* sigmaxx(i,j) + epsilon_yy .* sigmayy(i,j) + 2.0d0 .* epsilon_xy .* sigmaxy(i,j));
end; i = fix(nx-npoints_pml+1);
end; j = fix(ny-npoints_pml+1);
% output information
if(rem(it-1,it_display) == 0)
% print maximum of norm of velocity
tempmaxval=sqrt(vx.^2 + vy.^2);
velocnorm = max(tempmaxval(:));
fprintf(['%s %0.15g \n'], 'Time step # ',it);
fprintf(['%s %0.15g %s \n'], 'Time: ',((it-1).*deltat),' seconds');
fprintf(['%s %0.15g \n'], 'Max norm velocity vector V (m/s) = ',velocnorm);
fprintf(['%s %0.15g \n'], 'total energy = ',total_energy_kinetic(it+1) + total_energy_potential(it+1));
fprintf(['%0.15g \n']);
% check stability of the code, exit if unstable
if(velocnorm > stability_threshold)
error(['stop encountered in original Matlab code  ',char(10),' code became unstable and blew up;']);
end;

%imagesc(vx)
%surf(vx)
%plot(vx(:,fix(ny/2)))
%frame=getframe();


end;
% end of time loop
end; it =fix(nstep+1);

%%

%save the vx, vy and pf and plot them

%velocity X
    
vx_file=zeros(nstep,2);
fid_20=fopen('./vx1.dat','w+');
it=0;
for temp=1:temp_step:nstep;
it=it+1;
vx_file(it,1)= ((it-1).*deltat - t0);
vx_file(it,2)= sisvx(it,1);
fprintf(fid_20,['%0.15g %0.15g \n'],vx_file(it,1) , vx_file(it,2));
end; it=fix(nstep+1);
fclose(fid_20);

subplot(4,1,1);
plot(vx_file(:,1),vx_file(:,2));
xlabel('Time( in seconds)');
ylabel({'Velocity X','( in mm/s)'});
title(' Velocity X vs. Time AT SOURCE');
%saveas(gcf,'./Vx1.jpg');
%saveas(gca,'./Vx1.jpg');

 %%%%%%%%%%%%%%%%%%
%velocity X
    
vx_file=zeros(nstep,2);
fid_20=fopen('./vx2.dat','w+');
it=0;
for temp=1:temp_step:nstep;
it=it+1;
vx_file(it,1)= ((it-1).*deltat - t0);
vx_file(it,2)= sisvx(it,2);
fprintf(fid_20,['%0.15g %0.15g \n'],vx_file(it,1) , vx_file(it,2));
end; it=fix(nstep+1);
fclose(fid_20);


subplot(4,1,2);
plot(vx_file(:,1),vx_file(:,2));
xlabel('Time( in seconds)');
ylabel({'Velocity X','( in mm/s)'});
title(' Velocity X vs. Time AT EDGE');
%saveas(gcf,'./Vx2.jpg');
%saveas(gca,'./Vx2.jpg');

%%%%%%%%%%%%%%%%%%%%55
%velocity Y
vy_file=zeros(nstep,2);
fid_20=fopen('./vy1.dat','w+');
it=0;
for temp=1:temp_step:nstep;
it=it+1;
vy_file(it,1)= ((it-1).*deltat - t0);
vy_file(it,2)= sisvx(it,1);
fprintf(fid_20,['%0.15g %0.15g \n'],vy_file(it,1) , vy_file(it,2));
end; it=fix(nstep+1);
fclose(fid_20);

%plot it
subplot(4,1,3);
plot(vy_file(:,1),vy_file(:,2));
xlabel('Time( in seconds)');
ylabel({'Velocity Y','( in mm/s)'});
title('Velocity Y vs. Time AT SOURCE');
%saveas(gcf,'./Vy1.jpg');
%saveas(gca,'./Vy1.jpg');

%%%%%%%%%%%%%%%%%%%%%55
%velocity Y
vy_file=zeros(nstep,2);
fid_20=fopen('./vy2.dat','w+');
it=0;
for temp=1:temp_step:nstep;
it=it+1;
vy_file(it,1)= ((it-1).*deltat - t0);
vy_file(it,2)= sisvx(it,2);
fprintf(fid_20,['%0.15g %0.15g \n'],vy_file(it,1) , vy_file(it,2));
end; it =fix(nstep+1);
fclose(fid_20);

%plot it
subplot(4,1,4);
plot(vy_file(:,1),vy_file(:,2));
xlabel('Time( in seconds)');
ylabel({'Velocity Y','( in mm/s)'});
title('Velocity Y vs. Time AT EDGE');
%saveas(gcf,'./Vy2.jpg');
%saveas(gca,'./Vy2.jpg');
%%%%%%%%%%%%%%%%%%%%%%%%%


filename1= sprintf('2D_elastic_Velocity_type=%s.jpg',sourcetype)
filename2= sprintf('2D_elastic_Velocity_type=%s.fig',sourcetype)
saveas(gcf,filename1);
saveas(gca,filename1);
saveas(gcf,filename2);
saveas(gca,filename2);

%%%%%%%%%%%%%%%%%%%%%5555

% SOlid Stress

stressxx_file=zeros(nstep,2);
fid_20=fopen('./stressxx1.dat','w+');
for it=1:nstep;
stressxx_file(it,1)= ((it-1).*deltat - t0);
stressxx_file(it,2)= sissigmaxx(it,1);
%fprintf(fid_20,['%0.15g %0.15g \n'],pf_file(it,1) , pf_file(it,2));
end; it =fix(nstep+1);
fclose(fid_20);

%plot it
subplot(3,2,1);
plot(stressxx_file(:,1),stressxx_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Solid Stress XX','(in Pa)'});
title(' Solid Stress XX vs. Time AT SOURCE');
%saveas(gcf,'./pf2.jpg');
%saveas(gca,'./pf2.jpg');

%%%%%%%%%%%%%%%%%%55
% SOlid Stress

stressxx_file=zeros(nstep,2);
fid_20=fopen('./stressxx2.dat','w+');
for it=1:nstep;
stressxx_file(it,1)= ((it-1).*deltat - t0);
stressxx_file(it,2)= sissigmaxx(it,2);
%fprintf(fid_20,['%0.15g %0.15g \n'],pf_file(it,1) , pf_file(it,2));
end; it =fix(nstep+1);
fclose(fid_20);

%plot it
subplot(3,2,2);
plot(stressxx_file(:,1),stressxx_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Solid Stress XX','(in Pa)'});
title(' Solid Stress XX vs. Time AT EDGE');
%saveas(gcf,'./pf2.jpg');
%saveas(gca,'./pf2.jpg');

%%%%%%%%%%%%%%%%%%55

%filename1= sprintf('2D_elastic_StressXX_type=%s.jpg',sourcetype)
%filename2= sprintf('2D_elastic_StressXX_type=%s.fig',sourcetype)
%saveas(gcf,filename1);
%saveas(gca,filename1);
%saveas(gcf,filename2);
%saveas(gca,filename2);


%%%%%%%%%%%%%%%%%%%%55
% SOlid Stress

stressxy_file=zeros(nstep,2);
fid_20=fopen('./stressxy1.dat','w+');
for it=1:nstep;
stressxy_file(it,1)= ((it-1).*deltat - t0);
stressxy_file(it,2)= sissigmaxy(it,1);
%fprintf(fid_20,['%0.15g %0.15g \n'],pf_file(it,1) , pf_file(it,2));
end; it =fix(nstep+1);
fclose(fid_20);

%plot it
%subplot(2,1,1);
subplot(3,2,3);
plot(stressxy_file(:,1),stressxy_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Solid Stress XY','(in Pa)'});
title(' Solid Stress XY vs. Time AT SOURCE');
%saveas(gcf,'./pf2.jpg');
%saveas(gca,'./pf2.jpg');

%%%%%%%%%%%%%%%%%%55

% SOlid Stress

stressxy_file=zeros(nstep,2);
fid_20=fopen('./stressxy2.dat','w+');
for it=1:nstep;
stressxy_file(it,1)= ((it-1).*deltat - t0);
stressxy_file(it,2)= sissigmaxy(it,2);
%fprintf(fid_20,['%0.15g %0.15g \n'],pf_file(it,1) , pf_file(it,2));
end; it =fix(nstep+1);
fclose(fid_20);

%plot it
%subplot(2,1,2);
subplot(3,2,4);
plot(stressxy_file(:,1),stressxy_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Solid Stress XY','(in Pa)'});
title(' Solid Stress XY vs. Time AT EDGE');
%saveas(gcf,'./pf2.jpg');
%saveas(gca,'./pf2.jpg');

%%%%%%%%%%%%%%%%%%55

%filename1= sprintf('2D_elastic_StressXY_type=%s.jpg',sourcetype)
%filename2= sprintf('2D_elastic_StressXY_type=%s.fig',sourcetype)
%saveas(gcf,filename1);
%saveas(gca,filename1);
%saveas(gcf,filename2);
%saveas(gca,filename2);


%%%%%%%%%%%%%%%%%%%%55
% SOlid Stress

stressyy_file=zeros(nstep,2);
fid_20=fopen('./stressyy1.dat','w+');
for it=1:nstep;
stressyy_file(it,1)= ((it-1).*deltat - t0);
stressyy_file(it,2)= sissigmayy(it,1);
%fprintf(fid_20,['%0.15g %0.15g \n'],pf_file(it,1) , pf_file(it,2));
end; it =fix(nstep+1);
fclose(fid_20);

%plot it
%subplot(2,1,1);
subplot(3,2,5);
plot(stressyy_file(:,1),stressyy_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Solid Stress YY','(in Pa)'});
title(' Solid Stress YY vs. Time AT SOURCE');
%saveas(gcf,'./pf2.jpg');
%saveas(gca,'./pf2.jpg');

%%%%%%%%%%%%%%%%%%55


% SOlid Stress

stressyy_file=zeros(nstep,2);
fid_20=fopen('./stressyy2.dat','w+');
for it=1:nstep;
stressyy_file(it,1)= ((it-1).*deltat - t0);
stressyy_file(it,2)= sissigmayy(it,2);
%fprintf(fid_20,['%0.15g %0.15g \n'],pf_file(it,1) , pf_file(it,2));
end; it =fix(nstep+1);
fclose(fid_20);

%plot it
%subplot(2,1,2);
subplot(3,2,6);
plot(stressyy_file(:,1),stressyy_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Solid Stress YY','(in Pa)'});
title(' Solid Stress YY vs. Time AT EDGE');
%saveas(gcf,'./pf2.jpg');
%saveas(gca,'./pf2.jpg');

%%%%%%%%%%%%%%%%%%55

%filename1= sprintf('2D_elastic_StressYY_type=%s.jpg',sourcetype)
%filename2= sprintf('2D_elastic_StressYY_type=%s.fig',sourcetype)
%saveas(gcf,filename1);
%saveas(gca,filename1);
%saveas(gcf,filename2);
%saveas(gca,filename2);


%%%%%%%%%%%%%%%%%%%%55

filename1= sprintf('2D_elastic_Stress_type=%s.jpg',sourcetype)
filename2= sprintf('2D_elastic_Stress_type=%s.fig',sourcetype)
saveas(gcf,filename1);
saveas(gca,filename1);
saveas(gcf,filename2);
saveas(gca,filename2);


%%%%%%%%%%%%%%%%%%%%55
%  total energy
%energy_file=zeros(nstep,2) ;
%for it = 1:nstep;  
%energy_file(it,1)=((it-1).*deltat);   
%energy_file(it,2)=(total_energy_kinetic(it) + total_energy_potential(it));
%end; it =fix(nstep+1);

% savemlv ENergy
%fid_20=fopen('./energy.dat','w+');
%for it = 1:nstep;
%fprintf(fid_20,['%0.15g %0.15g \n'],energy_file(it,1) , energy_file(it,2));
%end; it =fix(nstep+1);
%fclose(fid_20);

%plot it
%subplot(3,1,3);
%plot(energy_file(:,1),energy_file(:,2));
%xlabel('time( in seconds)');
%ylabel('Energy ( in JOules)');
%title(' plot of Energy versus time');
%saveas(gcf,'./Energy.jpg');
%saveas(gca,'./Energy.jpg');









clear tempmaxval fid_20 ans

clear a a_x a_x_half a_y a_y_half abscissa_in_pml abscissa_normalized alpha_max_pml alpha_prime_x alpha_prime_x_half alpha_prime_y alpha_prime_y_half angle_force b_x b_x_half b_y b_y_half courant_number cp cs d0_x d0_y d_x d_x_half d_y d_y_half degrees_to_radians deltat deltax deltay density dist distval epsilon_xx epsilon_xy epsilon_yy f0 factor force_x force_y hugeval i irec isource it it_display ix_rec iy_rec j jsource k_max_pml k_x k_x_half k_y k_y_half lambda lambda_half_x lambda_plus_two_mu_half_x memory_dsigmaxx_dx memory_dsigmaxy_dx memory_dsigmaxy_dy memory_dsigmayy_dy memory_dvx_dx memory_dvx_dy memory_dvy_dx memory_dvy_dy mu mu_half_x mu_half_y npoints_pml npower nrec nstep nx ny pi rcoef rho rho_half_x_half_y sigmaxx sigmaxy sigmayy sisvx sisvy source_term stability_threshold t t0 thickness_pml_x thickness_pml_y total_energy_kinetic total_energy_potential use_pml_xmax use_pml_xmin use_pml_ymax use_pml_ymin value_dsigmaxx_dx value_dsigmaxy_dx value_dsigmaxy_dy value_dsigmayy_dy value_dvx_dx value_dvx_dy value_dvy_dx value_dvy_dy velocnorm vx vy xdeb xfin xoriginleft xoriginright xrec xsource xspacerec xval ydeb yfin yoriginbottom yorigintop yrec ysource yspacerec yval zero ; 


fprintf(['%0.15g \n']);
fprintf(['%s \n'], 'End of the simulation');
fprintf(['%0.15g \n']);
%----
%----  savemlv the seismograms in ASCII text format
%----
