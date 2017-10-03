
clear global; clear functions;
clear all;
global a a1_b a1_t a_outside a_inside a_x a_x_half_x a_y a_y_half_y abscissa_in_pml abscissa_normalized alp_sigma2 alpha alpha_outside alpha_max_pml alpha_prime_x alpha_prime_x_half_x alpha_prime_y alpha_prime_y_half_y alpha_inside angle_force b1_b b1_t b_x b_x_half_x b_y b_y_half_y c1 c1_b c1_t c2 c33_half_y co courant_number_outside courant_number_inside cp_outside cp_inside cps_outside cps_inside cs_outside cs_inside d0_x d0_y d_x d_x_half_x d_y d_y_half_y degrees_to_radians delta_b delta_t deltat deltax deltay dispersion_number_outside dispersion_number_inside dist distval epsilon_xx epsilon_xy epsilon_yy etaokappa etaokappa_outside etaokappa_half_x_half_y etaokappa_inside f0 factor ga_b ga_t gamma11 gamma12_1 gamma22 heterogeneous_model hugeval i ideb interface_height irec isource it it_display ix_rec iy_rec j jdeb jinterface jsource k_max_pml k_x k_x_half_x k_y k_y_half_y lambda_b lambda_t max_amplitude memory_dx_sigma2vx memory_dx_sigma2vxf memory_dx_sigmaxx memory_dx_sigmaxy memory_dx_sigmayy memory_dx_vx1 memory_dx_vx2 memory_dx_vy memory_dy_sigma2vy memory_dy_sigma2vyf memory_dy_sigmaxx memory_dy_sigmaxy memory_dy_sigmayy memory_dy_vx memory_dy_vy1 memory_dy_vy2 npoints_pml npower nrec nstep nx ny phi_outside phi_inside pi r_b r_t rbm rbm_outside rbm_inside rcoef rho rho_outside rho_half_x_half_y rho_inside rhof rhof_outside rhof_half_x_half_y rhof_inside rhos_outside rhos_inside rlambdac rlambdac_outside rlambdac_inside rlambdao rlambdao_outside rlambdao_inside rmu rmu_outside rmu_inside ro11_b ro11_t ro12_b ro12_t ro22_b ro22_t rsm rsm_outside rsm_half_x_half_y rsm_inside s_b s_t sigma2 sigmaxx sigmaxy sigmayy sisp sisvx sisvy source_term stability_threshold t t0 thickness_pml_x thickness_pml_y total_energy_kinetic total_energy_potential use_pml_outside use_pml_left use_pml_right use_pml_inside value_dx_sigma2vxf value_dx_sigmaxx value_dx_sigmaxy value_dx_vx1 value_dx_vx2 value_dx_vy value_dy_sigma2vyf value_dy_sigmaxx value_dy_sigmaxy value_dy_vx value_dy_vy1 value_dy_vy2 velocnorm_all vnorm vtemp vx vxf vy vyf xdeb xfin xi_1 xi_2 xoriginleft xoriginright xrec xsource xspacerec xval ydeb yfin yoriginoutside yorigininside yrec ysource yspacerec yval zero ; 



sourcetype= input(' [1] Gaussian   [2] First Gaussian derivative [3] Second Gaussian Derivative : ' ,'s')


%%
if isempty(nx), nx = 100; end;
if isempty(ny), ny = 100; end;
% size of a grid cell
if isempty(deltax), deltax = 1d0; end;
if isempty(deltay), deltay = deltax; end;

radius=30;

% flags to add PML layers to the edges of the grid
if isempty(use_pml_left), use_pml_left   = true; end;
if isempty(use_pml_right), use_pml_right  = true; end;
if isempty(use_pml_outside), use_pml_outside = true; end;
if isempty(use_pml_inside), use_pml_inside    = true; end;
% thickness of the PML layer in grid points
if isempty(npoints_pml), npoints_pml = 10; end;

%%
% total number of time steps
%if isempty(nstep), nstep = 100; end;

if isempty(nstep), nstep = 5000; end;
% time step in seconds
if isempty(deltat), deltat = 1.0d-07; end;
% parameters for the source
if isempty(f0), f0 = 80000.0d0; end;
if isempty(t0), t0 = 1.0d0./f0; end;
if isempty(factor), factor =1.0d02; end;
% source
if isempty(isource), isource = fix(nx./2)+1; end;
if isempty(jsource), jsource = fix(ny./2) +1; end;
if isempty(ideb), ideb =  fix(nx ./ 2) + 1; end;
if isempty(jdeb), jdeb =  fix(ny ./ 2) + 1; end;
if isempty(xsource), xsource = deltax .* isource; end;
if isempty(ysource), ysource = deltay .* jsource; end;
% angle of source force clockwise with respect to vertical (Y) axis
if isempty(angle_force), angle_force = 90.0d0; end;
% receivers
if isempty(nrec), nrec = 2; end;

% first receiver y in meters
if isempty(xdeb), xdeb =xsource; end;
% first receiver deb in meters
if isempty(ydeb), ydeb = ysource; end;

% second receiver y in meters
if isempty(xfin), xfin =deltax .* (isource-radius); end;
% second receiver fin in meters
if isempty(yfin), yfin = ysource; end;

% display information on the screen from time to time
if isempty(it_display), it_display = 5; end;

%%
% heterogeneous model and height of the interface
if isempty(interface_height), interface_height =105.0d0+npoints_pml.*deltay; end;
if isempty(jinterface), jinterface=fix(interface_height./deltay)+1; end;
if isempty(co), co=0; end;
if isempty(c1), c1=0; end;
if isempty(c2), c2=0; end;
if isempty(vtemp), vtemp=0; end;


% model mud saturated with water, see article by Martin and Komatitsch
if isempty(etaokappa_outside), etaokappa_outside=0.000001*3.38d05; end;
if isempty(rmu_outside), rmu_outside = 5.25d09; end;
if isempty(phi_outside), phi_outside =0.25d0; end;
if isempty(a_outside), a_outside = 2.49d0; end;
if isempty(rhos_outside), rhos_outside = 0.000001*2588.0d0; end;
if isempty(rhof_outside), rhof_outside = 0.000001*952.4d0; end;
if isempty(rho_outside), rho_outside = 0.000001*2179.1d0; end;
if isempty(rsm_outside), rsm_outside =0.000001*9486.0d0; end;
if isempty(alpha_outside), alpha_outside=0.89d0; end;
if isempty(rbm_outside), rbm_outside =7.71d09; end;
if isempty(rlambdao_outside), rlambdao_outside = 6.2d08; end;
if isempty(rlambdac_outside), rlambdac_outside =rlambdao_outside+alpha_outside.^2.*rbm_outside; end;
if isempty(ro11_b), ro11_b=rho_outside+phi_outside.*rhof_outside.*(a_outside-2.0d0); end;
if isempty(ro12_b), ro12_b=phi_outside.*rhof_outside.*(1.0d0-a_outside); end;
if isempty(ro22_b), ro22_b=a_outside.*phi_outside.*rhof_outside; end;
if isempty(lambda_b), lambda_b=rlambdao_outside+rbm_outside.*(alpha_outside-phi_outside).^2; end;
if isempty(r_b), r_b=rbm_outside.*phi_outside.^2; end;
if isempty(ga_b), ga_b=rbm_outside.*phi_outside.*(alpha_outside-phi_outside); end;
if isempty(s_b), s_b=lambda_b+2.*rmu_outside; end;
if isempty(c1_b), c1_b=s_b.*r_b-ga_b.^2; end;
if isempty(b1_b), b1_b=-s_b.*ro22_b-r_b.*ro11_b+2.*ga_b.*ro12_b; end;
if isempty(a1_b), a1_b=ro11_b.*ro22_b-ro12_b.^2; end;
if isempty(delta_b), delta_b=b1_b.^2-4.0d0.*a1_b.*c1_b; end;
if isempty(cp_outside), cp_outside=0; end;
if isempty(cps_outside), cps_outside=0; end;
if isempty(cs_outside), cs_outside=0; end;


if isempty(etaokappa_inside), etaokappa_inside=0.000001*3.33d06; end;
if isempty(rmu_inside), rmu_inside = 2.4d09; end;
if isempty(phi_inside), phi_inside =0.1d0; end;
if isempty(a_inside), a_inside = 2.42d0; end;
if isempty(rhos_inside), rhos_inside = 0.000001*2250.0d0; end;
if isempty(rhof_inside), rhof_inside = 0.000001*1040.0d0; end;
if isempty(rho_inside), rho_inside = 0.000001*2129.0d0; end;
if isempty(rsm_inside), rsm_inside =0.000001*25168.0d0; end;
if isempty(alpha_inside), alpha_inside=0.58d0; end;
if isempty(rbm_inside), rbm_inside = 7.34d09; end;
if isempty(rlambdao_inside), rlambdao_inside =6.0d08; end;
if isempty(rlambdac_inside), rlambdac_inside =rlambdao_inside+alpha_inside.^2.*rbm_inside;  end;
if isempty(ro11_t), ro11_t=rho_inside+phi_inside.*rhof_inside.*(a_inside-2.0d0); end;
if isempty(ro12_t), ro12_t=phi_inside.*rhof_inside.*(1.0d0-a_inside); end;
if isempty(ro22_t), ro22_t=a_inside.*phi_inside.*rhof_inside; end;
if isempty(lambda_t), lambda_t=rlambdao_inside+rbm_inside.*(alpha_inside-phi_inside).^2; end;
if isempty(r_t), r_t=rbm_inside.*phi_inside.^2; end;
if isempty(ga_t), ga_t=rbm_inside.*phi_inside.*(alpha_inside-phi_inside); end;
if isempty(s_t), s_t=lambda_t+2.*rmu_inside; end;
if isempty(c1_t), c1_t=s_t.*r_t-ga_t.^2; end;
if isempty(b1_t), b1_t=-s_t.*ro22_t-r_t.*ro11_t+2.*ga_t.*ro12_t; end;
if isempty(a1_t), a1_t=ro11_t.*ro22_t-ro12_t.^2; end;
if isempty(delta_t), delta_t=b1_t.^2-4.0d0.*a1_t.*c1_t; end;
if isempty(cp_inside), cp_inside=0; end;
if isempty(cps_inside), cps_inside=0; end;
if isempty(cs_inside), cs_inside=0; end;

fprintf('--------------------------------');
fprintf('rlambdac_inside =  %0.15g \n', rlambdac_inside);
fprintf('ro11_t =  %0.15g \n', ro11_t); 
fprintf('ro12_t =  %0.15g \n', ro12_t );
fprintf('ro22_t =  %0.15g \n', ro22_t );
 fprintf('lambda_t =  %0.15g \n', lambda_t ) ;
 fprintf('r_t =  %0.15g \n', r_t ) ;
 fprintf('ga_t =  %0.15g \n', ga_t ) ;
 fprintf('s_t =  %0.15g \n', s_t );
 fprintf('c1_t =  %0.15g \n ', c1_t ) ;
 fprintf('b1_t =  %0.15g \n', b1_t ) ;
 fprintf('a1_t =  %0.15g \n ', a1_t ) ;
 fprintf('delta_t =  %0.15g \n', delta_t ) ;
 fprintf('cp_inside =  %0.15g \n', cp_inside ) ; 
 fprintf('cps_inside =  %0.15g \n', cps_inside) ;
 fprintf('cs_inside =  %0.15g \n', cs_inside ); 
fprintf('------------------------');

%%

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
if isempty(vx), vx=zeros(nx+1+1,ny+1+1); end;
if isempty(vy), vy=zeros(nx+1+1,ny+1+1); end;
if isempty(sigmaxx), sigmaxx=zeros(nx+1+1,ny+1+1); end;
if isempty(sigma2), sigma2=zeros(nx+1+1,ny+1+1); end;
if isempty(alp_sigma2), alp_sigma2=zeros(nx+1+1,ny+1+1); end;
if isempty(sigmayy), sigmayy=zeros(nx+1+1,ny+1+1); end;
if isempty(sigmaxy), sigmaxy=zeros(nx+1+1,ny+1+1); end;
if isempty(vnorm), vnorm=zeros(nx+1+1,ny+1+1); end;
if isempty(vxf), vxf=zeros(nx+1+1,ny+1+1); end;
if isempty(vyf), vyf=zeros(nx+1+1,ny+1+1); end;
if isempty(rho), rho=zeros(nx+1+1,ny+1+1); end;
if isempty(rhof), rhof=zeros(nx+1+1,ny+1+1); end;
if isempty(rsm), rsm=zeros(nx+1+1,ny+1+1); end;
if isempty(rmu), rmu=zeros(nx+1+1,ny+1+1); end;
if isempty(rlambdac), rlambdac=zeros(nx+1+1,ny+1+1); end;
if isempty(rbm), rbm=zeros(nx+1+1,ny+1+1); end;
if isempty(alpha), alpha=zeros(nx+1+1,ny+1+1); end;
if isempty(etaokappa), etaokappa=zeros(nx+1+1,ny+1+1); end;
if isempty(rlambdao), rlambdao=zeros(nx+1+1,ny+1+1); end;
% to interpolate material parameters at the right location in the staggered grid cell
if isempty(rho_half_x_half_y), rho_half_x_half_y=0; end;
if isempty(rhof_half_x_half_y), rhof_half_x_half_y=0; end;
if isempty(rsm_half_x_half_y), rsm_half_x_half_y=0; end;
if isempty(etaokappa_half_x_half_y), etaokappa_half_x_half_y=0; end;
% for evolution of total energy in the medium
if isempty(epsilon_xx), epsilon_xx=0; end;
if isempty(epsilon_yy), epsilon_yy=0; end;
if isempty(epsilon_xy), epsilon_xy=0; end;
if isempty(total_energy_kinetic), total_energy_kinetic=zeros(1,nstep); end;
if isempty(total_energy_potential), total_energy_potential=zeros(1,nstep); end;
if isempty(c33_half_y), c33_half_y=0; end;
% power to compute d0 profile
if isempty(npower), npower = 2.0d0; end;
% doubleprecision, parameter :: K_MAX_PML = 7.0d0   ! from Gedney page 8.11
% from Gedney page 8.11
if isempty(k_max_pml), k_max_pml = 1.0d0; end;
% doubleprecision, parameter :: ALPHA_MAX_PML = 0.05d0   ! from Gedney page 8.22
% from festa and Vilotte
if isempty(alpha_max_pml), alpha_max_pml = 2.0d0.*pi.*(f0./2.0d0); end;
% 2D arrays for the memory variables
if isempty(gamma11), gamma11=zeros(nx+1+1,ny+1+1); end;
if isempty(gamma22), gamma22=zeros(nx+1+1,ny+1+1); end;
if isempty(gamma12_1), gamma12_1=zeros(nx+1+1,ny+1+1); end;
if isempty(xi_1), xi_1=zeros(nx+1+1,ny+1+1); end;
if isempty(xi_2), xi_2=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dx_vx1), memory_dx_vx1=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dx_vx2), memory_dx_vx2=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dy_vx), memory_dy_vx=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dx_vy), memory_dx_vy=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dy_vy1), memory_dy_vy1=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dy_vy2), memory_dy_vy2=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dx_sigmaxx), memory_dx_sigmaxx=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dx_sigmayy), memory_dx_sigmayy=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dx_sigmaxy), memory_dx_sigmaxy=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dx_sigma2vx), memory_dx_sigma2vx=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dx_sigma2vxf), memory_dx_sigma2vxf=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dy_sigma2vy), memory_dy_sigma2vy=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dy_sigma2vyf), memory_dy_sigma2vyf=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dy_sigmaxx), memory_dy_sigmaxx=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dy_sigmayy), memory_dy_sigmayy=zeros(nx+1+1,ny+1+1); end;
if isempty(memory_dy_sigmaxy), memory_dy_sigmaxy=zeros(nx+1+1,ny+1+1); end;
% 1D arrays for the damping profiles
if isempty(d_x), d_x=zeros(1,nx); end;
if isempty(k_x), k_x=zeros(1,nx); end;
if isempty(alpha_prime_x), alpha_prime_x=zeros(1,nx); end;
if isempty(a_x), a_x=zeros(1,nx); end;
if isempty(b_x), b_x=zeros(1,nx); end;
if isempty(d_x_half_x), d_x_half_x=zeros(1,nx); end;
if isempty(k_x_half_x), k_x_half_x=zeros(1,nx); end;
if isempty(alpha_prime_x_half_x), alpha_prime_x_half_x=zeros(1,nx); end;
if isempty(a_x_half_x), a_x_half_x=zeros(1,nx); end;
if isempty(b_x_half_x), b_x_half_x=zeros(1,nx); end;
if isempty(d_y), d_y=zeros(1,ny); end;
if isempty(k_y), k_y=zeros(1,ny); end;
if isempty(alpha_prime_y), alpha_prime_y=zeros(1,ny); end;
if isempty(a_y), a_y=zeros(1,ny); end;
if isempty(b_y), b_y=zeros(1,ny); end;
if isempty(d_y_half_y), d_y_half_y=zeros(1,ny); end;
if isempty(k_y_half_y), k_y_half_y=zeros(1,ny); end;
if isempty(alpha_prime_y_half_y), alpha_prime_y_half_y=zeros(1,ny); end;
if isempty(a_y_half_y), a_y_half_y=zeros(1,ny); end;
if isempty(b_y_half_y), b_y_half_y=zeros(1,ny); end;
if isempty(thickness_pml_x), thickness_pml_x=0; end;
if isempty(thickness_pml_y), thickness_pml_y=0; end;
if isempty(xoriginleft), xoriginleft=0; end;
if isempty(xoriginright), xoriginright=0; end;
if isempty(yoriginoutside), yoriginoutside=0; end;
if isempty(yorigininside), yorigininside=0; end;
if isempty(rcoef), rcoef=0; end;
if isempty(d0_x), d0_x=0; end;
if isempty(d0_y), d0_y=0; end;
if isempty(xval), xval=0; end;
if isempty(yval), yval=0; end;
if isempty(abscissa_in_pml), abscissa_in_pml=0; end;
if isempty(abscissa_normalized), abscissa_normalized=0; end;
if isempty(value_dx_vx1), value_dx_vx1=0; end;
if isempty(value_dx_vx2), value_dx_vx2=0; end;
if isempty(value_dx_vy), value_dx_vy=0; end;
if isempty(value_dx_sigmaxx), value_dx_sigmaxx=0; end;
if isempty(value_dx_sigmaxy), value_dx_sigmaxy=0; end;
if isempty(value_dy_vy1), value_dy_vy1=0; end;
if isempty(value_dy_vy2), value_dy_vy2=0; end;
if isempty(value_dy_vx), value_dy_vx=0; end;
if isempty(value_dy_sigmaxx), value_dy_sigmaxx=0; end;
if isempty(value_dy_sigmaxy), value_dy_sigmaxy=0; end;
if isempty(value_dx_sigma2vxf), value_dx_sigma2vxf=0; end;
if isempty(value_dy_sigma2vyf), value_dy_sigma2vyf=0; end;
% for the source
if isempty(a), a=0; end;
if isempty(t), t=0; end;
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
if isempty(sisvx), sisvx=zeros(nstep,nrec); end;
if isempty(sisvy), sisvy=zeros(nstep,nrec); end;
if isempty(sisp), sisp=zeros(nstep,nrec); end;
if isempty(i), i=0; end;
if isempty(j), j=0; end;
if isempty(it), it=0; end;
if isempty(irec), irec=0; end;
if isempty(courant_number_outside), courant_number_outside=0; end;
if isempty(courant_number_inside), courant_number_inside=0; end;
if isempty(velocnorm_all), velocnorm_all=0; end;
if isempty(max_amplitude), max_amplitude=0; end;
if isempty(dispersion_number_outside), dispersion_number_outside=0; end;
if isempty(dispersion_number_inside), dispersion_number_inside=0; end;

%%
%---
%--- program starts here
%---
cp_outside=(-b1_b+sqrt(delta_b))./(2.0d0.*a1_b);
cps_outside=(-b1_b-sqrt(delta_b))./(2.0d0.*a1_b);
cp_outside=sqrt(cp_outside);
cps_outside=sqrt(cps_outside);
cs_outside=sqrt(rmu_outside./(ro11_b-ro12_b.^2./ro22_b));
cp_inside=(-b1_t+sqrt(delta_t))./(2.0d0.*a1_t);
cps_inside=(-b1_t-sqrt(delta_t))./(2.0d0.*a1_t);
cp_inside=sqrt(cp_inside);
cps_inside=sqrt(cps_inside);
cs_inside=sqrt(rmu_inside./(ro11_t-ro12_t.^2./ro22_t));
fprintf(['%s %0.15g \n'], 'cp_outside= ',cp_outside);
fprintf(['%s %0.15g \n'], 'cps_outside=',cps_outside);
fprintf(['%s %0.15g \n'], 'cs_outside= ',cs_outside);
fprintf(['%s %0.15g \n'], 'cp_inside= ',cp_inside);
fprintf(['%s %0.15g \n'], 'cps_inside=',cps_inside);
fprintf(['%s %0.15g \n'], 'cs_inside= ',cs_inside);
fprintf(['%s %0.15g \n'], 'rho_outside= ',rho_outside);
fprintf(['%s %0.15g \n'], 'rsm_outside= ',rsm_outside);
fprintf(['%s %0.15g \n'], 'rho_inside= ',rho_inside);
fprintf(['%s %0.15g \n'], 'rsm_inside= ',rsm_inside);
fprintf(['%s %0.15g \n'], 'rmu_outside= ',rmu_outside);
fprintf(['%s %0.15g \n'], 'rlambdac_outside= ',rlambdac_outside);
fprintf(['%s %0.15g \n'], 'rlambdao_outside= ',rlambdao_outside);
fprintf(['%s %0.15g \n'], 'alpha_outside= ',alpha_outside);
fprintf(['%s %0.15g \n'], 'rbM_outside= ',rbm_outside);
fprintf(['%s %0.15g \n'], 'etaokappa_outside= ',etaokappa_outside);
fprintf(['%s %0.15g \n'], 'rmu_inside= ',rmu_inside);
fprintf(['%s %0.15g \n'], 'rlambdac_inside= ',rlambdac_inside);
fprintf(['%s %0.15g \n'], 'rlambdao_inside= ',rlambdao_inside);
fprintf(['%s %0.15g \n'], 'alpha_inside= ',alpha_inside);
fprintf(['%s %0.15g \n'], 'rbM_inside= ',rbm_inside);
fprintf(['%s %0.15g \n'], 'etaokappa_inside= ',etaokappa_inside);
fprintf(['%s %0.15g \n'],  'DELTAT CPML=', deltat);
fprintf(['%s \n'], '2D poroelastic finite-difference code in velocity and stress formulation with C-PML');
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
error(['sinside encountered in original fortran code  ',char(10),' NPOWER must be greater than 1;']);
end;
% compute d0 from INRIA report section 6.1
if(heterogeneous_model)
d0_x = -(npower + 1) .* max(cp_outside,cp_inside) .* log(rcoef) ./(2.0d0 .* thickness_pml_x);
d0_y = -(npower + 1) .* max(cp_outside,cp_inside) .* log(rcoef) ./(2.0d0 .* thickness_pml_y);
else;
d0_x = -(npower + 1) .* cp_outside .* log(rcoef) ./(2.0d0 .* thickness_pml_x);
d0_y = -(npower + 1) .* cp_outside .* log(rcoef) ./(2.0d0 .* thickness_pml_y);
end;
fprintf(['%s %0.15g \n'], 'd0_x = ',d0_x);
fprintf(['%s %0.15g \n'], 'd0_y = ',d0_y);
d_x(:) = zero;
d_x_half_x(:) = zero;
d_y(:) = zero;
d_y_half_y(:) = zero;
k_x(:) = 1.0d0;
k_x_half_x(:) = 1.0d0;
k_y(:) = 1.0d0;
k_y_half_y(:) = 1.0d0;
alpha_prime_x(:) = zero;
alpha_prime_x_half_x(:) = zero;
alpha_prime_y(:) = zero;
alpha_prime_y_half_y(:) = zero;
a_x(:) = zero;
a_x_half_x(:) = zero;
a_y(:) = zero;
a_y_half_y(:) = zero;
% origin of the PML layer (position of right edge minus thickness, in meters)
xoriginleft = thickness_pml_x;
xoriginright =(nx-1).*deltax - thickness_pml_x;
for i = 1:nx;
% abscissa of current grid point along the damping profile
xval = deltax .* (i-1);
%!!! ---------- left edge
if(use_pml_left)
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
d_x_half_x(i) = d0_x .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_x_half_x(i) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_x_half_x(i) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
end;
%!!! ---------- right edge
if(use_pml_right)
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
d_x_half_x(i) = d0_x .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_x_half_x(i) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_x_half_x(i) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
end;
% just in case, for -5 at the end
if(alpha_prime_x(i) < zero)
alpha_prime_x(i) = zero;
end;
if(alpha_prime_x_half_x(i) < zero)
alpha_prime_x_half_x(i) = zero;
end;
b_x(i) = exp(-(d_x(i) ./ k_x(i) + alpha_prime_x(i)) .* deltat);
b_x_half_x(i) = exp(-(d_x_half_x(i) ./ k_x_half_x(i) + alpha_prime_x_half_x(i)) .* deltat);
% this to avoid division by zero inside the PML
if(abs(d_x(i)) > 1.0d-6)
a_x(i) = d_x(i) .*(b_x(i) - 1.0d0) ./(k_x(i) .*(d_x(i) + k_x(i) .* alpha_prime_x(i)));
end;
if(abs(d_x_half_x(i)) > 1.0d-6)
a_x_half_x(i) = d_x_half_x(i).*(b_x_half_x(i) - 1.0d0) ./(k_x_half_x(i) .*(d_x_half_x(i) + k_x_half_x(i) .* alpha_prime_x_half_x(i)));
end;
end; i =fix(nx+1);
%!!!!!!!!!!!! added Y damping profile
% origin of the PML layer (position of right edge minus thickness, in meters)
yoriginoutside = thickness_pml_y;
yorigininside = ny.*deltay - thickness_pml_y;
for j = 1:ny;
% abscissa of current grid point along the damping profile
yval = deltay .* (j-1);
%!!! ---------- outside edge
if(use_pml_outside)
% define damping profile at the grid points
abscissa_in_pml = yoriginoutside - yval;
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_y;
d_y(j) = d0_y .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_y(j) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_y(j) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
% define damping profile at half the grid points
abscissa_in_pml = yoriginoutside -(yval + deltay./2.0d0);
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_y;
d_y_half_y(j) = d0_y .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_y_half_y(j) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_y_half_y(j) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
end;
%!!! ---------- inside edge
if(use_pml_inside)
% define damping profile at the grid points
abscissa_in_pml = yval - yorigininside;
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_y;
d_y(j) = d0_y .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_y(j) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_y(j) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
% define damping profile at half the grid points
abscissa_in_pml = yval + deltay./2.0d0 - yorigininside;
if(abscissa_in_pml >= zero)
abscissa_normalized = abscissa_in_pml ./ thickness_pml_y;
d_y_half_y(j) = d0_y .* abscissa_normalized.^npower;
% this taken from Gedney page 8.2
k_y_half_y(j) = 1.0d0 +(k_max_pml - 1.0d0) .* abscissa_normalized.^npower;
alpha_prime_y_half_y(j) = alpha_max_pml .*(1.0d0 - abscissa_normalized);
end;
end;
% just in case, for -5 at the end
%   if(alpha_prime_y(j) < ZERO) alpha_prime_y(j) = ZERO
%   if(alpha_prime_y_half_y(j) < ZERO) alpha_prime_y_half_y(j) = ZERO
b_y(j) = exp(-(d_y(j) ./ k_y(j) + alpha_prime_y(j)) .* deltat);
b_y_half_y(j) = exp(-(d_y_half_y(j) ./ k_y_half_y(j) + alpha_prime_y_half_y(j)) .* deltat);
% this to avoid division by zero inside the PML
if(abs(d_y(j)) > 1.0d-6)
a_y(j) = d_y(j) .*(b_y(j) - 1.0d0)./(k_y(j) .*(d_y(j) + k_y(j) .* alpha_prime_y(j)));
end;
if(abs(d_y_half_y(j)) > 1.0d-6)
a_y_half_y(j) = d_y_half_y(j).*(b_y_half_y(j) - 1.0d0) ./(k_y_half_y(j) .*(d_y_half_y(j) + k_y_half_y(j) .* alpha_prime_y_half_y(j)));
end;
end; j =fix(ny+1);
% compute the Lame parameters and density
for j = 0:ny+1;
for i = 0:nx+1;
  
rho(i+1,j+1)= rho_outside;
rhof(i+1,j+1) = rhof_outside;
rsm(i+1,j+1) = rsm_outside;
rmu(i+1,j+1)= rmu_outside;
rlambdac(i+1,j+1) = rlambdac_outside;
rbm(i+1,j+1) = rbm_outside;
alpha(i+1,j+1)=alpha_outside;
etaokappa(i+1,j+1)=etaokappa_outside;
rlambdao(i+1,j+1) = rlambdao_outside;      

if(((i-isource).^2+(j-jsource).^2)<=(radius.^2))

    
rho(i+1,j+1)= rho_inside;
rhof(i+1,j+1) = rhof_inside;
rsm(i+1,j+1) = rsm_inside;
rmu(i+1,j+1)= rmu_inside;
rlambdac(i+1,j+1) = rlambdac_inside;
rbm(i+1,j+1) = rbm_inside;
alpha(i+1,j+1)=alpha_inside;
etaokappa(i+1,j+1)=etaokappa_inside;
rlambdao(i+1,j+1) = rlambdao_inside;
  

end;



end; i =fix(nx+1+1);
end; j =fix(ny+1+1);
% print position of the source
fprintf(['%0.15g \n']);
fprintf(['%s \n'], 'Position of the source:');
fprintf(['%0.15g \n']);
fprintf(['%s %0.15g \n'], 'x = ',xsource);
fprintf(['%s %0.15g \n'], 'y = ',ysource);
fprintf(['%0.15g \n']);
% define location of receivers
fprintf(['%0.15g \n']);
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
courant_number_outside = cp_outside .* deltat ./ min(deltax,deltay);
dispersion_number_outside=min(cs_outside,cps_outside)./(2.5d0.*f0.*max(deltax,deltay));
fprintf(['%s %0.15g \n'], 'Courant number at the outside is ',courant_number_outside);
fprintf(['%s %0.15g \n'], 'Dispersion number at the outside is ',dispersion_number_outside);
fprintf(['%0.15g \n']);
if(courant_number_outside > 1.0d0./sqrt(2.0d0))
error(['sinside encountered in original fortran code  ',char(10),' time step is too large, simulation will be unstable;']);
end;
if(heterogeneous_model)
courant_number_inside = max(cp_inside,cp_outside) .* deltat ./ min(deltax,deltay);
dispersion_number_inside=min([cs_inside,cs_outside,cps_outside,cps_inside])./(2.5d0.*f0.*max(deltax,deltay));
fprintf(['%s %0.15g \n'], 'Courant number at the inside is ',courant_number_inside);
fprintf(['%0.15g \n']);
fprintf(['%s %0.15g \n'], 'Dispersion number at the inside is ',dispersion_number_inside);
if(courant_number_inside > 6.0d0./7.0d0./sqrt(2.0d0))
error(['sinside encountered in original fortran code  ',char(10),' time step is too large, simulation will be unstable;']);
end;
end;
% suppress old files
% call system('rm -f Vx_*.dat Vy_*.dat vect*.ps image*.pnm image*.gif')
% initialize arrays
vx(:,:) = zero;
vy(:,:) = zero;
sigmaxx(:,:) = zero;
sigmayy(:,:) = zero;
sigmaxy(:,:) = zero;
sigma2(:,:) = zero;
alp_sigma2(:,:) = zero;
gamma11(:,:)=0.0d0;
gamma22(:,:)=0.0d0;
gamma12_1(:,:)=0.0d0;
gamma12_1(:,:)=0.0d0;
xi_1(:,:)=0.0d0;
xi_2(:,:)=0.0d0;
vxf(:,:) = zero;
vyf(:,:) = zero;
memory_dx_vx1(:,:)=0.0d0;
memory_dx_vx2(:,:)=0.0d0;
memory_dy_vx(:,:)=0.0d0;
memory_dx_vy(:,:)=0.0d0;
memory_dy_vy1(:,:)=0.0d0;
memory_dy_vy2(:,:)=0.0d0;
memory_dx_sigmaxx(:,:)=0.0d0;
memory_dx_sigmayy(:,:)=0.0d0;
memory_dx_sigmaxy(:,:)=0.0d0;
memory_dx_sigma2vx(:,:)=0.0d0;
memory_dx_sigma2vxf(:,:)=0.0d0;
memory_dy_sigmaxx(:,:)=0.0d0;
memory_dy_sigmayy(:,:)=0.0d0;
memory_dy_sigmaxy(:,:)=0.0d0;
memory_dy_sigma2vy(:,:)=0.0d0;
memory_dy_sigma2vyf(:,:)=0.0d0;
% initialize seismograms
sisvx(:,:) = zero;
sisvy(:,:) = zero;
sisp(:,:) = zero;
% initialize total energy
total_energy_kinetic(:) = zero;
total_energy_potential(:) = zero;

%%
%---
%---  beginning of time loop
%---
for it = 1:nstep;
%----------------------
% compute stress sigma
%----------------------
%-----------------------------------
% update memory variables for C-PML
%-----------------------------------
for j = 2:ny;
for i = 1:nx-1;
%  memory of sigmaxx
value_dx_sigmaxx =(27.0d0.*vx(i+1+1,j+1)-27.0d0.*vx(i+1,j+1)-vx(i+2+1,j+1)+vx(i-1+1,j+1))./deltax./24.0d0;
value_dy_sigmaxx =(27.0d0.*vy(i+1,j+1)-27.0d0.*vy(i+1,j-1+1)-vy(i+1,j+1+1)+vy(i+1,j-2+1))./deltay./24.0d0;
memory_dx_sigmaxx(i+1,j+1) = b_x_half_x(i) .* memory_dx_sigmaxx(i+1,j+1) + a_x_half_x(i) .* value_dx_sigmaxx;
memory_dy_sigmaxx(i+1,j+1) = b_y(j) .* memory_dy_sigmaxx(i+1,j+1) + a_y(j) .* value_dy_sigmaxx;
gamma11(i+1,j+1) = gamma11(i+1,j+1)+deltat.*(value_dx_sigmaxx ./ k_x_half_x(i) + memory_dx_sigmaxx(i+1,j+1));
gamma22(i+1,j+1) = gamma22(i+1,j+1)+deltat.*(value_dy_sigmaxx ./ k_y(j) + memory_dy_sigmaxx(i+1,j+1));
% sigma2
value_dx_sigma2vxf=(27.0d0.*vxf(i+1+1,j+1)-27.0d0.* vxf(i+1,j+1)-vxf(i+2+1,j+1)+vxf(i-1+1,j+1)) ./ deltax./24.0d0;
value_dy_sigma2vyf=(27.0d0.*vyf(i+1,j+1)-27.0d0.*vyf(i+1,j-1+1)-vyf(i+1,j+1+1)+vyf(i+1,j-2+1)) ./ deltay./24.0d0;
memory_dx_sigma2vxf(i+1,j+1) = b_x_half_x(i) .* memory_dx_sigma2vxf(i+1,j+1) + a_x_half_x(i) .* value_dx_sigma2vxf;
memory_dy_sigma2vyf(i+1,j+1) = b_y(j) .* memory_dy_sigma2vyf(i+1,j+1) + a_y(j) .* value_dy_sigma2vyf;
xi_1(i+1,j+1) = xi_1(i+1,j+1) -(value_dx_sigma2vxf./ k_x_half_x(i) + memory_dx_sigma2vxf(i+1,j+1)).*deltat;
xi_2(i+1,j+1) = xi_2(i+1,j+1) -(value_dy_sigma2vyf./k_y(j)+memory_dy_sigma2vyf(i+1,j+1)).*deltat;
sigma2(i+1,j+1)=-alpha(i+1,j+1).*rbm(i+1,j+1).*(gamma11(i+1,j+1)+gamma22(i+1,j+1))+rbm(i+1,j+1).*(xi_1(i+1,j+1)+xi_2(i+1,j+1));
end; i =fix(nx-1+1);
end; j =fix(ny+1);
% add the source (point source located at a given grid point)
a = pi.*pi.*f0.*f0;
t = (it-1).*deltat;

% Gaussian
if sourcetype=='1'
source_term = factor .* exp(-a.*(t-t0).^2)./(-2.0d0.*a);
end

% first derivative of a Gaussian
if sourcetype=='2'
source_term =  factor * 2.0d0*a*(t-t0)*exp(-a*(t-t0).^2);
source_term =  factor *(t-t0)*exp(-a*(t-t0).^2);
end

% Ricker source time function (second derivative of a Gaussian)
if sourcetype=='3'
source_term = factor * (1.0d0 - 2.0d0*a*(t-t0).^2)*exp(-a*(t-t0).^2);
end



% define location of the source
i = fix(isource);
j = fix(jsource);
% add the source term
sigma2(i+1,j+1) = sigma2(i+1,j+1) + source_term.*rbm(i+1,j+1);
for j = 1:ny-1;
for i = 2:nx;
% interpolate material parameters at the right location in the staggered grid cell
c33_half_y = 2.0d0./(1.0d0./rmu(i+1,j+1)+1.0d0./rmu(i+1,j+1+1));
c33_half_y = rmu(i+1,j+1+1);
value_dx_sigmaxy =(27.0d0.*vy(i+1,j+1) - 27.0d0.*vy(i-1+1,j+1)-vy(i+1+1,j+1)+vy(i-2+1,j+1)) ./ deltax./24.0d0;
value_dy_sigmaxy =(27.0d0.*vx(i+1,j+1+1) - 27.0d0.*vx(i+1,j+1)-vx(i+1,j+2+1)+vx(i+1,j-1+1)) ./ deltay./24.0d0;
memory_dx_sigmaxy(i+1,j+1) = b_x(i) .* memory_dx_sigmaxy(i+1,j+1) + a_x(i) .* value_dx_sigmaxy;
memory_dy_sigmaxy(i+1,j+1) = b_y_half_y(j) .* memory_dy_sigmaxy(i+1,j+1) + a_y_half_y(j) .* value_dy_sigmaxy;
sigmaxy(i+1,j+1) = sigmaxy(i+1,j+1) +c33_half_y./1.0d0 .*(value_dx_sigmaxy ./ k_x(i) + memory_dx_sigmaxy(i+1,j+1) +value_dy_sigmaxy ./ k_y(j) + memory_dy_sigmaxy(i+1,j+1)) .* deltat;
end; i =fix(nx+1);
end; j =fix(ny-1+1);
for j = 2:ny;
for i = 1:nx-1;
sigmaxx(i+1,j+1)=(rlambdao(i+1,j+1)+2.0d0.*rmu(i+1,j+1)).*gamma11(i+1,j+1)+rlambdao(i+1,j+1).*gamma22(i+1,j+1) -alpha(i+1,j+1).*sigma2(i+1,j+1);
sigmayy(i+1,j+1)=rlambdao(i+1,j+1).*gamma11(i+1,j+1)+(rlambdao(i+1,j+1)+2.0d0.*rmu(i+1,j+1)).*gamma22(i+1,j+1) -alpha(i+1,j+1).*sigma2(i+1,j+1);
end; i =fix(nx-1+1);
end; j =fix(ny+1);
%------------------
% compute velocity
%------------------
%-----------------------------------
% update memory variables for C-PML
%-----------------------------------
for j = 2:ny;
for i = 2:nx;
co=(rho(i+1,j+1).*rsm(i+1,j+1)-rhof(i+1,j+1).*rhof(i+1,j+1))./deltat;
c1=co+rho(i+1,j+1).*etaokappa(i+1,j+1).*0.5d0;
c2=co-rho(i+1,j+1).*etaokappa(i+1,j+1).*0.5d0;
vtemp=vxf(i+1,j+1);
value_dx_vx1 =(27.0d0.*sigmaxx(i+1,j+1) - 27.0d0.*sigmaxx(i-1+1,j+1)-sigmaxx(i+1+1,j+1)+sigmaxx(i-2+1,j+1)) ./ deltax./24.0d0;
value_dx_vx2 =(27.0d0.*sigma2(i+1,j+1) - 27.0d0.*sigma2(i-1+1,j+1)-sigma2(i+1+1,j+1)+sigma2(i-2+1,j+1)) ./ deltax./24.0d0;
value_dy_vx =(27.0d0.*sigmaxy(i+1,j+1) - 27.0d0.*sigmaxy(i+1,j-1+1)-sigmaxy(i+1,j+1+1)+sigmaxy(i+1,j-2+1)) ./ deltay./24.0d0;
memory_dx_vx1(i+1,j+1) = b_x(i) .* memory_dx_vx1(i+1,j+1) + a_x(i) .* value_dx_vx1;
memory_dx_vx2(i+1,j+1) = b_x(i) .* memory_dx_vx2(i+1,j+1) + a_x(i) .* value_dx_vx2;
memory_dy_vx(i+1,j+1) = b_y(j) .* memory_dy_vx(i+1,j+1) + a_y(j) .* value_dy_vx;
vxf(i+1,j+1) =(c2.*vxf(i+1,j+1) +(-rhof(i+1,j+1).*(value_dx_vx1./ k_x(i) + memory_dx_vx1(i+1,j+1)+ value_dy_vx ./ k_y(j) + memory_dy_vx(i+1,j+1))-rho(i+1,j+1).*(value_dx_vx2./ k_x(i) + memory_dx_vx2(i+1,j+1)))) ./c1;
vtemp=(vtemp+vxf(i+1,j+1)).*0.5d0;
vx(i+1,j+1) = vx(i+1,j+1) +(rsm(i+1,j+1).*(value_dx_vx1./ k_x(i) + memory_dx_vx1(i+1,j+1)+value_dy_vx ./ k_y(j) + memory_dy_vx(i+1,j+1))+rhof(i+1,j+1).*(value_dx_vx2./ k_x(i) + memory_dx_vx2(i+1,j+1)) +rhof(i+1,j+1).*etaokappa(i+1,j+1).*vtemp)./co;
end; i =fix(nx+1);
end; j =fix(ny+1);
for j = 1:ny-1;
for i = 1:nx-1;
rho_half_x_half_y = rho(i+1,j+1+1);
rsm_half_x_half_y = rsm(i+1,j+1+1);
rhof_half_x_half_y = rhof(i+1,j+1+1);
etaokappa_half_x_half_y = etaokappa(i+1,j+1+1);
co=(rho_half_x_half_y.*rsm_half_x_half_y-rhof_half_x_half_y.^2)./deltat;
c1=co+rho_half_x_half_y.*etaokappa_half_x_half_y.*0.5d0;
c2=co-rho_half_x_half_y.*etaokappa_half_x_half_y.*0.5d0;
vtemp=vyf(i+1,j+1);
value_dx_vy =(27.0d0.*sigmaxy(i+1+1,j+1) - 27.0d0.*sigmaxy(i+1,j+1)-sigmaxy(i+2+1,j+1)+sigmaxy(i-1+1,j+1)) ./ deltax./24.0d0;
value_dy_vy1 =(27.0d0.*sigmayy(i+1,j+1+1)- 27.0d0.*sigmayy(i+1,j+1)-sigmayy(i+1,j+2+1)+sigmayy(i+1,j-1+1)) ./ deltay./24.0d0;
value_dy_vy2 =(27.0d0.*sigma2(i+1,j+1+1) - 27.0d0.*sigma2(i+1,j+1)-sigma2(i+1,j+2+1)+sigma2(i+1,j-1+1)) ./ deltay./24.0d0;
memory_dx_vy(i+1,j+1)  = b_x_half_x(i) .* memory_dx_vy(i+1,j+1) + a_x_half_x(i) .* value_dx_vy;
memory_dy_vy1(i+1,j+1) = b_y_half_y(j) .* memory_dy_vy1(i+1,j+1) + a_y_half_y(j) .* value_dy_vy1;
memory_dy_vy2(i+1,j+1) = b_y_half_y(j) .* memory_dy_vy2(i+1,j+1) + a_y_half_y(j) .* value_dy_vy2;
vyf(i+1,j+1) =(c2.*vyf(i+1,j+1) +(-rhof_half_x_half_y.*(value_dx_vy ./ k_x_half_x(i) + memory_dx_vy(i+1,j+1)+value_dy_vy1 ./ k_y_half_y(j) + memory_dy_vy1(i+1,j+1))-rho_half_x_half_y.*(value_dy_vy2 ./ k_y_half_y(j) + memory_dy_vy2(i+1,j+1)))) ./c1;
vtemp=(vtemp+vyf(i+1,j+1)).*0.5d0;
vy(i+1,j+1) = vy(i+1,j+1) +(rsm_half_x_half_y.*(value_dx_vy ./ k_x_half_x(i) + memory_dx_vy(i+1,j+1)+ value_dy_vy1 ./ k_y_half_y(j) + memory_dy_vy1(i+1,j+1))+ rhof_half_x_half_y.*(value_dy_vy2 ./ k_y_half_y(j) + memory_dy_vy2(i+1,j+1))+ rhof_half_x_half_y.*etaokappa_half_x_half_y.*vtemp)./co;
end; i =fix(nx-1+1);
end; j =fix(ny-1+1);
% Dirichlet conditions (rigid boundaries) on the edges or at the outside of the PML layers
vx(1+1,:) = zero;
vx(nx+1,:) = zero;
vx(:,1+1) = zero;
vx(:,ny+1) = zero;
vy(1+1,:) = zero;
vy(nx+1,:) = zero;
vy(:,1+1) = zero;
vy(:,ny+1) = zero;
vxf(1+1,:) = zero;
vxf(nx+1,:) = zero;
vxf(:,1+1) = zero;
vxf(:,ny+1) = zero;
vyf(1+1,:) = zero;
vyf(nx+1,:) = zero;
vyf(:,1+1) = zero;
vyf(:,ny+1) = zero;
% store seismograms
for irec = 1:nrec;
% solid x component
sisvx(it,irec) = vx(ix_rec(irec)+1,iy_rec(irec)+1);
% solid y component
sisvy(it,irec) = vy(ix_rec(irec)+1,iy_rec(irec)+1);
% fluid pressure
sisp(it,irec) = sigma2(ix_rec(irec)+1,iy_rec(irec)+1);
% fluid velocity x
sisvxf(it,irec) = vxf(ix_rec(irec)+1,iy_rec(irec)+1);
% fluid velocity y
sisvyf(it,irec) = vyf(ix_rec(irec)+1,iy_rec(irec)+1);
% solid stress xx
sissigmaxx(it,irec) = sigmaxx(ix_rec(irec)+1,iy_rec(irec)+1);
% solid stress xy
sissigmaxy(it,irec) = sigmaxy(ix_rec(irec)+1,iy_rec(irec)+1);
% solid stress yy
sissigmayy(it,irec) = sigmayy(ix_rec(irec)+1,iy_rec(irec)+1);

end; irec =fix(nrec+1);

% compute total energy
% compute kinetic energy first, defined as 1/2 rho ||v||^2
% in principle we should use rho_half_x_half_y instead of rho for vy
% in order to interpolate density at the right location in the staggered grid cell
% but in a homogeneous medium we can safely ignore it
total_energy_kinetic(it) = 0.5d0 .*sum(sum(rho([npoints_pml:nx-npoints_pml+1]+1,[npoints_pml:ny-npoints_pml+1]+1).*(vx([npoints_pml:nx-npoints_pml+1]+1,[npoints_pml:ny-npoints_pml+1]+1).^2+vy([npoints_pml:nx-npoints_pml+1]+1,[npoints_pml:ny-npoints_pml+1]+1).^2))).*deltax .* deltay+0.5d0.*sum(sum(rsm([npoints_pml:nx-npoints_pml+1]+1,[npoints_pml:ny-npoints_pml+1]+1).*(vxf([npoints_pml:nx-npoints_pml+1]+1,[npoints_pml:ny-npoints_pml+1]+1).^2+vyf([npoints_pml:nx-npoints_pml+1]+1,[npoints_pml:ny-npoints_pml+1]+1).^2))).*deltax.*deltay;
% add potential energy, defined as 1/2 epsilon_ij sigma_ij
% in principle we should interpolate the medium parameters at the right location
% in the staggered grid cell but in a homogeneous medium we can safely ignore it
total_energy_potential(it) = zero;
for j = npoints_pml:ny-npoints_pml+1;
for i = npoints_pml:nx-npoints_pml+1;
epsilon_xx =((rlambdao(i+1,j+1) + 2.0d0.*rmu(i+1,j+1)) .* sigmaxx(i+1,j+1) - rlambdao(i+1,j+1) .* sigmayy(i+1,j+1)) ./(4.0d0 .* rmu(i+1,j+1) .*(rlambdao(i+1,j+1) + rmu(i+1,j+1)));
epsilon_yy =((rlambdao(i+1,j+1) + 2.0d0.*rmu(i+1,j+1)) .* sigmayy(i+1,j+1) - rlambdao(i+1,j+1) .* sigmaxx(i+1,j+1)) ./(4.0d0 .* rmu(i+1,j+1) .*(rlambdao(i+1,j+1) + rmu(i+1,j+1)));
epsilon_xy = sigmaxy(i+1,j+1) ./(2.0d0 .* rmu(i+1,j+1));
total_energy_potential(it) = total_energy_potential(it) +0.5d0 .*(epsilon_xx .* sigmaxx(i+1,j+1) + epsilon_yy .* sigmayy(i+1,j+1) + 2.0d0 .* epsilon_xy .* sigmaxy(i+1,j+1)+sigma2(i+1,j+1).^2./rbm(i+1,j+1)+2.0d0.*rhof(i+1,j+1).*(vx(i+1,j+1).*vxf(i+1,j+1)+vy(i+1,j+1).*vyf(i+1,j+1))).*deltax .* deltay;
end; i =fix(nx-npoints_pml+1+1);
end; j =fix(ny-npoints_pml+1+1);
% output information


if(rem(it,it_display) == 0 || it == 1)
% print maximum of norm of velocity
tempmaxval=sqrt(vx(:,:).^2 + vy(:,:).^2);
velocnorm_all = max(tempmaxval(:));
fprintf(['%s %0.15g %0.15g \n'], 'time step, time = ',it,(it-1).*deltat);
fprintf(['%s %0.15g \n'], 'maximum of norm of velocity is ',velocnorm_all);
fprintf(['%s %0.15g \n'], 'total energy = ',total_energy_kinetic(it) + total_energy_potential(it));
fprintf(['%0.15g \n']);
% check stability of the code, exit if unstable
if(velocnorm_all > stability_threshold)
error(['sinside encountered in original fortran code  ',char(10),' code became unstable and blew up;']);
end;
vnorm(:,:)=sqrt(vx(:,:).^2+vy(:,:).^2);


%imagesc(vx);
%frame=getframe();
%gray = mat2gray(vx);
%X = gray2ind(gray, 256);
%rgb = ind2rgb(X, hot(256));
%filename = sprintf('vx%d.jpg', it)
%imwrite(rgb, filename)

end;

% end of time loop
end; it =fix(nstep+1);



%%
% savemlv seismograms

%save the vx, vy and pf and plot them


% Solid velocity X
    
vx_file=zeros(nstep,2);
fid_20=fopen('./vx1.dat','w+');
for it=1:nstep;
vx_file(it,1)= ((it-1).*deltat - t0);
vx_file(it,2)= sisvx(it,1);
%fprintf(fid_20,['%0.15g %0.15g \n'],vx_file(it,1) , vx_file(it,2));
end; it=fix(nstep+1);
fclose(fid_20);

subplot(4,1,1);
plot(vx_file(:,1),vx_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Solid Velocity X','(in mm/s)'});
title('Solid Velocity X vs. Time AT SOURCE');
%saveas(gcf,'./Vx1.jpg');
%saveas(gca,'./Vx1.jpg');

 %%%%%%%%%%%%%%%%%%
% SOlid velocity X
    
vx_file=zeros(nstep,2);
fid_20=fopen('./vx2.dat','w+');
for it=1:nstep;
vx_file(it,1)= ((it-1).*deltat - t0);
vx_file(it,2)= sisvx(it,2);
%fprintf(fid_20,['%0.15g %0.15g \n'],vx_file(it,1) , vx_file(it,2));
end; it=fix(nstep+1);
fclose(fid_20);


subplot(4,1,2);
plot(vx_file(:,1),vx_file(:,2));
xlabel('time( in seconds)');
ylabel({'Solid Velocity X','( in mm/s)'});
title('Solid Velocity X vs. Time AT EDGE');
%saveas(gcf,'./Vx2.jpg');
%saveas(gca,'./Vx2.jpg');

%%%%%%%%%%%%%%%%%%%%55
% SOlid velocity Y
vy_file=zeros(nstep,2);
fid_20=fopen('./vy1.dat','w+');
for it=1:nstep;
vy_file(it,1)= ((it-1).*deltat - t0);
vy_file(it,2)= sisvy(it,1);
%fprintf(fid_20,['%0.15g %0.15g \n'],vy_file(it,1) , vy_file(it,2));
end; it=fix(nstep+1);
fclose(fid_20);

%plot it
subplot(4,1,3);
plot(vy_file(:,1),vy_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Solid Velocity Y','(in mm/s)'});
title(' Solid Velocity Y vs. Time AT SOURCE');
%saveas(gcf,'./Vy1.jpg');
%saveas(gca,'./Vy1.jpg');

%%%%%%%%%%%%%%%%%%%%%55
% SOlid velocity Y
vy_file=zeros(nstep,2);
fid_20=fopen('./vy2.dat','w+');
for it=1:nstep;
vy_file(it,1)= ((it-1).*deltat - t0);
vy_file(it,2)= sisvy(it,2);
%fprintf(fid_20,['%0.15g %0.15g \n'],vy_file(it,1) , vy_file(it,2));
end; it =fix(nstep+1);
fclose(fid_20);

%plot it
subplot(4,1,4);
plot(vy_file(:,1),vy_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Solid Velocity Y','(in mm/s)'});
title(' Solid Velocity Y vs. Time AT EDGE');
%saveas(gcf,'./Vy2.jpg');
%saveas(gca,'./Vy2.jpg');

%%%%%%%%%%%%%%%%%%%%%5555

filename1= sprintf('2D_Poroelastic_VelocitySolid_type=%s.jpg',sourcetype)
filename2= sprintf('2D_Poroelastic_VelocitySolid_type=%s.fig',sourcetype)
saveas(gcf,filename1);
saveas(gca,filename1);
saveas(gcf,filename2);
saveas(gca,filename2);



%%%%%%%%%%%%%%%%%%%%55



% fluid pressure
pf_file=zeros(nstep,2);
fid_20=fopen('./pf1.dat','w+');
for it=1:nstep;
pf_file(it,1)= ((it-1).*deltat - t0);
pf_file(it,2)= sisp(it,1);
%fprintf(fid_20,['%0.15g %0.15g \n'],pf_file(it,1) , pf_file(it,2));
end; it =fix(nstep+1);
fclose(fid_20);

%plot it
subplot(2,1,1);
plot(pf_file(:,1),pf_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Fluid Pressure','(in Pa)'});
title(' Fluid Pressure vs. Time AT SOURCE' );
%saveas(gcf,'./pf1.jpg');
%saveas(gca,'./pf1.jpg');

%%%%%%%%%%%%%%%%%%%%


% fluid pressure

pf_file=zeros(nstep,2);
fid_20=fopen('./pf2.dat','w+');
for it=1:nstep;
pf_file(it,1)= ((it-1).*deltat - t0);
pf_file(it,2)= sisp(it,2);
%fprintf(fid_20,['%0.15g %0.15g \n'],pf_file(it,1) , pf_file(it,2));
end; it =fix(nstep+1);
fclose(fid_20);

%plot it
subplot(2,1,2);
plot(pf_file(:,1),pf_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Fluid Pressure','(in Pa)'});
title(' Fluid Pressure vs. Time AT EDGE');
%saveas(gcf,'./pf2.jpg');
%saveas(gca,'./pf2.jpg');

%%%%%%%%%%%%%%%%%%55

filename1= sprintf('2D_Poroelastic_Pressure_type=%s.jpg',sourcetype)
filename2= sprintf('2D_Poroelastic_Pressure_type=%s.fig',sourcetype)
saveas(gcf,filename1);
saveas(gca,filename1);
saveas(gcf,filename2);
saveas(gca,filename2);


%%%%%%%%%%%%%%%%%%%%55



% FLuid velocity X
    
vxf_file=zeros(nstep,2);
fid_20=fopen('./vxf1.dat','w+');
for it=1:nstep;
vxf_file(it,1)= ((it-1).*deltat - t0);
vxf_file(it,2)= sisvxf(it,1);
%fprintf(fid_20,['%0.15g %0.15g \n'],vx_file(it,1) , vx_file(it,2));
end; it=fix(nstep+1);
fclose(fid_20);

subplot(4,1,1);
plot(vxf_file(:,1),vxf_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Fluid Velocity X','(in mm/s)'});
title('Fluid Velocity X vs. Time AT SOURCE');
%saveas(gcf,'./Vx1.jpg');
%saveas(gca,'./Vx1.jpg');

 %%%%%%%%%%%%%%%%%%
% Fluid velocity X
    
vxf_file=zeros(nstep,2);
fid_20=fopen('./vxf2.dat','w+');
for it=1:nstep;
vxf_file(it,1)= ((it-1).*deltat - t0);
vxf_file(it,2)= sisvxf(it,2);
%fprintf(fid_20,['%0.15g %0.15g \n'],vx_file(it,1) , vx_file(it,2));
end; it=fix(nstep+1);
fclose(fid_20);


subplot(4,1,2);
plot(vxf_file(:,1),vxf_file(:,2));
xlabel('Time( in seconds)');
ylabel({'Fluid Velocity X','( in mm/s)'});
title('Fluid Velocity X vs. Time AT EDGE');
%saveas(gcf,'./Vx2.jpg');
%saveas(gca,'./Vx2.jpg');

%%%%%%%%%%%%%%%%%%%%55
% Fluid velocity Y
vyf_file=zeros(nstep,2);
fid_20=fopen('./vyf1.dat','w+');
for it=1:nstep;
vyf_file(it,1)= ((it-1).*deltat - t0);
vyf_file(it,2)= sisvyf(it,1);
%fprintf(fid_20,['%0.15g %0.15g \n'],vy_file(it,1) , vy_file(it,2));
end; it=fix(nstep+1);
fclose(fid_20);

%plot it
subplot(4,1,3);
plot(vyf_file(:,1),vyf_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Fluid Velocity Y','(in mm/s)'});
title(' FLuid Velocity Y vs. Time AT SOURCE');
%saveas(gcf,'./Vy1.jpg');
%saveas(gca,'./Vy1.jpg');

%%%%%%%%%%%%%%%%%%%%%55
% Fluid velocity Y
vyf_file=zeros(nstep,2);
fid_20=fopen('./vyf2.dat','w+');
for it=1:nstep;
vyf_file(it,1)= ((it-1).*deltat - t0);
vyf_file(it,2)= sisvyf(it,2);
%fprintf(fid_20,['%0.15g %0.15g \n'],vy_file(it,1) , vy_file(it,2));
end; it =fix(nstep+1);
fclose(fid_20);

%plot it
subplot(4,1,4);
plot(vyf_file(:,1),vyf_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Fluid Velocity Y','(in mm/s)'});
title(' Fluid Velocity Y vs. Time AT EDGE');
%saveas(gcf,'./Vy2.jpg');
%saveas(gca,'./Vy2.jpg');

%%%%%%%%%%%%%%%%%%%%%5555

filename1= sprintf('2D_Poroelastic_VelocityFluid_type=%s.jpg',sourcetype)
filename2= sprintf('2D_Poroelastic_VelocityFluid_type=%s.fig',sourcetype)
saveas(gcf,filename1);
saveas(gca,filename1);
saveas(gcf,filename2);
saveas(gca,filename2);

%%%%%%%%%%%%%%%%%%%%55


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
subplot(3,2,4);
plot(stressxy_file(:,1),stressxy_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Solid Stress XY','(in Pa)'});
title(' Solid Stress XY vs. Time AT EDGE');
%saveas(gcf,'./pf2.jpg');
%saveas(gca,'./pf2.jpg');

%%%%%%%%%%%%%%%%%%55

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
subplot(3,2,6);
plot(stressyy_file(:,1),stressyy_file(:,2));
xlabel('Time (in seconds)');
ylabel({'Solid Stress YY','(in Pa)'});
title(' Solid Stress YY vs. Time AT EDGE');
%saveas(gcf,'./pf2.jpg');
%saveas(gca,'./pf2.jpg');

%%%%%%%%%%%%%%%%%%55

filename1= sprintf('2D_Poroelastic_Stress_type=%s.jpg',sourcetype)
filename2= sprintf('2D_Poroelastic_Stress_type=%s.fig',sourcetype)
saveas(gcf,filename1);
saveas(gca,filename1);
saveas(gcf,filename2);
saveas(gca,filename2);


%%%%%%%%%%%%%%%%%%%%55


fprintf(['%0.15g \n']);
fprintf(['%s \n'], 'End of the simulation');
fprintf(['%0.15g \n']);

clear global a a1_b a1_t a_bottom a_top a_x a_x_half_x a_y a_y_half_y abscissa_in_pml abscissa_normalized alp_sigma2 alpha alpha_bottom alpha_max_pml alpha_top alpha_x alpha_x_half_x alpha_y alpha_y_half_y angle_force b1_b b1_t b_x b_x_half_x b_y b_y_half_y c1 c1_b c1_t c2 c33_half_y co courant_number_bottom courant_number_top cp_bottom cp_top cps_bottom cps_top cs_bottom cs_top d0_x d0_y d_x d_x_half_x d_y d_y_half_y degrees_to_radians delta_b delta_t deltat deltax deltay dispersion_number_bottom dispersion_number_top dist distval epsilon_xx epsilon_xy epsilon_yy etaokappa etaokappa_bottom etaokappa_half_x_half_y etaokappa_top f0 factor ga_b ga_t gamma11 gamma12_1 gamma22 heterogeneous_model hugeval i ideb interface_height irec isource it it_display ix_rec iy_rec j jdeb jinterface jsource k_max_pml k_x k_x_half_x k_y k_y_half_y lambda_b lambda_t max_amplitude memory_dx_sigma2vx memory_dx_sigma2vxf memory_dx_sigmaxx memory_dx_sigmaxy memory_dx_sigmayy memory_dx_vx1 memory_dx_vx2 memory_dx_vy memory_dy_sigma2vy memory_dy_sigma2vyf memory_dy_sigmaxx memory_dy_sigmaxy memory_dy_sigmayy memory_dy_vx memory_dy_vy1 memory_dy_vy2 npoints_pml npower nrec nstep nx ny phi_bottom phi_top pi r_b r_t rbm rbm_bottom rbm_top rcoef rho rho_bottom rho_half_x_half_y rho_top rhof rhof_bottom rhof_half_x_half_y rhof_top rhos_bottom rhos_top rlambdac rlambdac_bottom rlambdac_top rlambdao rlambdao_bottom rlambdao_top rmu rmu_bottom rmu_top ro11_b ro11_t ro12_b ro12_t ro22_b ro22_t rsm rsm_bottom rsm_half_x_half_y rsm_top s_b s_t sigma2 sigmaxx sigmaxy sigmayy sisp sisvx sisvy source_term stability_threshold system_command1 t t0 temp thickness_pml_x thickness_pml_y total_energy_kinetic total_energy_potential use_pml_bottom use_pml_left use_pml_right use_pml_top value_dx_sigma2vxf value_dx_sigmaxx value_dx_sigmaxy value_dx_vx1 value_dx_vx2 value_dx_vy value_dy_sigma2vyf value_dy_sigmaxx value_dy_sigmaxy value_dy_vx value_dy_vy1 value_dy_vy2 velocnorm_all vnorm vtemp vx vxf vy vyf xdeb xfin xi_1 xi_2 xoriginleft xoriginright xrec xsource xspacerec xval ydeb yfin yoriginbottom yorigintop yrec ysource yspacerec yval zero ; 

clear tempmaxval fid_20 ans 

%%


