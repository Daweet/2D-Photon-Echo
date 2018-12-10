% Concepts and Methods of 2D Infrared Spectroscopy
% Peter Hamm and Martin T. Zanni
% This matlab program duplicates Figs. 4.14 and 10.8

% the number of excited states, their frequencies, and dipoles, all
% taken from Gaussian output. An offset frequency is subtracted
% (and added back after the Fourier transform) to save computation
% time. The two transition dipoles are in the x- and y-directions,
% respectively. 
% define system parameters
Delta_omega = 1; %linewidth (sigma) in *radians* / ps ! 
tau = .1; % correlation time in ps
w_0 = 0; %the center frequency (0 in the ``rotating frame'')
Delta = 5; % anharmonicity
h=1.054571726E-34*5.03445E22;  %hbar in cm-1 s

Lambda = 1/tau;

% define the kubo lineshape function
g = @(t) Delta_omega^2/Lambda^2.*(exp(-Lambda.*t)-1+Lambda.*t);


[T1,T3] = meshgrid(t,t);
R1_r=exp(-1i*w_0.*(-T1+T3)).*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
R1_nr=exp(-1i*w_0.*(T1+T3)).*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
%hamiltonian
J=300;
e1=16200;
e2=15800;
d1=1;
d2=-0.23;
th=0.5*atan(2*J/(e1-e2));
k=cos(th)*sin(th);
E1=e1*cos(th)^2+e2*sin(th)^2+2*J*k; E2=e2*cos(th)^2+e1*sin(th)^2-2*J*k; E3=e1+e2;
wge1=-E1/h; wge2=-E2/h; we1e2=(E1-E2)/h; w3e1=(E3-E1)/h; w3e2=(E3-E2)/h;
%wge1=E1; wge2=E2; we1e2=(E1-E2); wfe1=(Ef-E1); wfe2=(Ef-E2);
Ut=[cos(th) sin(th); -sin(th) cos(th)];
dip=Ut*[d1;d2];
%mu1=dip(1);mu2=dip(2);mu13=sin(th)*d1+cos(th)*d2;mu23=cos(th)*d1-sin(th)*d2;
n = 2;
w = [E1 E2];
mu = [sqrt(0.908) 0 0 ; 0 sqrt(0.671) 0];
w_off = 1915;
w = w - w_off;

% The number of doubly excited states, their frequencies, and their
% transition dipoles taken from the Gaussian output. We use the
% harmonic approximation for the transition dipoles.

n2 = 3;
w2 = [E1 E2 E3];
w2 = w2 - 2*w_off;
mu2 = zeros(n2,3,n);
mu2(:,:,1) = [sqrt(2).*mu(1,:) ; 0,0,0 ; mu(2,:)];
mu2(:,:,2) = [ 0,0,0 ; sqrt(2).*mu(2,:) ; mu(1,:)];


% Parameters of the simulation
n_t = 128;
n_zp = 2*n_t; % the zeropadded length
dt = 0.25;
T2 = 2; %dephasing time
c = 0.188; %unit conversion factor
t2 = 700; %the population time

% Calculate response functions for the rephasing and non-rephasing
% diagrams (Eqs. 10.11, 10.12, and 10.13). The first time points
% need to be halved (Sec. 9.5.3).

R_r = zeros(n_t,n_t);
R_nr = zeros(n_t,n_t);
t=0:dt:(n_t-1)*dt;
[T1,T3] = meshgrid(t,t);

% first calculate all rephasing diagrams

 
    
    mu1 = sqrt(mu(1,:)*mu(1,:)');
    mu2 = sqrt(mu(2,:)*mu(2,:)');
    dipole = mu1^2*mu2^2;
    
    % rephasing diagram R1
    %R_r = R_r - dipole*exp(+ 1i*w(j)*c.*T1 ...
				 %- 1i*w(i)*c.*T3 ...
				 %+ 1i*(w(j)-w(i))*c*t2 ...
				 %- (T1+T3)/T2);
    R_r = R_r - dipole*exp(+ 1i*w(2)*c.*T1 ...
				 - 1i*w(1)*c.*T3 ...
				 + 1i*(w(2)-w(1))*c*t2 ...
				 - (T1+T3)/T2).*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
    R_r = R_r - dipole*exp(+ 1i*w(1)*c.*T1 ...
				 - 1i*w(2)*c.*T3 ...
				 + 1i*(w(1)-w(2))*c*t2 ...
				 - (T1+T3)/T2).*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
    % rephasing diagram R2
    R_r = R_r - dipole*exp(+ 1i*w(2)*c.*T1 ...
				 - 1i*w(1)*c.*T3 ...
				 - (T1+T3)/T2).*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
    
    R_r = R_r - dipole*exp(+ 1i*w(1)*c.*T1 ...
				 - 1i*w(2)*c.*T3 ...
				 - (T1+T3)/T2).*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
      %mujk = sqrt(mu2(k,:,j)*mu2(k,:,j)');
      %muik = sqrt(mu2(k,:,i)*mu2(k,:,i)');
      mu23 = sqrt(mu2(1,:,1)*mu2(1,:,1)');
      mu13 = sqrt(mu2(1,:,1)*mu2(1,:,1)');
      
      dipole = mu1*mu2*mu13*mu23;
      
      
      %rephasing diagram R3
      R_r = R_r + dipole*exp(+ 1i*w(1)*c.*T1 ...
				   - 1i*(w2(3)-w(1))*c.*T3 ...
				   + 1i*(w(1)-w(2))*c.*t2 ...
				   - (T1+T3)/T2).*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
      R_r = R_r + dipole*exp(+ 1i*w(2)*c.*T1 ...
				   - 1i*(w2(3)-w(2))*c.*T3 ...
				   + 1i*(w(2)-w(1))*c.*t2 ...
				   - (T1+T3)/T2).*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
    
 


% now non-rephasing diagrams

    %mu1 = sqrt(mu(1,:)*mu(1,:)');
    %mu2 = sqrt(mu(2,:)*mu(2,:)');
    dipole = mu1^2*mu2^2;
      
    % non-rephasing diagram R4
    R_nr = R_nr - dipole*exp(- 1i*w(1)*c.*T1 ...
				   - 1i*w(1)*c.*T3 ... %?
				   - 1i*(w(1)-w(2))*c*t2 ...
				   - (T1+T3)/T2).*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
    R_nr = R_nr - dipole*exp(- 1i*w(2)*c.*T1 ...
				   - 1i*w(2)*c.*T3 ... %?
				   - 1i*(w(2)-w(1))*c*t2 ...
				   - (T1+T3)/T2).*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
    % non-rephasing diagram R5
    R_nr = R_nr - dipole*exp(- 1i*w(2)*c.*T1 ...
				   - 1i*w(1)*c.*T3 ...
				   - (T1+T3)/T2).*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
    
    R_nr = R_nr - dipole*exp(- 1i*w(1)*c.*T1 ...
				   - 1i*w(2)*c.*T3 ...
				   - (T1+T3)/T2).*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
     %mu23 = sqrt(mu2(1,:,1)*mu2(1,:,1)');
      %mu13 = sqrt(mu2(1,:,1)*mu2(1,:,1)');
      
      dipole = mu1*mu2*mu13*mu23; 
      %non-rephasing diagram R6
      R_nr = R_nr + dipole*exp(- 1i*w(2)*c.*T1 ...
				     - 1i*(w2(3)-w(1))*c.*T3 ...
				     - 1i*(w(2)-w(1))*c.*t2 ...
				     - (T1+T3)/T2).*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
      R_nr = R_nr + dipole*exp(- 1i*w(1)*c.*T1 ...
				     - 1i*(w2(3)-w(2))*c.*T3 ...
				     - 1i*(w(1)-w(2))*c.*t2 ...
				     - (T1+T3)/T2).*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
    

% divide first points (by row and column) by 2
R_r(:,1) = R_r(:,1)./2;
R_r(1,:) = R_r(1,:)./2;
R_nr(:,1) = R_nr(:,1)./2;
R_nr(1,:) = R_nr(1,:)./2;

%what we have so far in the time domain
figure(1),clf
subplot(1,2,1)
contour(real(R_r'),10); 
axis equal tight
subplot(1,2,2)
contour(real(R_nr'),10);
axis equal tight

% do the fft
R_r = ifft2(R_r,n_zp,n_zp); %given the frequency definitions used
                            %above, use the ifft to get the
                            %frequencies right (Mathematica has the
                            %opposite definition of the fft by default)
R_nr = ifft2(R_nr,n_zp,n_zp);

%now frequency domain
figure(2),clf
subplot(1,2,1)
contour(fftshift(real(R_r')),20); %pump-probe axis convention
%contourf(fftshift(real(R_r)),20; % the (omega_1, omega_3) axis convention
axis equal tight
subplot(1,2,2)
contour(fftshift(real(R_nr')),20)
%contourf(fftshift(real(R_nr)),20)
axis equal tight

% flip R_r (being careful to keep zero frequency as the first time
% point), add the response functions, take the real part, and
% finally reorganize so that the 0 frequency is in the center
R = fftshift(real(fliplr(circshift(R_r,[0 -1]))+R_nr));
    

figure(3),clf
n_contours = 40;
MAX = max(abs(R(:)));
level_list = linspace(-MAX,MAX,n_contours+2);
dl = level_list(2)-level_list(1);
cmin =level_list(1)-dl/2;
cmax =level_list(end);

contour(R',level_list) %use R' to display the pump-probe axis convention
%contourf(R,level_list) %use R to display the (omega_1, omega_3) axis convention
caxis([cmin cmax]);			
axis equal tight
