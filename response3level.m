% Concepts and Methods of 2D Infrared Spectroscopy
% Peter Hamm and Martin T. Zanni
% This matlab program duplicates Figs. 7.8

% define system parameters
Delta_omega = 5; %linewidth (sigma) in *radians* / ps ! 
tau = 1; % correlation time in ps
w_0 = 0; %the center frequency (0 in the ``rotating frame'')
Delta = 5; % anharmonicity

% define desired population time(s) t_2 with stepsize dt and nt
% data points
t2 = 0.3;
dt = 0.1;
n_t = 64;
n_zp = n_t*2; %zero-padded length (optimal is 2*nt, see section ???)
t=0:dt:(n_t-1)*dt;

Lambda = 1/tau;
disp(['Delta_omega / Lambda = ',num2str(Delta_omega/Lambda)]);
disp(['Delta_omega = ',num2str(Delta_omega),' ps^-1']);
disp(['Lambda = ',num2str(Lambda),' ps^-1']);

% define the kubo lineshape function
g = @(t) Delta_omega^2/Lambda^2.*(exp(-Lambda.*t)-1+Lambda.*t);


[T1,T3] = meshgrid(t,t);
R_r=exp(-1i*w_0.*(-T1+T3)).*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));
R_nr=exp(-1i*w_0.*(T1+T3)).*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*Delta.*T3));

%the first time points need to be divided by 2 (Sect 9.5.3)
R_r(:,1) = R_r(:,1)/2;
R_r(1,:) = R_r(1,:)/2;
R_nr(:,1) = R_nr(:,1)/2;
R_nr(1,:) = R_nr(1,:)/2;

%this is what we have calculated so far (response functions in the
%time domain)
figure(1),clf
subplot(1,2,1)
contourf(real(R_r),12);
axis equal tight
subplot(1,2,2)
contourf(real(R_nr),12)
axis equal tight

% do the fft
R_r = fft2(R_r,n_zp,n_zp);
R_nr = fft2(R_nr,n_zp,n_zp);

%now frequency domain
figure(2),clf
subplot(1,2,1)
contourf(fftshift(real(R_r')),12); %pump-probe axis convention
%contourf(fftshift(real(R_r)),12); % the (omega_1, omega_3) axis convention
axis equal tight
subplot(1,2,2)
contourf(fftshift(real(R_nr')),12)
%contourf(fftshift(real(R_nr)),12)
axis equal tight

% flip R1 (being careful to keep zero frequency as the first time
% point), add the response functions, take the real part, and
% finally reorganize so that the 0 frequency is in the center
R = fftshift(real(fliplr(circshift(R_r,[0 -1]))+R_nr));
    

figure(3),clf
contourf(R',12) %use R' to display the pump-probe axis convention
%contourf(R,12) %use R to display the (omega_1, omega_3) axis convention
axis equal tight
 
