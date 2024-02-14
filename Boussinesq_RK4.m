%Solve Boussinesq equation with exponentially decaying initial conditions:
% u_tt - u_xx + a(u^2)_xx + g u_xxxx = 0, adapted to the system: 
%
% u_t = v_x
% v_t = -gam*u_xxx + u_x - alpha*(u^2)_xx
%
% on [-L,L] with initial conditions u(x,0) and v(x,0). Can integrate u_t(x,0) in x
% using, e.g., Mathematica, to obtain initial data v(x,0). The script then
% plots the results in x and t.
%
% This is a Fourier spectral method. Once the linear terms
% are eliminated via integrating factor, the resulting system is integrated
% using a slightly adapted Runge-Kutta fourth-order time-stepping scheme.
%
% Based on the KdV solver in Trefethen, Pogram 27 (p27.m)

clc
clear
close all

tic
% Set up grid:
 N = 2^8; dt = 1e-3; L = 20*pi; x = (L*2/N)*(-N/2:N/2-1)';
 clf, drawnow, set(gcf,'renderer','zbuffer')

  %two-soliton initial data
%   A = 0.3; B = 0.2; x0 = -15; y0 = 20; b = 0; c = sqrt((3-2*A)/3); d = -sqrt((3-2*B)/3);
%   u = A*sech(sqrt(A/6)*(x-x0)).^2 + B*sech(sqrt(B/6)*(x-y0)).^2;
%   v = -A*sqrt((6-4*A)/6)*sech(sqrt(A/6)*(x-x0)).^2 + B*sqrt((6-4*B)/6)*sech(sqrt(B/6)*(x-y0)).^2;
%   uh = fft(u); vh = fft(v);

%single soliton initial data
    A = .5; x0 = 0; C = sqrt((3-2*A)/3);
    u = A*sech(sqrt(A/6)*(x-x0)).^2;
    v = -A*C*sech(sqrt(A/6)*(x-x0)).^2;
    
% Exact soliton solution:
uex = @(x,t) A*sech(sqrt(A/6)*(x-x0-C*t)).^2;
        
%Take Fourier transform and set up grid for k  
uh = fft(u); vh = fft(v);
dk = pi/L;                            %Frequencies should be integer multiples of the
k = dk*[0:N/2-1 0 -N/2+1:-1]';        %fundamental frequency, 2pi/T, T the period. Here T = 2L.
  

%Good Boussinesq:
alpha = 1; gam = 1;

%Preallocating the matrix exponential:
m = sqrt(gam*k.^2+ones(size(k)));
km = k.*m;

E1 = .5*(exp(km*.5i*dt) + exp(-km*.5i*dt));
E2 = .5./m.*(exp(km*.5i*dt) - exp(-km*.5i*dt));
E3 = .5*m.*(exp(km*.5i*dt) - exp(-km*.5i*dt));
E4 = .5*(exp(km*.5i*dt) + exp(-km*.5i*dt));

E12 = .5*(exp(km.*1i*dt) + exp(-km*1i*dt));
E22 = .5./m.*(exp(km.*1i*dt) - exp(-km.*1i*dt));
E32 = .5*m.*(exp(km.*1i*dt) - exp(-km.*1i*dt));
E42 = .5*(exp(km.*1i*dt) + exp(-km.*1i*dt));


% Solve PDE and plot results:
  nmax = round(50/dt); tmax = nmax*dt; nplt = floor((tmax/25)/dt); %nmax = round(tmax/dt);
  udata = u; tdata = 0; h = waitbar(0,'please wait...');
  errvec = 0;
  for n = 1:nmax
      t = n*dt; g = -dt*1i*alpha*k;
      a = g.*fft(real(ifft(uh)).^2);
      b = g.*fft(real(ifft(E1.*uh + E2.*vh + E2.*a/2)).^2);     %RK4 with 
      c = g.*fft(real(ifft(E1.*uh + E2.*vh )).^2);              %integrating factor
      d = g.*fft(real(ifft(E12.*uh + E22.*vh + E2.*c)).^2);

      uh0 = uh;
      vh0 = vh;

      uh = E12.*uh0 + E22.*vh0 + 1/6*(E22.*a + 2*E2.*(b+c));
      vh = E32.*uh0 + E42.*vh0 + 1/6*(E42.*a + 2*E4.*(b+c)+d);

    if mod(n,nplt) == 0 
      u = real(ifft(uh)); waitbar(n/nmax)
      udata = [udata u]; tdata = [tdata t];
      errvec = [errvec norm(udata(:,end) - uex(x,t),inf)];
    end
  end
  waterfall(x,tdata,udata'), colormap(1e-6*[1 1 1]); view(-20,25)
  xlabel x, ylabel t, axis([-L L 0 tmax 0 1.5]), grid off
  set(gca,'ztick',[0 1.5]), close(h), pbaspect([1 1 .13]) 
  
  %Plot uh vs Fourier nodes; should be decaying rapidly
  figure(2);
  clf();
  uhat = fft(udata(:,end));
  semilogy(k,abs(uhat)/max(abs(uhat)),'k-');
  xlabel('$k$','interpreter','latex');
  ylabel('$|\hat{u}|$','interpreter','latex');
%%
  figure(3)
  clf()
  scatter(tdata, errvec)
  set(gca,'YScale','log')

tend = toc;
disp(tend)


  

