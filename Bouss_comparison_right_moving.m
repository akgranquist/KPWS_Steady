% Solve Boussinesq equation with Riemann initial conditions:
% u_tt - u_xx + a(u^2)_xx + g u_xxxx = 0, adapted to the system: 
%
% w_t = v_x
% v_t = w_x - gam*w_xxx - 2*alpha*(w^2 + antideriv(w)*w_x
%
% where w = u_x and v_x = w_t (ie v = u_t), and IC
%            
%           / ul  x<0
% u(x,0) = {
%           \ ur  x>0   
%            
% This is a Fourier spectral method. Once the linear terms
% are eliminated via integrating factor, the resulting system is integrated
% using a slightly adapted Runge-Kutta fourth-order time-stepping scheme.
%
% ref: mkdvb_nodir, mkdvb_run_kdv_dsw, Trefethen p27.
%
% This code evaluates a RW + soliton BC starting on the left and going to
% the right.

clc
clear
close all

tic
% Set up grid:
N = 2^12; dt = 5e-3; L = 240*pi; x = (L*2/N)*(-N/2:N/2-1)';
tmax =150;

% Grid for k:
dk = pi/L;                            %Frequencies should be integer multiples of the
k = fftshift(dk*[-N/2:N/2-1])';        %fundamental frequency, 2pi/T, T the period. Here T = 2L.
kinv = 1./(1i*k);
kinv(1) = 0;

%Initial data options below:

% Riemann initial data (RW only)

% %%%%%%%%%%%%%%%%%%%%%%%%%
% % Riemann initial data
% %%%%%%%%%%%%%%%%%%%%%%%
% % Rarefaction wave only - s3 constant
% x0 = 0; ur = 3*(-1)+1/2; ul = 3*(-5/2)+1/2;  B = .1; %A = 0.5;
% vrkp = 0; vr = 3*vrkp;
% 
% urkp = 1/3*(ur - 1/2); ulkp = 1/3*(ul - 1/2); %The KP right and left values
% s30 = 2/3*sqrt(-6*urkp)*urkp - vrkp; %s3=s30=const; 1-wave. To get
% % 3-wave, make both 2/3's negative (this is setting s1=s10=const) and
% % will give only the DSW
% 
% vlkp = 2/3*sqrt(-6*ulkp)*ulkp - s30;
% vl = 3*vlkp;
% %vl = -0.1; %Use this if you don't want to isolate the RW
% 
% ap = (ur-ul)/(2*L);
% uin = 1/2*(tanh(-(x-x0)*B)+1)*(ul-ur)+ur;    % u(x,0) = tanh function that goes from ul to ur
% uinkp = 1/3*(uin-1/2);
% 
% win = -1/2*sech(-(x-x0)*B).^2*(ul-ur)*B;     % w(x,0) = u'(x,0) bc w = u_x. So we integrate the pde in w
% winkp = 1/3*win;
% % v_x = w_t aka v = u_t
% %vinp = 1/2*(tanh(-(x-x0)*B)+1)*(vl-vr)+vr;  % This is the v st v_x = u_t from
% %the evolution form. But v in this code is defined as v = u_t. So we want
% %the x derviative of this tanh function
% %vin = -1/2*sech(-(x-x0)*B).^2*(vl-vr)*B;
% %vin = zeros(size(win));
% 
% %using the Riemann invariant - do this, it works better!
% vinkp = 2/3*( 1/2*(-6*uinkp).^(-1/2).*(-6*winkp).*uinkp + sqrt(-6*uinkp).*winkp ); %this is vbar_x
% vin = 3*vinkp;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%

% Soliton only:

% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % Just a soliton
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0 = 0; ur = -1/2; ul = -1/2; %for now: need ul-ur=0 or we get weird
% % double soliton
% x0sol = -30; A = .5; C = sqrt(-2*ul + 1-2/3*A); 
% % need the -2u_l for solitons starting on left!!!! This comes from the
% % constraint on the relationship between a and q from KPII, converted to
% % Boussinesq. This ensures that we are modulating an exact solution.
% vrkp = 0; vr = 3*vrkp;
% 
% urkp = 1/3*(ur - 1/2); ulkp = 1/3*(ul - 1/2);
% s30 = 2/3*sqrt(-6*urkp)*urkp - vrkp; %s3=s30=const; 1-wave
% 
% vlkp = 2/3*sqrt(-6*ulkp)*ulkp - s30;
% vl = 3*vlkp;
% 
% ap = (ur-ul)/(2*L);
% uin = A*sech(sqrt(A/6)*(x-x0sol)).^2 + ul;    % u(x,0) = sech^2 function (soliton)
% win = -sqrt(2/3)*A^(3/2)*sech(sqrt(A/6)*(x-x0sol)).^2.*tanh(sqrt(A/6)*(x-x0sol));     % w(x,0) = u'(x,0) bc w = u_x. So we integrate the pde in w
% % v_x = w_t aka v = u_t
% vin = sqrt(2/3)*C*A^(3/2)*sech(sqrt(A/6)*(x-x0sol)).^2.*tanh(sqrt(A/6)*(x-x0sol));
% %vin = zeros(size(win));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Riemann initial data + soliton (RW/soliton interaction):

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Riemann initial data, RW mean flow (s3 const) + soliton (starting on
% % left and moving right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = 500; ur = 3*(-1)+1/2; ul = 3*(-5/2)+1/2; B = .3;
x0sol = -250; A = 3/2; C = sqrt(-2*ul + 1 - 2/3*A); %C = -sqrt(-2*ur + 1 - 2/3*A) %for left-moving soliton
% need the -2u_l for solitons starting on left!!!! This comes from the
% constraint on the relationship between a and q from KPII, converted to
% Boussinesq. Do the same for solitons starting on right with ur.
vrkp = 0; vr = 3*vrkp;

urkp = 1/3*(ur - 1/2); ulkp = 1/3*(ul - 1/2);
s30 = 2/3*sqrt(-6*urkp)*urkp - vrkp; %s3=s30=const; 1-wave

vlkp = 2/3*sqrt(-6*ulkp)*ulkp - s30;
vl = 3*vlkp;

ap = (ur-ul)/(2*L);
uin = 1/2*(tanh(-(x-x0)*B)+1)*(ul-ur)+ur + A*sech(sqrt(A/6)*(x-x0sol)).^2;    
% u(x,0) = tanh function that goes from ul to ur and rightward soliton
uinrw = 1/2*(tanh(-(x-x0)*B)+1)*(ul-ur)+ur; %just the rw part of uin
uinkprw = 1/3*(uinrw-1/2); %converted to kp data
win = -1/2*sech(-(x-x0)*B).^2*(ul-ur)*B - sqrt(2/3)*A^(3/2)*sech(sqrt(A/6)*(x-x0sol)).^2.*tanh(sqrt(A/6)*(x-x0sol));  
winrw = -1/2*sech(-(x-x0)*B).^2*(ul-ur)*B;
winkprw = 1/3*winrw;
% w(x,0) = u'(x,0) bc w = u_x. So we integrate the pde in w
% v_x = w_t aka v = u_t = vtwiddle_x
% old (exact) v(x,0):
%vin = -1/2*sech(-(x-x0)*B).^2*(vl-vr)*B + sqrt(2/3)*C*A^(3/2)*sech(sqrt(A/6)*(x-x0sol)).^2.*tanh(sqrt(A/6)*(x-x0sol));
%vin = zeros(size(win));

%Using the Riemann invariant to determine v(x,0) - do this, it works better!!
vinkprw = 2/3*( 1/2*(-6*uinkprw).^(-1/2).*(-6*winkprw).*uinkprw + sqrt(-6*uinkprw).*winkprw );
vin = 3*vinkprw + sqrt(2/3)*C*A^(3/2)*sech(sqrt(A/6)*(x-x0sol)).^2.*tanh(sqrt(A/6)*(x-x0sol));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
% To calculate the analytical soliton trajectory:
alkp = 1/3*A;   %starting amp on left, for kp soliton
qlkp = -sqrt(-6*ulkp - 2*alkp); %=C     %startin q on left - determined by alkp
s20 = 2*qlkp*ulkp + 4/9*qlkp^3 - vlkp;
thetafun = @(x,y) abs(y./(x)).^3.*(((x)./y).^3 - 9*(s30-s20));  %second Riemann invariant (for kp)
qfun = @(x,y) abs((x)./y).*cos( 1/3*acos(thetafun(x,y)) - 4*pi/3); %analytical solution q(x/y) (for kp)
xsol = x0sol-x0;
dx = 2*L/N;
xtraj = xsol;
ttraj = 0:dt:tmax;

zetamin = -sqrt(-6*urkp);   % x/y at right edge
zetamax = -sqrt(-6*ulkp);   % x/y at left edge

qex = NaN(1,N);
qexdata = [];

yint = (x0sol-x0)/(zetamax + qlkp);
xint = zetamax*yint + x0;       %point where soliton intersects RW
% yout = (x0sol-x0)/(zetamin + qrkp);
% xout = zetamin*you + x0; 
xleftvec = x0;
xrightvec = x0;
yleftvec = 0;
yrightvec = 0;

% The soliton trajectory is calculated inside the loop using the above.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = uin; w = win; v = vin;

%Take Fourier transform  
wh = fft(w); vh = fft(v);

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



% Solve PDE:
 %tpts = 25;
 tpts = 300;
 %tpts = 80;
  nmax = round(tmax/dt);  nplt = floor((tmax/tpts)/dt); %nmax = round(tmax/dt); tmax = nmax*dt;
  udata = u; tdata = 0; wait1 = waitbar(0,'please wait...'); 
  errl = abs(u(1)-ul);
  for n = 1:nmax
      t = n*dt; g = -dt*2*alpha;

      ua = ul + ifft([-sum(x.*w); wh(2:end).*kinv(2:end)],...
          'symmetric') + ap*(x+L);
      a = g.*fft( ifft(wh,'symmetric').^2  + ifft(1i*k.*wh,'symmetric').*ua );

            bwh = E1.*wh + E2.*vh + E2.*a/2;
      ub = ul + ifft([-sum(x.*w); bwh(2:end).*kinv(2:end)],...
          'symmetric') + ap*(x+L);  
      b = g.*fft( ifft(bwh,'symmetric').^2  +  ifft(1i*k.*bwh,'symmetric').*ub ); 


        cwh = E1.*wh + E2.*vh;
      uc = ul + ifft([-sum(x.*w);cwh(2:end).*kinv(2:end)],...
          'symmetric') + ap*(x+L);
      c = g.*fft( ifft(cwh,'symmetric').^2 + ifft(1i*k.*cwh,'symmetric').*uc );
       

        dwh = E12.*wh + E22.*vh + E2.*c;
      ud = ul + ifft([-sum(x.*w); dwh(2:end).*kinv(2:end)],...
          'symmetric') + ap*(x+L);
      d = g.*fft( ifft(dwh,'symmetric').^2 + ifft(1i*k.*dwh,'symmetric').*ud );

      wh0 = wh;
      vh0 = vh;

      wh = E12.*wh0 + E22.*vh0 + 1/6*(E22.*a + 2*E2.*(b+c));
      vh = E32.*wh0 + E42.*vh0 + 1/6*(E42.*a + 2*E4.*(b+c)+d);
    
      w = ifft(wh, 'symmetric');

    %soliton trajectory
      xright = ttraj(n)*zetamin+x0;
      xleft = ttraj(n)*zetamax+x0;
      yright = (xright-x0)/zetamin;
      yleft = (xleft-x0)/zetamax;
      yleftvec = [yleftvec yleft]; yrightvec = [yrightvec yright];
      xleftvec = [xleftvec xleft]; xrightvec = [xrightvec xright];

      xsol_old = xsol;
      if xsol_old < xleft-x0
          qex(n) = qlkp;
          xsol = xsol - dt*qlkp;      
      elseif xsol_old > xright-x0
          qrkp = qex(n-1);
          qex(n) = qrkp;
          xsol = xsol-dt*qrkp;
      else
          qex(n) = qfun(xsol,t);
          %xsol = xsol - dt*qfun(xsol,t);       %forward Euler
          c1 = -qfun(xsol,t);
          c2 = -qfun(xsol + c1*dt/2, t + dt/2);
          c3 = -qfun(xsol + c2*dt/2, t + dt/2);
          c4 = -qfun(xsol + c3*dt, t + dt);
          xsol = xsol + dt/6*(c1 + 2*c2 + 2*c3 + c4);       %RK4
          xout = xright;
      end

      xtraj = [xtraj xsol];

    %collect data at select time points
    if mod(n,nplt) == 0 
      u = ul + ifft([-sum(x.*w); wh(2:end).*kinv(2:end)],...
          'symmetric') + ap*(x+L); 
      waitbar(n/nmax)
      udata = [udata u]; tdata = [tdata t]; errl = [errl abs(u(1)-ul)];
     

    end

  end
  close(wait1)
  tend = toc;
  disp(['Computation time: ',num2str(floor(tend/60)), ' min, ',...
      num2str( tend-floor(tend/60)*60 ), ' s' ])

 %Convert udata to KP:
 kpudata = 1/3*(udata-1/2);


%% Plot results:

  % plot solutions to Boussinesq:
  figure(1)
  %waterfall(x,tdata,udata'), colormap(1e-6*[1 1 1]); view(-20,25)
  ss = surf(x,tdata,kpudata');
  set(ss,'EdgeColor','none')
  xlabel x, ylabel t, grid off
  %set(gca,'ztick',[-1 0 1]), pbaspect([1 1 .13]) 
  title('Boussinesq solution')
  
  figure(2)
  plot(tdata,errl','k')
  set(gca,'YScale','log')
  xlabel('$t$','interpreter','latex');
  ylabel('$|u_l - u(1)|$','interpreter','latex');
  title('error at left boundary')

  % Convert to solution of KPII and plot:
  figure(3)
  pp=pcolor(x,tdata,1/3*(udata'-1/2));
  colorbar
  xlabel x, ylabel y, %clim([1/3*ul-1/2 1/3*ur-1/2]), grid off
  set(pp, 'EdgeColor', 'none');
  title('KPII solution')

  figure(4)
  %waterfall(x,tdata,kpudata'); colormap(1e-6*[1 1 1]); view(-20,25)
  ss = surf(x,tdata,kpudata');
  set(ss,'EdgeColor','none')
  xlabel x, ylabel y, axis([-L L 0 tmax 1/3*(ul-1/2)-1/2 1/3*(ur-1/2)+1/2]), grid off
  %set(gca,'ztick',[-1 0 1]), pbaspect([1 1 .13]) 
  title('KPII solution')

%overlay exact soliton trajectory and 
  figure(5) 
  %plot(xtraj+x0,ttraj)
  hold on
  %plot(xint*ones(size(xtraj)), ttraj)
  pp=pcolor(x,tdata,1/3*(udata'-1/2));
  set(pp, 'EdgeColor', 'none');
  colorbar
  plot(xleftvec,yleftvec,'k','LineWidth',1,'LineStyle','--');
  plot(xrightvec,yrightvec,'k','LineWidth',1,'LineStyle','--');
  plot(xtraj+x0,ttraj,'r','LineWidth',1,'LineStyle','--');

