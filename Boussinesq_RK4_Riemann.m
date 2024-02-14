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

clc
clear
close all

tic
% Set up grid:
N = 2^11; dt = 5e-3; L = 120*pi; x = (L*2/N)*(-N/2:N/2-1)';

% Grid for k:
dk = pi/L;                            %Frequencies should be integer multiples of the
k = fftshift(dk*[-N/2:N/2-1])';        %fundamental frequency, 2pi/T, T the period. Here T = 2L.
kinv = 1./(1i*k);
kinv(1) = 0;

%Initial data options below:

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Riemann initial data, RW mean flow (s3 const) + soliton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = 250; ur = 3*(-1)+1/2; ul = 3*(-5/2)+1/2; B = .3; %for now: need ul=0 if starting on left or else we end up with weird double soliton
x0sol = -300; A = 3/2; C = sqrt(-2*ul + 1 - 2/3*A); %C = -sqrt(-2*ur + 1 - 2/3*A)
% need the -2u_l for solitons starting on left!!!! This comes from the
% constraint on the relationship between a and q from KPII, converted to
% Boussinesq. Do the same for solitons starting on right with ur.
vrkp = 0; vr = 3*vrkp;

urkp = 1/3*(ur - 1/2); ulkp = 1/3*(ul - 1/2);
s30 = 2/3*sqrt(-6*urkp)*urkp - vrkp; %s3=s30=const; 1-wave

vlkp = 2/3*sqrt(-6*ulkp)*ulkp - s30;
vl = 3*vlkp;

ap = (ur-ul)/(2*L);
uin = 1/2*(tanh(-(x-x0)*B)+1)*(ul-ur)+ur + A*sech(sqrt(A/6)*(x-x0sol)).^2;    % u(x,0) = tanh function that goes from ul to ur and rightward soliton
uinrw = 1/2*(tanh(-(x-x0)*B)+1)*(ul-ur)+ur;
uinkprw = 1/3*(uinrw-1/2);
%uin doesn't do anything. It's just for personal reference
win = -1/2*sech(-(x-x0)*B).^2*(ul-ur)*B - sqrt(2/3)*A^(3/2)*sech(sqrt(A/6)*(x-x0sol)).^2.*tanh(sqrt(A/6)*(x-x0sol));  
winrw = -1/2*sech(-(x-x0)*B).^2*(ul-ur)*B;
winkprw = 1/3*winrw;
% w(x,0) = u'(x,0) bc w = u_x. So we integrate the pde in w
% v_x = w_t aka v = u_t = vtwiddle_x
%vin = -1/2*sech(-(x-x0)*B).^2*(vl-vr)*B + sqrt(2/3)*C*A^(3/2)*sech(sqrt(A/6)*(x-x0sol)).^2.*tanh(sqrt(A/6)*(x-x0sol));
%vin = zeros(size(win));

%Using the Riemann invariant - do this, it works better!!
vinkprw = 2/3*( 1/2*(-6*uinkprw).^(-1/2).*(-6*winkprw).*uinkprw + sqrt(-6*uinkprw).*winkprw );
vin = 3*vinkprw + sqrt(2/3)*C*A^(3/2)*sech(sqrt(A/6)*(x-x0sol)).^2.*tanh(sqrt(A/6)*(x-x0sol));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%
% % For the soliton trajectory:
% alkp = 1/3*A;
% qlkp = -sqrt(-6*ulkp - 2*alkp); %=C
% s20 = 2*qlkp*ulkp + 4/9*qlkp^3 - vlkp;
% thetafun = @(x,y) abs(y./(x-x0)).^3.*(((x-x0)./y).^3 - 9*(s30-s20));
% qfun = @(z) abs((x-x0)./y).*cos( 1/3*acos(theta((x-x0)./y)) - 4*pi/3);


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
  nmax = round(150/dt); tmax = nmax*dt; nplt = floor((tmax/tpts)/dt); %nmax = round(tmax/dt);
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
  xlabel x, ylabel t, axis([-L L 0 tmax ul-1/2 ur+1/2]), grid off
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


  %% Write video file for specific data
close all
vid = VideoWriter('RWmoviept3','Motion JPEG AVI');

for i = 1:length(tdata)
    pframe = plot(x,kpudata(:,i),'k','LineWidth',1);
    xlim([x(1) x(end)]); ylim([1/3*(ul-1/2)-1/3 1/3*(ur-1/2)+1/3])
    xlabel('$x$','interpreter','latex'); ylabel('$u(x,t)$','interpreter','latex');
    title("$y=$"+ tdata(i),'interpreter','latex')

    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    currframe = getframe(ax,rect);
    open(vid)
    writeVideo(vid,currframe);
end
close(vid)
%% Comparison to exact solution - 1-RW
close all

zetamin = -sqrt(-6*urkp);
zetamax = -sqrt(-6*ulkp);

tlast = 1;
ylast = tdata(tlast);
xright = ylast*zetamin+x0;
xleft = ylast*zetamax+x0;

xslopeind = xleft <= x & x <= xright;
xslope = x(xslopeind);
zetaslope = (xslope-x0)./ylast;

uex = NaN(size(kpudata));
qex = NaN(size(kpudata));

alkp = 1/3*A;
qlkp = -sqrt(-6*ulkp - 2*alkp); %=C
s20 = 2*qlkp*ulkp + 4/9*qlkp^3 - vlkp;
theta = abs(ylast./(x-x0)).^3.*(((x-x0)./ylast).^3 - 9*(s30-s20));

ind = abs(x-xright) < 2*L/N;
qrkp = max(abs((x(ind)-x0)/ylast).*( cos( 1/3*acos(theta(ind))-4*pi/3) ));

for i = 1:size(uex,1)
    if x(i) < xleft
        uex(i) = ulkp;
        qex(i) = qlkp;
    elseif x(i) > xright
        uex(i) = urkp;
        qex(i) = qrkp;
    else 
        uex(i) = -1/6*((x(i)-x0)/ylast)^2;
        qex(i) = abs((x(i)-x0)/ylast).*cos( 1/3*acos(theta(i))-4*pi/3) ;
    end
end

aex = -1/2*qex.^2 - 3*uex;

figure(5)
hold on
col = [ 0.0078    0.3569    0.5882];
p1 = plot(x,uex,'k','LineWidth',1.5,'LineStyle','-');
p2 = plot(x,kpudata(:,tlast),'b','LineWidth',0.8,'LineStyle','-');
p3 = plot(x,uex + aex,'k','LineWidth',1,'LineStyle','-.');
xlabel('$x$','interpreter','latex'); ylabel('$u(x,y)$','Interpreter','latex')
title("t = "+ tdata(tlast), 'interpreter','latex')
legend([p1(1),p2(1),p3(1)],'exact RW','numerical RW/soliton interaction','soliton amplitude','Location','southeast')

%% Video comparison to exact solution - 1-RW
close all
vid = VideoWriter('RWcompwide','Motion JPEG AVI');

for i = 1:length(tdata)
    zetamin = -sqrt(-6*urkp);
    zetamax = -sqrt(-6*ulkp);

    ylast = tdata(i);
    xright = ylast*zetamin+x0;
    xleft = ylast*zetamax+x0;

    uex = NaN(size(kpudata));
    qex = NaN(size(kpudata));

    alkp = 1/3*A;
    qlkp = -sqrt(-6*ulkp - 2*alkp); %=C
    s20 = 2*qlkp*ulkp + 4/9*qlkp^3 - vlkp;
    theta = abs(ylast./(x-x0)).^3.*(((x-x0)./ylast).^3 - 9*(s30-s20));

    ind = abs(x-xright) < 2*L/N;
    qrkp = max(abs((x(ind)-x0)/ylast).*( cos( 1/3*acos(theta(ind))-4*pi/3) ));

    for j = 1:size(uex,1)
        if x(j) < xleft
            uex(j) = ulkp;
            qex(j) = qlkp;
        elseif x(j) > xright
            uex(j) = urkp;
            qex(j) = qrkp;
        else
            uex(j) = -1/6*((x(j)-x0)/ylast)^2;
            qex(j) = abs((x(i)-x0)/ylast).*cos( 1/3*acos(theta(i))-4*pi/3) ;
        end
    end

    pnum = plot(x,kpudata(:,i),'k','LineWidth',1,'LineStyle','-');
    xlabel('$x$','interpreter','latex'); ylabel('$u(x,y)$','interpreter','latex');
    title("$y=$"+ tdata(i),'interpreter','latex')
    hold on
    pex = plot(x,uex,'k','LineWidth',1,'LineStyle','-.');
    pampl = plot(x,uex + aex,'k','LineWidth',1,'LineStyle','--');
    legend([pnum(1),pex(1)],'numerical','exact','amplitude','Location','southeast')
    %xlim([x(1) x(end)]); 
    %ylim([1/3*(ul-1/2)-1/3 1/3*(ur-1/2)+1/3])
    ylim([-3 0.5])
    hold off

    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    currframe = getframe(ax,rect);
    open(vid)
    writeVideo(vid,currframe);

end
close(vid)


