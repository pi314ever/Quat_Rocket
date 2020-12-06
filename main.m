%% Main
% 
%% Setup
clc;clear;close all
L = 2; % Length of rocket [m]
r = 0.2; % Rocket radius [m]
m = 5; % Mass of rocket [kg]
zcm = L/2; % Distance from Thrust to cm [m]
I = [1/4*m*r^2+1/3*m*L^2 1/4*m*r^2+1/3*m*L^2 1/2*m*r^2]; % Mass properties

ref = [pi/4 pi/6 0]; % Reference attitude [rad]

% Gimbal angles
% omega = 2*pi*1; 
% psg = @(t) 0.2*sin(omega*t); 
% thg = @(t) -0.5*cos(omega*1.3*t);
psg = @(t) 0.5;
thg = @(t) 0.5;

% Thrust
T = @(t) 15*(t<4); % Thrust for 4 seconds [N]

tspan = [0 10]; 
x0 = zeros(6,1); % Initial conditions
h = 1E-2; % Time step
%% Simulate
sol = RK4(@(t,x)EoM(t,x,I,zcm,psg,thg,T),tspan,x0,h);
%sol = RK4(@(t,x)EoMfb(t,x,I,zcm,ref,T),tspan,x0,h);
%% Extract quat info
th1 = sol.y(1,:);
th2 = sol.y(2,:);
th3 = sol.y(3,:);
e = [cos(th1/2).*cos(th2/2).*cos(th3/2)-sin(th1/2).*sin(th2/2).*sin(th3/2);
    sin(th1/2).*cos(th2/2).*cos(th3/2)+sin(th2/2).*sin(th3/2).*cos(th1/2);
    sin(th2/2).*cos(th1/2).*cos(th3/2)-sin(th1/2).*sin(th3/2).*cos(th2/2);
    sin(th1/2).*sin(th2/2).*cos(th3/2)+sin(th3/2).*cos(th1/2).*cos(th2/2)];
th = 2*acos(e(1,:));
lam = e(2:4,:)./sin(th/2);

%% Rocket attitude
na1 = @(psi) [1 0 0;0 cos(psi) -sin(psi);0 sin(psi) cos(psi)];
a1a2 = @(theta) [cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)];
a2b = @(phi) [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0;0 0 1];
tr = @(psi,theta,phi)na1(psi)*a1a2(theta)*a2b(phi)*[0;0;1];
bz = zeros(3,length(th1));
for ii = 1:length(th1)
    bz(:,ii) = tr(th1(ii),th2(ii),th3(ii));
end

%% Plot 
figure
plot(sol.x,e);
title('Euler parameters vs Time')
xlabel('Time [s]');ylabel('Magnitude')
legend('e_0','e_1','e_2','e_3','Location','Bestoutside')
%% Animation
refpos = tr(ref(1),ref(2),ref(3));
figure
[x,y,z] = sphere(128);
h = surfl(x, y, z); 
h.FaceAlpha = 0.10;
%h.CData = 0.05*ones(129);
shading interp
hold on
% plot3(refpos(1),refpos(2),refpos(3),'Om','markersize',8) % Reference attitude
plot3([-1 1],[0 0],[0 0],'k--','Handlevisibility','off')
plot3([0 0],[-1 1],[0 0],'k--','Handlevisibility','off')
plot3([0 0],[0 0],[-1 1],'k--','Handlevisibility','off')
RA = plot3([1 -1]*bz(1,8),[1 -1]*bz(2,4),[1 -1]*bz(3,8),'color',[0 1 1],'Linewidth',2);
RA1 = plot3([1 0]*bz(1,5),[1 0]*bz(2,3),[1 0]*bz(3,5),'color',[0.4 1 1],'Linewidth',1.6);
RA2 = plot3([1 0]*bz(1,3),[1 0]*bz(2,2),[1 0]*bz(3,3),'color',[0.6 1 1],'Linewidth',1.4);
RA3 = plot3([1 0]*bz(1,1),[1 0]*bz(2,1),[1 0]*bz(3,1),'color',[0.8 1 1],'Linewidth',1.2);
RH = plot3(bz(1,8),bz(2,8),bz(3,8),'bO','Markersize',5);
LAM = plot3([0 lam(1,8)]*th(1)/(2*pi),[0 lam(2,8)]*th(1)/(2*pi),[0 lam(3,8)]*th(1)/(2*pi),'rx-','Linewidth',1.3);
tstr = title('Attitude animation @ t = 0 s');
for ii = 9:1:length(th1)
    RA.XData = [1 -1]*bz(1,ii);
    RA.YData = [1 -1]*bz(2,ii);
    RA.ZData = [1 -1]*bz(3,ii);
    RA1.XData = [1 0]*bz(1,ii-3);
    RA1.YData = [1 0]*bz(2,ii-3);
    RA1.ZData = [1 0]*bz(3,ii-3);
    RA2.XData = [1 0]*bz(1,ii-5);
    RA2.YData = [1 0]*bz(2,ii-5);
    RA2.ZData = [1 0]*bz(3,ii-5);
    RA3.XData = [1 0]*bz(1,ii-7);
    RA3.YData = [1 0]*bz(2,ii-7);
    RA3.ZData = [1 0]*bz(3,ii-7);
    RH.XData = bz(1,ii);
    RH.YData = bz(2,ii);
    RH.ZData = bz(3,ii);
    LAM.XData = [0 lam(1,ii)]*th(ii)/(2*pi);
    LAM.YData = [0 lam(2,ii)]*th(ii)/(2*pi);
    LAM.ZData = [0 lam(3,ii)]*th(ii)/(2*pi);
    tstr.String = sprintf('Attitude animation @ t = %.2f s',sol.x(ii));
    axis equal
    drawnow
    pause(0.01)
end

%% EoM
function dx = EoM(t,x,I,zcm,psg,thg,T)
%    ps th ph
% f  1  2  3  
% f' 4  5  6  
dx = zeros(6,1);
ps = x(1);
th = x(2);
ph = x(3);
psd = x(4);
thd = x(5);
phd = x(6);
psg = psg(t);
thg = thg(t);
T = T(t);

Ixx = I(1);
Iyy = I(2);
Izz = I(3);
P = thd*sin(ph)+psd*cos(ph)*cos(th);
Q = thd*cos(ph)-psd*cos(th)*sin(ph);
R = phd+psd*sin(th);
for ii = 1:3
    dx(ii) = x(ii+3);
end

% Transformations G -> B
bc = [1 0 0;0 cos(psg) -sin(psg);0 sin(psg) cos(psg)];
cg = [cos(thg) 0 sin(thg);0 1 0;-sin(thg) 0 cos(thg)];

Thrustb = bc*cg*[0;0;T];

A = [Ixx*cos(ph)*cos(th) Ixx*sin(ph) 0;
    -Iyy*cos(th)*sin(ph) Iyy*cos(ph) 0;
    Izz*sin(th) 0 Izz];
b = [zcm*Thrustb(2)+Q*R*(Iyy-Izz)-Ixx*(thd*phd*cos(ph)-psd*phd*sin(ph)*cos(th)-psd*thd*cos(ph)*sin(th));
    -zcm*Thrustb(1)+R*P*(Izz-Ixx)-Iyy*(-thd*phd*sin(ph)+psd*thd*sin(th)*sin(ph)-psd*phd*cos(th)*cos(ph));
    P*Q*(Ixx-Izz)-Izz*psd*thd*cos(th)];

dx(4:6) = A\b;
end
function dx = EoMfb(t,x,I,zcm,ref,T)
%    ps th ph
% f  1  2  3  
% f' 4  5  6  
dx = zeros(6,1);
ps = x(1);
th = x(2);
ph = x(3);
psd = x(4);
thd = x(5);
phd = x(6);
psg = (ps-ref(1))*0.4-psd*0.1;
thg = (th-ref(2))*0.5-thd*0.1;
T = T(t);

Ixx = I(1);
Iyy = I(2);
Izz = I(3);
P = thd*sin(ph)+psd*cos(ph)*cos(th);
Q = thd*cos(ph)-psd*cos(th)*sin(ph);
R = phd+psd*sin(th);
for ii = 1:3
    dx(ii) = x(ii+3);
end

% Transformations G -> B
bc = [1 0 0;0 cos(psg) -sin(psg);0 sin(psg) cos(psg)];
cg = [cos(thg) 0 sin(thg);0 1 0;-sin(thg) 0 cos(thg)];

Thrustb = bc*cg*[0;0;T];

A = [Ixx*cos(ph)*cos(th) Ixx*sin(ph) 0;
    -Iyy*cos(th)*sin(ph) Iyy*cos(ph) 0;
    Izz*sin(th) 0 Izz];
b = [zcm*Thrustb(2)+Q*R*(Iyy-Izz)-Ixx*(thd*phd*cos(ph)-psd*phd*sin(ph)*cos(th)-psd*thd*cos(ph)*sin(th));
    -zcm*Thrustb(1)+R*P*(Izz-Ixx)-Iyy*(-thd*phd*sin(ph)+psd*thd*sin(th)*sin(ph)-psd*phd*cos(th)*cos(ph));
    P*Q*(Ixx-Izz)-Izz*psd*thd*cos(th)];

dx(4:6) = A\b;
end
%% RK4
function [sol] = RK4(func,tspan,y0,h)
% func = function of x and y
x0 = tspan(1);
xf = tspan(end);
n = floor((xf-x0)/h)+1;
x = zeros(1,n);
y = zeros(length(y0),n);

% Initial conditions
y(:,1) = y0;
x(1) = x0;

for ii = 1:n-1
    x(ii+1) = x0+ii*h;
    y(:,ii+1) = RK4_calc(x(ii),h,y(:,ii));
end

% Hit final time value straight on
if xf > x(end)
    h_end = xf-x(n);
    x(n+1) = xf;
    y(:,n+1) = RK4_calc(x(n),h_end,y(:,n));
end

% Save to solution 
sol.x = x;
sol.y = y;

function yi2 = RK4_calc(xi,dx,yi1)
    k1 = func(xi,yi1);
    k2 = func(xi+.5*dx,yi1+.5*k1*dx);
    k3 = func(xi+.5*dx,yi1+.5*k2*dx);
    k4 = func(xi+dx,yi1+k3*dx);
    yi2 = yi1 + (1/6)*(k1+2*k2+2*k3+k4)*dx; 
end
end