%% Main
% 
%% Setup







%% EoM
function dx = EoM(t,x,I,t_burn_start,m,g,zcm,zac)
%    x  y  z  ps th ph
% f  1  2  3  4  5  6
% f' 7  8  9  10 11 12
dx = zeros(12,1);

T = @(t) ((t > t_burn_start) & (t < t_burn_start+3.5))*15;
% Assume 0 gimbal angles
psg = 0;
thg = 0;
% Set aero forces to 0
Abx = 0;
Aby = 0;
Abz = 0;

Ixx = I(1);
Iyy = I(2);
Izz = I(3);
P = x(11)*sin(x(6))+x(10)*cos(x(6))*cos(x(5));
Q = x(11)*cos(x(6))-x(10)*cos(x(5))*sin(x(6));
R = x(12)+x(10)*sin(x(5));
for ii = 1:6
    dx(ii) = x(ii+6);
end
% Transformations B -> N
na1 = [1 0 0;0 cos(x(4)) -sin(x(4));0 sin(x(4)) cos(x(4))];
a1a2 = [cos(x(5)) 0 sin(x(5));0 1 0;-sin(x(5)) 0 cos(x(5))];
a2b = [cos(x(6)) -sin(x(6)) 0;sin(x(6)) cos(x(6)) 0;0 0 1];
% Transformations G -> B
bc = [1 0 0;0 cos(psg) -sin(psg);0 sin(psg) cos(psg)];
cg = [cos(thg) 0 sin(thg);0 1 0;-sin(thg) 0 cos(thg)];

Thrustn = na1*a1a2*a2b*bc*cg*[0;0;T(t)];
Thrustb = na1*a1a2*a2b*[0;0;T(t)];
An = na1*a1a2*a2b*[Abx;Aby;Abz];

dx(7) = (Thrustn(1)+An(1))/m;
dx(8) = (Thrustn(2)+An(2))/m;
dx(9) = (Thrustn(3)+An(3)-m*g)/m;

A = [Ixx*cos(x(6))*cos(x(5)) Ixx*sin(x(6)) 0;
    -Iyy*cos(x(5))*sin(x(6)) Iyy*cos(x(6)) 0;
    Izz*sin(x(5)) 0 Izz];
b = [zcm*Thrustb(2)-zac*Aby+Q*R*(Iyy-Izz)-Ixx*(x(11)*x(12)*cos(x(6))-x(10)*x(12)*sin(x(6))*cos(x(5))-x(10)*x(11)*cos(x(6))*sin(x(5)));
    -zcm*Thrustb(1)+zac*Abx+R*P*(Izz-Ixx)-Iyy*(-x(11)*x(12)*sin(x(6))+x(10)*x(11)*sin(x(5))*sin(x(6))-x(10)*x(12)*cos(x(5))*cos(x(6)));
    P*Q*(Ixx-Izz)-Izz*x(10)*x(11)*cos(x(5))];

dx(10:12) = A\b;
end