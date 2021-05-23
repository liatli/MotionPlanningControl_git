function dy=Equation2(t,y)
dy=zeros(5,1);

f= y(4);

mb=0.800;
e=1-(0.06/(0.19+0.058+0.050+0.152))^2;
a=2*(1-e^2)*(0.5*log((1+e)/(1-e))-e)/(e^3);
b=1/(e^2)-0.5*(1-e^2)*log((1+e)/(1-e))/(e^3);
X=-mb*a/(2-a);
Y=-mb*b/(2-b);
N=-0.2*mb*((0.03^2-0.225^2)^2)*(a-b)/(2*(0.03^2-0.225^2)+(0.03^2+0.225^2)*(b-a));

x_set = y(5);
PI=deg2rad(180);
X1=deg2rad(x_set);
X2=deg2rad(x_set);
X3=deg2rad(x_set);
R1=deg2rad(5);
R2=deg2rad(10);
R3=deg2rad(15);
phi1=1.396;
phi2=2.094;

L1=0.058;
L2=0.050;
L3=0.052;
L4=0.152;

d=0.067;
D=0.535;

c=0.044;
p=1000;
s=0.02668582;
% j3=388680.61e-9-N;
% j3=388680.61e-8-N;  % the ratio will be lower
% j3=388680.61e-10-N;  % the ratio will be lower
% j3=158680.61e-7-N;  % the ratio will be lower
j3= y(6) - N; % Orientation
m1=mb-X;
m2=mb-Y;

C1=sign(y(1))*0.4357/abs(y(1))^0.4;
C2=45;
K=0.0078;
c1=1/2*p*s*C1;
c2=1/2*p*s*C2;
c4=K/j3;

syms l;
m=1/4*PI*p*d*d;
mm=1/4*PI*p*D*D;
mmm=1/4*PI*p*(0.7*(l-0.052)+0.065)^2;

th1=X1+R1*sin(2*PI*t*f);
th2=X2+R2*sin(2*PI*t*f+phi1);
th3=X3+R3*sin(2*PI*t*f+phi2);

dth1=2*PI*R1*cos(2*PI*t*f)*f;
dth2=2*PI*R2*cos(2*PI*t*f+phi1)*f;
dth3=2*PI*R3*cos(2*PI*t*f+phi2)*f;

ddth1=-2*PI*2*PI*R1*sin(2*PI*t*f)*f*f;
ddth2=-2*PI*2*PI*R2*sin(2*PI*t*f+phi1)*f*f;
ddth3=-2*PI*2*PI*R3*sin(2*PI*t*f+phi2)*f*f;

dw1=l*ddth1;
dw2=l*ddth2+L1*ddth1*cos(th2)-L1*dth1*dth2*sin(th2);
dw3=l*ddth3-dth3*dth2*L2*sin(th3)+ddth2*L2*cos(th3)-dth1*L1*(dth2+dth3)*sin(th2+th3)+L1*ddth1*cos(th2+th3);
dw4=(L1*dth1*cos(th2+th3)+L2*dth2*cos(th3)+L3*dth3)^2;
dw5=(L1*dth1*sin(th2+th3)+L2*dth2*sin(th3))*(L1*dth1*cos(th2+th3)+L2*dth2*cos(th3)+L3*dth3);

f1=-m*(int(dw1*sin(th1),0,L1)+int(dw2*sin(th1+th2),0,L2)+int(dw3*sin(th1+th2+th3),0,L3)+...
    int(mmm*dw3*sin(th1+th2+th3),L3,L4)/m)+0.5*mm*cos(th1+th2+th3)*dw4+mm*sin(th1+th2+th3)*dw5;
f2=-m*(int(dw1*cos(th1),0,L1)+int(dw2*cos(th1+th2),0,L2)+int(dw3*cos(th1+th2+th3),0,L3)+...
    int(mmm*dw3*cos(th1+th2+th3),L3,L4)/m)-0.5*mm*sin(th1+th2+th3)*dw4+mm*cos(th1+th2+th3)*dw5;
f3=-m*(int(dw1*cos(th1)*(c+l*cos(th1))+dw1*sin(th1)*l*sin(th1),0,L1)+...
    int(dw2*cos(th1+th2)*(c+L1*cos(th1)+l*cos(th1+th2))+dw2*sin(th1+th2)*(L1*sin(th1)+l*sin(th1+th2)),0,L2)+...
    int(dw3*cos(th1+th2+th3)*(c+L1*cos(th1)+L2*cos(th1+th2)+l*cos(th1+th2+th3))+...
        dw3*sin(th1+th2+th3)*(l*sin(th1+th2+th3)+L1*sin(th1)+L2*sin(th1+th2)),0,L3)+...
    int(mmm*dw3*cos(th1+th2+th3)*(c+L1*cos(th1)+L2*cos(th1+th2)+l*cos(th1+th2+th3))+...
        mmm*dw3*sin(th1+th2+th3)*(l*sin(th1+th2+th3)+L1*sin(th1)+L2*sin(th1+th2)),L3,L4)/m)+...
0.5*mm*(-sin(th1+th2+th3)*(c+L1*cos(th1)+L2*cos(th1+th2)+L3*cos(th1+th2+th3))+...
    cos(th1+th2+th3)*(L1*sin(th1)+L2*sin(th1+th2)+L3*sin(th1+th2+th3)))*dw4+...
mm*(cos(th1+th2+th3)*(c+L1*cos(th1)+L2*cos(th1+th2)+L3*cos(th1+th2+th3))+...
    sin(th1+th2+th3)*(L1*sin(th1)+L2*sin(th1+th2)+L3*sin(th1+th2+th3)))*dw5;


dy(1) = m2/m1*y(2)*y(3)-c1/m1*y(1)*y(1)+c2/m1*y(2)*y(2)+f1/m1;
dy(2) = -m1/m2*y(1)*y(3)-(c1+c2)/m2*y(1)*y(2)+f2/m2;
dy(3) = (m1-m2)/j3*y(1)*y(2)-c4*y(3)*y(3)*sign(y(3))+f3/j3;
dy(4) = 0;
dy(5) = 0;
dy(6) = 0;




