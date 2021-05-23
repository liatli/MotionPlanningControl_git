clc
clear all
close all

%%%%%%%%%%%%%%%%
% ADRC control and MNPL control
%%%%%%%%%%%%%%
% omega/X = K/(Ts + 1)
% dt_omega = - 1/T (omega) + K/T X
K = 1.15;
T = 0.5;

% ADRC paramters
b = 2.3;
beta1 = 30;
beta2 = 300;
beta3 = 15;
kp  = 0.6;
kd  = 0.05;


% Control targets
m2pix = 400;

%%%%% Posture Control Parameters %%%%%%%
p0.x = 196/m2pix ;
p0.y = 106/m2pix;
p0.w = -0.27*pi;

pe.x = 400/m2pix;
pe.y = 100/m2pix;
pe.w = pi/3;

% initialisation
k = 0.112;
ba = 2.3;
bp = 2.3;
R = 0;
f = 0.0;

beta1a = 30;
beta2a = 300;
beta3a = 15;
beta1p = 30;
beta2p = 300;
beta3p = 15;
z1a = 0;
z2a = 0;
z3a = 0;
z1p = 0;
z2p = 0;
z3p = 0;
kpa = 0.04;
kda = 0.02;
kpp = 0.1;
kdp = 0.01;
%%%%% End Posture Control Parameters %%%%%%%


dt = 0.1;
p.x = p0.x;
p.y = p0.y;
p.w = p0.w;


theta_v = atan2(pe.y-p0.y,pe.x-p0.x);


z1 = 0;
z2 = 0;
z3 = 0;
X = 0;
u = 0;
u0 = 0;

dbug.y = [];
dbug.u = [];
dbug.X = [];
dbug.ew = [];

for t = 0:dt:50
    
    % solve the model
    if abs(X)<5
    sol=ode45(@Equation2 ,[t ,t + dt] ,...
        [0.12, 0 ,0, f, X,108680.61e-7]);
    else
    sol=ode45(@Equation2 ,[t ,t + dt] ,...
        [0.12, 0 ,0, 1.0, X,108680.61e-7]);     
    end
    ndt = size(sol.y,2);
    for ii = 1:size(sol.y,2)
        p.w(end+1)=p.w(end)+dt/ndt*deval(sol,sol.x(ii),3);
        p.x(end+1)=p.x(end)+dt/ndt*cos(p.w(end))*(deval(sol,sol.x(ii),1) - 0.12)...
            -dt/ndt*sin(p.w(end))*deval(sol,sol.x(ii),2);
        p.y(end+1)=p.y(end)+dt/ndt*sin(p.w(end))*(deval(sol,sol.x(ii),1) - 0.12)...
            +dt/ndt*cos(p.w(end))*deval(sol,sol.x(ii),2);
    end
    
    
    

    nT = 168;
    if length(p.y)<=nT
        nl = length(p.y);
        nl2 = int16(nl/2);
        y = atan2(mean(p.y(end-nl2:end))-mean(p.y(1:end-nl2)),mean(p.x(end-nl2:end))-mean(p.x(1:end-nl2)));
        theta_v = atan2(pe.y-mean(p.y(end-nl2:end)), pe.x-mean(p.x(end-nl2:end)));
    else
        y = atan2(mean(p.y(end-nT/2:end)-p.y(end-nT:end-nT/2)),mean(p.x(end-nT/2:end)-p.x(end-nT:end-nT/2)));
        theta_v = atan2(pe.y-mean(p.y(end-nT/2:end)), pe.x-mean(p.x(end-nT/2:end)));
    end
    
    
    e.x = cos(pe.w)*(p.x(end)-pe.x)+sin(pe.w)*(p.y(end)-pe.y);
    e.y = -sin(pe.w)*(p.x(end)-pe.x)+cos(pe.w)*(p.y(end)-pe.y);
    e.w = y - pe.w;
    
    v = k*sqrt(e.x^2+e.y^2);
    
    ea   = z1a - e.w;
    dz1a = z2a - beta1a*ea;
    dz2a = z3a - beta2a*ea;
    dz3a = -beta3a*ea;
    z1a = z1a + dt*dz1a;
    z2a = z2a + dt*dz2a;
    z3a = z3a + dt*dz3a;
    u0a  =- kpa*(e.w);
    ua   = u0a - z3a/ba;
    
    
    ep   = z1p - e.y;
    dz1p = z2p - beta1p*ep;
    dz2p = z3p - beta2p * ep;
    dz3p = -beta3p * ep;
    z1p = z1p + dt*dz1p;
    z2p = z2p + dt*dz2p;
    z3p = z3p + dt*dz3p;
    u0p  = -kpp*e.y;
    up   = u0p - z3p/bp;
    
    
    u = ua + up;
    
    
    
    P = [ -53.93  -3.265];
    if u > 0.5
        u = 0.5;
    elseif u<-0.5
        u = -0.5;
    end

    X = polyval(P,u);
    
    P = [2.385, 0.2412];
    if v >0.5
        v = 0.5;
    elseif v<0
        v = 0;
    end
    f = polyval(P,v);
        
    
    
    dbug.y(end+1) = y;
    dbug.u(end+1) = u;
    dbug.X(end+1) = X;
    dbug.ew(end+1) = e.w;
    
    
end

plot(p0.x,p0.y,'<')
hold on
plot(pe.x,pe.y,'o')
plot(p.x,p.y)









