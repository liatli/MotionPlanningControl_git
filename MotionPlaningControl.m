clc
clear all
close all

%%%%%%%%%%%%%%%%
% ADRC control and MNPL control
%%%%%%%%%%%%%%
K = 1.15;
T = 0.5;

% ADRC parameters:
b = 2.3;
beta1 = 30;
beta2 = 300;
beta3 = 15;
kp  = 0.6;
kd  = 0.05;


% Control targets
m2pix = 400;

% Initial position
p0.x = 150/m2pix ;
p0.y = 125/m2pix;
p0.w = -(1/3)*pi;

% target position
pe.x = 600/m2pix;
pe.y = 350/m2pix;


% initialisation
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

for t = 0:dt:80
    
    % solve the model
    sol=ode45(@Equation2 ,[t ,t + dt] ,...
        [0.12, 0 ,0, 1.0, X,108680.61e-7]);
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
    
    
    e = z1 - y;
    dz1 = z2 - beta1*e;
    dz2 = z3 - beta2*e + b*u;
    dz3 = -beta3*e;
    z1 = z1 + dt*dz1;
    z2 = z2 + dt*dz2;
    z3 = z3 + dt*dz3;
    
    %%%%%%%%%%%%%%%% Peridic the angle differences %%%%%%%%%%%%%%%
    dy_theta =  theta_v - y;    
    u0 = kp*(dy_theta);
    u = u0 - z3/b;
 
    
    P = [ -53.93  -3.265];
    if u > 0.5
        u = 0.5;
    elseif u<-0.5
        u = -0.5;
    end

    X = polyval(P,u);
    
    dbug.y(end+1) = y;
    dbug.u(end+1) = u;
    dbug.X(end+1) = X;

    
end

plot(p0.x,p0.y,'<')
hold on
plot(pe.x,pe.y,'o')
plot(p.x,p.y)









