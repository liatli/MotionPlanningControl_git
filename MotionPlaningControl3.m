clc
clear all
close all

%%%%%%%%%%%%%%%%
% ADRC control and MNPL control 3:  target pose control with a temporal point
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

%%%%% Posture Control Parameters %%%%%%%

p0.x = 608/m2pix ;
p0.y = 394/m2pix;
p0.w = -0.056*pi;

pe.x = 200/m2pix;
pe.y = 300/m2pix;
pe.w = 2*pi/3;

pt.x = 275/m2pix;
pt.y = 170/m2pix;


% initialisation

f = 1.0;

beta1 = 30;
beta2 = 300;
beta3 = 15;
kp  = 0.96;

k = 0.112;
ba = 2.3;
bp = 2.3;
b = 2.3;
R = 0;
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
kpa = 0.9;
kda = 0.01;
kpp = 8;
kdp = 0.05;
v = 0.1;
Flag = 'Point';
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

for t = 0:dt:90
    
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
        theta_v = atan2(pt.y-mean(p.y(end-nl2:end)), pt.x-mean(p.x(end-nl2:end)));
    else
        y = atan2(mean(p.y(end-nT/2:end)-p.y(end-nT:end-nT/2)),mean(p.x(end-nT/2:end)-p.x(end-nT:end-nT/2)));
        theta_v = atan2(pt.y-mean(p.y(end-nT/2:end)), pt.x-mean(p.x(end-nT/2:end)));
    end
    
    
    %%%%% Target Pose Control With Temp Point%%%%%%%
    
    if norm([p.x(end)-pt.x,p.y(end)-pt.y]) > 0.07 && strcmp(Flag,'Point')
        
        e1 = z1 - y;
        dz1 = z2 - beta1*e1;
        dz2 = z3 - beta2*e1 + b*u;
        dz3 = -beta3*e1;
        z1 = z1 + dt*dz1;
        z2 = z2 + dt*dz2;
        z3 = z3 + dt*dz3;
        
        %%%%%%%%%%%%%%%% Peridic the angle differences %%%%%%%%%%%%%%%
        dy_theta =  theta_v - y;
        if dy_theta > pi
            while (dy_theta > pi)
                dy_theta = dy_theta - 2*pi ;
            end
        end
        
        if dy_theta < -pi
            while (dy_theta < -pi)
                dy_theta = dy_theta + 2*pi ;
            end
        end
        
        u0 = kp*(dy_theta);
        u = u0 - z3/b;
    else
        Flag = 'Pose';
        theta_v = atan2(pe.y-mean(p.y(end-nT/2:end)), pe.x-mean(p.x(end-nT/2:end)));
        e.x = cos(pe.w)*(p.x(end)-pe.x)+sin(pe.w)*(p.y(end)-pe.y);
        e.y = -sin(pe.w)*(p.x(end)-pe.x)+cos(pe.w)*(p.y(end)-pe.y);
        e.w = y - pe.w;
        if e.w > pi
            while (e.w > pi)
                e.w = e.w - 2*pi ;
            end
        end
        if e.w < -pi
            while (e.w < -pi)
                e.w = e.w + 2*pi ;
            end
        end
        
        v = k*sqrt(e.x^2+e.y^2);
        
        ea   = z1a - e.w;
        if ea > pi
            while (ea > pi)
                ea = ea - 2*pi ;
            end
        end
        
        if ea < -pi
            while (ea < -pi)
                ea = ea + 2*pi ;
            end
        end
        dz1a = z2a - beta1a*ea;
        dz2a = z3a - beta2a*ea;
        dz3a = -beta3a*ea;
        z1a = z1a + dt * dz1a;
        z2a = z2a + dt * dz2a;
        z3a = z3a + dt * dz3a;
        u0a  = -kpa*(e.w);
        ua   = u0a - z3a/ba;
        
        
        
        ep   = z1p - e.y;
        dz1p = z2p - beta1p*ep;
        dz2p = z3p - beta2p * ep;
        dz3p = -beta3p * ep;
        z1p = z1p + dt * dz1p;
        z2p = z2p + dt * dz2p;
        z3p = z3p + dt * dz3p;
        u0p  = -kpp*e.y;
        up   = u0p - z3p/bp;
        
        u = ua + up;
    end
   
    
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

    
    
    
end

plot(p0.x,p0.y,'<')
hold on
plot(pt.x,pt.y,'*')
plot(pe.x,pe.y,'o')
plot(p.x,p.y)









