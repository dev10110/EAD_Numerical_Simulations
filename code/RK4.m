function [ x,y,t,dxdt,dydt,xtemp,ytemp ] = RK4( x,y,Ex1,Ey1,t,dt,mu,uxinterp,uyinterp,r,forward )
        
%% Step 1
% [Ex1, Ey1] = tri2gridDEX(x,y,uxinterp,uyinterp);
% [ux1, uy1] = getPF( x, y, r, U );

if forward
    kx1 = mu*Ex1 ;
    ky1 = mu*Ey1;
else
    kx1 = -mu*Ex1;
    ky1 = -mu*Ey1;
end

%% Step 2
[Ex2, Ey2] = tri2gridDEX(x + dt/2*kx1,y + dt/2*ky1,uxinterp,uyinterp);

    ux2 = 0;
    uy2 = 0;


if forward
    kx2 = mu*Ex2 + ux2;
    ky2 = mu*Ey2 + uy2;
else
    kx2 = -mu*Ex2 - ux2;
    ky2 = -mu*Ey2 - uy2;
end

%% Step 3
[Ex3, Ey3] = tri2gridDEX(x + dt/2*kx2,y + dt/2*ky2,uxinterp,uyinterp);

    ux3 = 0;
    uy3 = 0;

if forward
    kx3 = mu*Ex3 + ux3;
    ky3 = mu*Ey3 + uy3;
else
    kx3 = -mu*Ex3 - ux3;
    ky3 = -mu*Ey3 - uy3;
end

%% Step 4
[Ex4, Ey4] = tri2gridDEX(x + dt*kx3,y + dt*ky3,uxinterp,uyinterp);

    ux4 = 0;
    uy4 = 0;

if forward
    kx4 = mu*Ex4 + ux4;
    ky4 = mu*Ey4 + uy4;
else
    kx4 = -mu*Ex4 - ux4;
    ky4 = -mu*Ey4 - uy4;
end

%% Final output
xtemp = x + 1/6*dt*(kx1 + 2*kx2 + 2*kx3 + kx4);
ytemp = y + 1/6*dt*(ky1 + 2*ky2 + 2*ky3 + ky4);
t = t + dt;

dxdt = (xtemp-x)/dt;
dydt = (ytemp-y)/dt;

x = xtemp;
y = ytemp;
end

