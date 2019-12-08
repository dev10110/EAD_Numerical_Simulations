function [ x,y,t,dxdt,dydt ] = Eulerf( x,y,Ex,Ey,t,dt,mu,ux,uy )
        
dxdt = mu*Ex + ux;
dydt = mu*Ey + uy;

x = x + dxdt*dt;
y = y + dydt*dt;

t = t + dt;

end

