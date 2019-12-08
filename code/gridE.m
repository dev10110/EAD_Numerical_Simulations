function [ Ex_grid,Ey_grid ] = gridE( ng, meshdomain,p,t,temp )

x = linspace(meshdomain(1),meshdomain(2),ng);
y = linspace(meshdomain(3),meshdomain(4),ng);
Ex_grid = tri2grid(p,t,-temp(:,1),x,y);
Ey_grid = tri2grid(p,t,-temp(:,2),x,y);

end

