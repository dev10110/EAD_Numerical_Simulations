function [ f_new,F] = interpf( p,t,rho_store,epsi)

% Triangle point indices
it1=t(1,:);
it2=t(2,:);
it3=t(3,:);

% Find centroids of triangles
xpts=(p(1,it1)+p(1,it2)+p(1,it3))/3;
ypts=(p(2,it1)+p(2,it2)+p(2,it3))/3;

% get all query points to define f
xq = xpts;
yq = ypts;

% use built-in interpolation MATLAB function
F = scatteredInterpolant(rho_store(1,:)',rho_store(2,:)',rho_store(3,:)','natural');
f_new = F(xq,yq)/epsi; %since rho actually the space charge, f_new is the right side of Poissons Eq, rho/eps0

end

