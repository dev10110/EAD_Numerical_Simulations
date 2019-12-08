function [ uxinterp,uyinterp ] = getEInterpFull( p,t,temp)

% Triangle point indices
it1=t(1,:);
it2=t(2,:);
it3=t(3,:);

% Find centroids of triangles
xpts=(p(1,it1)+p(1,it2)+p(1,it3))/3;
ypts=(p(2,it1)+p(2,it2)+p(2,it3))/3;

uxinterp = scatteredInterpolant(xpts',ypts',-temp(1,:)');
uyinterp = scatteredInterpolant(xpts',ypts',-temp(2,:)');

end

