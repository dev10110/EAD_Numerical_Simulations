function [ dA ] = getdAFEM( p,t )

% Triangle point indices
it1=t(1,:);
it2=t(2,:);
it3=t(3,:);

% Find centroids of triangles
x1 = p(1,it1);
x2 = p(1,it2);
x3 = p(1,it3);
y1 = p(2,it1);
y2 = p(2,it2);
y3 = p(2,it3);

dA = abs(x1.*(y2 - y3) + x2.*(y3 - y1) + x3.*(y1 - y2))/2;

end

