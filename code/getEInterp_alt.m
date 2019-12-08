function [ uxinterp,uyinterp,interp_store ] = getEInterp_alt( p,t,temp,r_c,npts,meshdomain)

% Triangle point indices
it1=t(1,:);
it2=t(2,:);
it3=t(3,:);

% Find centroids of triangles
xpts=(p(1,it1)+p(1,it2)+p(1,it3))/3;
ypts=(p(2,it1)+p(2,it2)+p(2,it3))/3;

% only need half the points
inds = xpts<=2*r_c;
xpts = xpts(inds);
ypts = ypts(inds);
temp1 = temp(1,inds);
temp2 = temp(2,inds);

temp = [temp1;temp2];

% split points into grids so as to reduce size of overall interpolating
% surface

xmin = meshdomain(1);
xmax = meshdomain(2);
ymin = meshdomain(3);
ymax = meshdomain(4);

x = linspace(xmin,xmax,npts+1);
y = linspace(ymin,ymax,npts+1);

for i = 1:npts
    for j = 1:npts
        count = 1;
        for k = 1:length(xpts)
            if xpts(k) <= x(i+1) && xpts(k) >= x(i) && ...
                    ypts(k) <= y(j+1) && ypts(k) >= y(j)
                xtemp(count) = xpts(k);
                ytemp(count) = ypts(k);
                temptemp(:,count) = temp(:,k);
                count = count + 1;
            end
        end
        name1 = ['ux',num2str(i),num2str(j)];
        name2 = ['uy',num2str(i),num2str(j)];
        interp_store.(name1) = scatteredInterpolant(xtemp',ytemp',-temptemp(1,:)');
        interp_store.(name2) = scatteredInterpolant(xtemp',ytemp',-temptemp(2,:)');
    end
end

uxinterp = scatteredInterpolant(xpts',ypts',-temp(1,:)');
uyinterp = scatteredInterpolant(xpts',ypts',-temp(2,:)');

end

