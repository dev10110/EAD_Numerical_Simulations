function [ uxinterp,uyinterp ] = getEInterp( p,t,temp,r_c)

% Triangle point indices
it1=t(1,:);
it2=t(2,:);
it3=t(3,:);

% Find centroids of triangles
xpts=(p(1,it1)+p(1,it2)+p(1,it3))/3;
ypts=(p(2,it1)+p(2,it2)+p(2,it3))/3;

% only need half the points
if r_c == 'All'
    temp1 = temp(1,:);
    temp2 = temp(2,:);
else
    inds = xpts<=2*r_c;
    xpts = xpts(inds);
    ypts = ypts(inds);
    temp1 = temp(1,inds);
    temp2 = temp(2,inds);    
end


temp = [temp1;temp2];

uxinterp = scatteredInterpolant(xpts',ypts',-temp(1,:)','natural');
uyinterp = scatteredInterpolant(xpts',ypts',-temp(2,:)','natural');

end

