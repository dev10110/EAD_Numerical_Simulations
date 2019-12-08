function [ Ex,Ey ] = getE_alt( inds, temp )

% if isempty(uxinterp)   
%     
%     uxinterp = TriScatteredInterp(p(1,:)',p(2,:)',-temp(:,1));
%     uyinterp = TriScatteredInterp(p(1,:)',p(2,:)',-temp(:,2));
%     
% end
% 
% [Ex, Ey] = tri2gridDEX(p,t,-temp,xe,ye,uxinterp,uyinterp);

Ex = -temp(1,inds);
Ey = -temp(2,inds);

end