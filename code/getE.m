function [ Ex,Ey ] = getE( xe, ye, temp, p, t, uxinterp,uyinterp )

if isempty(uxinterp)   
    
    uxinterp = TriScatteredInterp(p(1,:)',p(2,:)',-temp(:,1));
    uyinterp = TriScatteredInterp(p(1,:)',p(2,:)',-temp(:,2));
    
end

[Ex, Ey] = tri2gridDEX(p,t,-temp,xe,ye,uxinterp,uyinterp);

end