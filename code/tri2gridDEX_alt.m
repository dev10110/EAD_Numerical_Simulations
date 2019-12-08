function [ uxi, uyi ] = tri2gridDEX_alt(xc,yc,interp_store,meshdomain,varargin)

% uxi = varargin{1}(x,y);
% uyi = varargin{2}(x,y);

xmin = meshdomain(1);
xmax = meshdomain(2);
ymin = meshdomain(3);
ymax = meshdomain(4);

x = linspace(xmin,xmax,20);
y = linspace(ymin,ymax,20);

for i = 1:19
    for j = 1:19
        if xc <= x(i+1) && xc >= x(i) && ...
                yc <= y(j+1) && yc >= y(j)
            name1 = ['ux',num2str(i),num2str(j)];
            name2 = ['uy',num2str(i),num2str(j)];
            uxtemp = interp_store.(name1);
            uytemp = interp_store.(name2);
            uxi = uxtemp(xc,yc);
            uyi = uytemp(xc,yc);
        end
    end
end

end