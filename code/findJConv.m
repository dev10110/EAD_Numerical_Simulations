function [Jin, Jout, Jnet, Jdiff] = findJConv(r_c,finterp,mu,ni,center,s,uxinterp,uyinterp)

global do_flow

%% determine Gaussian surface
% for now, let's make it a square

% y_top = 3*d/4;
% y_bot = d/4;
% x_right = d/4;
% x_left = -d/4;

y_top = center(2) + s/2;
y_bot = center(2) - s/2;
x_right = center(1) + s/2;
x_left = center(1) - s/2;

x1 = linspace(x_left,x_right,ni*2);
x2 = x_right*ones(1,ni*2);
x3 = linspace(x_right,x_left,ni*2);
x4 = x_left*ones(1,ni*2);

y1 = y_bot*ones(1,ni*2);
y2 = linspace(y_bot,y_top,ni*2);
y3 = y_top*ones(1,ni*2);
y4 = linspace(y_top,y_bot,ni*2);

x = [x1 x2 x3 x4];
y = [y1 y2 y3 y4];

dA = s/(2*ni-1);

Jnet = 0;
Jin = 0;
Jout = 0;

for i = 1:length(x)-1
    
    % skip redundant corner points
    if x(i+1) == x(i) && y(i+1) == y(i)
        continue
    end
    
    x_mid = (x(i+1) + x(i))/2;
    y_mid = (y(i+1) + y(i))/2;
    
    %% find space charge
    rho = finterp(x_mid,y_mid);
    
    %% find electric field    
    Ex = uxinterp(x_mid,y_mid);
    Ey = uyinterp(x_mid,y_mid);
    
    %% find potential flow contribution
    
        u = 0;
        v = 0;
    
    
    %% determine normal vector
    if y(i+1) == y(i)
        nx = 0;
        if y(i) == y_bot
            ny = -1;
        else
            ny = 1;
        end
    else
        ny = 0;
        if x(i) == x_left
            nx = -1;
        else
            nx = 1;
        end
    end
    
    %% find flux
    dJ = rho*((mu*Ex+u)*nx + (mu*Ey+v)*ny)*dA;
    if dJ<0
        Jin = Jin + dJ;
    else
        Jout = Jout + dJ;
    end    
    Jnet = Jnet + dJ;
    
    Jdiff = abs(Jnet/Jout)*100;
    
end

end

















