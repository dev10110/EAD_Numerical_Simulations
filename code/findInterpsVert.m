function [ xi, yi ] = findInterpsVert( ni,meshdomain,x,y )

pert = 0.000001;

y_temp = linspace(y,meshdomain(4)-pert,ni*3);
x_temp = x;

xi = x_temp;
yi = y_temp(2:end-1);

end

