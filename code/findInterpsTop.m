function [ xi, yi ] = findInterpsTop( ni,x,y )

pert = 0.000001;

y_temp = y; 
x_temp = linspace(x,0-pert,ni*4);

xi = x_temp(2:end);
yi = y_temp;

end

