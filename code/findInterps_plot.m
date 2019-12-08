function [ xi, yi ] = findInterps_plot( ni, d, meshdomain, x_max )

pert = 0.000001;
offset = 3;

y_temp = d/2; 
x_temp = linspace(meshdomain(1),-x_max-pert,ni*4);

xi = x_temp(1+5*offset:end);
yi = y_temp;

end

