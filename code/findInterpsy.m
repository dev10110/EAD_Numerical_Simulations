function [ xi, yi ] = findInterpsy( ni, d, meshdomain, x_max )

pert = 0.000001;
offset = 3;

y_temp = linspace(meshdomain(3)+pert,meshdomain(4) - pert,ni*4);
x_temp = 3/4*meshdomain(1);

xi = x_temp;
yi = y_temp((1+2*offset):(end-2*offset));

end

