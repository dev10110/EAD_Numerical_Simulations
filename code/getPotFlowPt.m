function [ u, v ] = getPotFlowPt(r_c,U,x,y)

%% set constants

% set problem geometry in bipolar coordinates
xi2_val = 2.82017;
xi1_val = -6.26456;
a = sqrt(r_c^2*sinh(xi2_val)^2);

%% initialize bipolar variables

xi_val = acoth((a^2 + x^2 + y^2)/(2*a*x));
eta_val = acot((-a^2 + x^2 + y^2)/(2*a*y));

if xi_val < xi1_val
    xi_val = xi1_val;
end
if xi_val > xi2_val
    xi_val = xi2_val;
end

num_n = 7;

total_dxi = 0;
total_deta = 0;

%% term 1

temp_dxi = -2*sin(eta_val)*sinh(xi_val)/(2*cosh(xi_val) - 2*cos(eta_val))^2;
temp_deta = cos(eta_val)/( 2*cosh(xi_val) - 2*cos(eta_val)) ...
    - 2*sin(eta_val)^2/(2*cosh(xi_val) - 2*cos(eta_val))^2;

total_dxi = total_dxi + temp_dxi;
total_deta = total_deta + temp_deta;

%% term 2

temp_dxi = 0;
temp_deta = 0;
for k = 1:num_n
    
    temp_dxi = temp_dxi - cosh((xi_val - xi1_val)*k)*k*exp(-xi2_val*k)*sin(eta_val*k)/sinh((xi2_val - xi1_val)*k);
    temp_deta = temp_deta - sinh((xi_val - xi1_val)*k)*exp(-xi2_val*k)*cos(eta_val*k)*k/sinh((xi2_val - xi1_val)*k);
    
end

total_dxi = total_dxi + temp_dxi;
total_deta = total_deta + temp_deta;

%% term 3

temp_dxi = 0;
temp_deta = 0;
for k = 1:num_n
    
    temp_dxi = temp_dxi + cosh((xi_val - xi2_val)*k)*k*exp(xi1_val*k)*sin(eta_val*k)/sinh((xi2_val - xi1_val)*k);
    temp_deta = temp_deta + sinh((xi_val - xi2_val)*k)*exp(xi1_val*k)*cos(eta_val*k)*k/sinh((xi2_val - xi1_val)*k);
    
end

total_dxi = total_dxi + temp_dxi;
total_deta = total_deta + temp_deta;

%% determine transformation matrix

temp_dxdxi = a*cosh(xi_val)/(cosh(xi_val)-cos(eta_val)) ...
    - a*sinh(xi_val)^2/(cosh(xi_val)-cos(eta_val))^2;
temp_dydxi = -a*sin(eta_val)*sinh(xi_val)/(cosh(xi_val)-cos(eta_val))^2;
temp_dxdeta = -a*sin(eta_val)*sinh(xi_val)/(cosh(xi_val)-cos(eta_val))^2;
temp_dydeta = a*cos(eta_val)/(cosh(xi_val)-cos(eta_val)) ...
    - a*sin(eta_val)^2/(cosh(xi_val)-cos(eta_val))^2;

A = [temp_dxdxi,temp_dydxi;temp_dxdeta,temp_dydeta];
b = [total_dxi;total_deta];
c = A\b;

u = 2*a*U*c(2);
v = -2*a*U*c(1);


end




























