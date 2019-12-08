function [ u_interp, v_interp ] = getPotFlow(r_c,d,meshdomain,U)

syms xi eta n xi2 xi1 

%% set constants

% set problem geometry in bipolar coordinates
xi2_val = 2.82017;
xi1_val = -6.26456;
a = sqrt(r_c^2*sinh(xi2_val)^2);

%% initialize bipolar variables

[xi_val,eta_val] = meshgrid(linspace(xi1_val,xi2_val,100),linspace(0,2*pi(),20));

num_n = 7;

total_dxi = zeros(size(xi_val,1),size(xi_val,2));
total_deta = zeros(size(xi_val,1),size(xi_val,2));
x_store = zeros(size(xi_val,1),size(xi_val,2));
y_store = zeros(size(xi_val,1),size(xi_val,2));
u_store = zeros(size(xi_val,1),size(xi_val,2));
v_store = zeros(size(xi_val,1),size(xi_val,2));

%% term 1

for i = 1:size(xi_val,1)
    for j = 1:size(xi_val,2)
        
        temp_dxi = -2*sin(eta_val(i,j))*sinh(xi_val(i,j))/(2*cosh(xi_val(i,j)) - 2*cos(eta_val(i,j)))^2;
        temp_deta = cos(eta_val(i,j))/( 2*cosh(xi_val(i,j)) - 2*cos(eta_val(i,j))) ...
            - 2*sin(eta_val(i,j))^2/(2*cosh(xi_val(i,j)) - 2*cos(eta_val(i,j)))^2;
        
        total_dxi(i,j) = total_dxi(i,j) + temp_dxi;
        total_deta(i,j) = total_deta(i,j) + temp_deta;
        
    end
end



%% term 2

for i = 1:size(xi_val,1)
    for j = 1:size(xi_val,2)
        
        temp_dxi = 0;
        temp_deta = 0;
        for k = 1:num_n
            
            temp_dxi = temp_dxi - cosh((xi_val(i,j) - xi1_val)*k)*k*exp(-xi2_val*k)*sin(eta_val(i,j)*k)/sinh((xi2_val - xi1_val)*k);
            temp_deta = temp_deta - sinh((xi_val(i,j) - xi1_val)*k)*exp(-xi2_val*k)*cos(eta_val(i,j)*k)*k/sinh((xi2_val - xi1_val)*k);
            
        end
        
        total_dxi(i,j) = total_dxi(i,j) + temp_dxi;
        total_deta(i,j) = total_deta(i,j) + temp_deta;
        
    end
end

%% term 3

for i = 1:size(xi_val,1)
    for j = 1:size(xi_val,2)
        
        temp_dxi = 0;
        temp_deta = 0;
        for k = 1:num_n

            temp_dxi = temp_dxi + cosh((xi_val(i,j) - xi2_val)*k)*k*exp(xi1_val*k)*sin(eta_val(i,j)*k)/sinh((xi2_val - xi1_val)*k);
            temp_deta = temp_deta + sinh((xi_val(i,j) - xi2_val)*k)*exp(xi1_val*k)*cos(eta_val(i,j)*k)*k/sinh((xi2_val - xi1_val)*k);
            
        end
        
        total_dxi(i,j) = total_dxi(i,j) + temp_dxi;
        total_deta(i,j) = total_deta(i,j) + temp_deta;
        
    end
end

%% determine transformation matrix

x = a*sinh(xi)/(cosh(xi)-cos(eta));
y = a*sin(eta)/(cosh(xi)-cos(eta));

for i = 1:size(xi_val,1)
    for j = 1:size(xi_val,2)
        
        temp_dxdxi = a*cosh(xi_val(i,j))/(cosh(xi_val(i,j))-cos(eta_val(i,j))) ...
            - a*sinh(xi_val(i,j))^2/(cosh(xi_val(i,j))-cos(eta_val(i,j)))^2;
        temp_dydxi = -a*sin(eta_val(i,j))*sinh(xi_val(i,j))/(cosh(xi_val(i,j))-cos(eta_val(i,j)))^2;
        temp_dxdeta = -a*sin(eta_val(i,j))*sinh(xi_val(i,j))/(cosh(xi_val(i,j))-cos(eta_val(i,j)))^2;
        temp_dydeta = a*cos(eta_val(i,j))/(cosh(xi_val(i,j))-cos(eta_val(i,j))) ...
            - a*sin(eta_val(i,j))^2/(cosh(xi_val(i,j))-cos(eta_val(i,j)))^2;
        
        A = [temp_dxdxi,temp_dydxi;temp_dxdeta,temp_dydeta];
        b = [total_dxi(i,j);total_deta(i,j)];
        c = A\b;
        
        u_store(i,j) = 2*a*U*c(2);
        v_store(i,j) = -2*a*U*c(1);
        x_store(i,j) = a*sinh(xi_val(i,j))/(cosh(xi_val(i,j))-cos(eta_val(i,j)));
        y_store(i,j) = a*sin(eta_val(i,j))/(cosh(xi_val(i,j))-cos(eta_val(i,j)));
        
    end
end

%% build interpolation surface

count = 1;

for i = 1:size(xi_val,1)
    for j = 1:size(xi_val,2)
        
        x_all(count) = x_store(i,j);
        y_all(count) = y_store(i,j);
        u_all(count) = u_store(i,j);
        v_all(count) = v_store(i,j);
        
        count = count + 1;
        
    end
end

last_count = count;

% add in boundary values
ni = 50;

x_min = meshdomain(1);
x_max = meshdomain(2);
y_min = meshdomain(3);
y_max = meshdomain(4);

x1 = linspace(x_min,x_max,ni*2);
x2 = x_max*ones(1,ni*2);
x3 = linspace(x_min,x_max,ni*2);
x4 = x_min*ones(1,ni*2);

y1 = y_min*ones(1,ni*2);
y2 = linspace(y_min,y_max,ni*2);
y3 = y_max*ones(1,ni*2);
y4 = linspace(y_min,y_max,ni*2);

x = [x1, x2(2:end), x3(1:end-1), x4(2:end-1)];
y = [y1, y2(2:end), y3(1:end-1), y4(2:end-1)];

x_all(last_count:last_count + length(x) - 1) = x;
y_all(last_count:last_count + length(y) - 1) = y;
u_all(last_count:last_count + length(x) - 1) = U;
v_all(last_count:last_count + length(x) - 1) = 0;

u_interp = TriScatteredInterp(x_all',y_all',u_all');
v_interp = TriScatteredInterp(x_all',y_all',v_all'); 

eval(['save flow_interp_U', num2str(U),'_d',num2str(100*d),' u_interp v_interp'])

end




























