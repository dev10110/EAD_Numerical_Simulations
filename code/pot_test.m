clear
clc;
close all;

r_c = 0.00635/2;
r_e = 0.202/2*10^-3;
d = 0.05;

U = 10;

syms xi eta n xi2 xi1 

%% set constants

% set problem geometry in bipolar coordinates
xi2_val = 2.82017;
xi1_val = -6.26456;
a = sqrt(r_c^2*sinh(xi2_val)^2);

%% initialize bipolar variables

%[xi_val,eta_val] = meshgrid(linspace(xi1_val,xi2_val,100),linspace(0,2*pi(),20));
[xi_val,eta_val] = meshgrid(xi2_val,linspace(0,2*pi(),20));

num_n = 7;

total_dxi = zeros(size(xi_val,1),size(xi_val,2));
total_deta = zeros(size(xi_val,1),size(xi_val,2));
x_store = zeros(size(xi_val,1),size(xi_val,2));
y_store = zeros(size(xi_val,1),size(xi_val,2));
u_store = zeros(size(xi_val,1),size(xi_val,2));
v_store = zeros(size(xi_val,1),size(xi_val,2));

%% term 1

t1 = sin(eta)/(2*(cosh(xi)-cos(eta)));
t1_dxi = diff(t1,xi);
t1_deta = diff(t1,eta);


for i = 1:size(xi_val,1)
    for j = 1:size(xi_val,2)
        
        %         temp_dxi = eval(subs(t1_dxi,[xi,eta],[xi_val(i,j),eta_val(i,j)]));
        %         temp_deta = eval(subs(t1_deta,[xi,eta],[xi_val(i,j),eta_val(i,j)]));
        
        temp_dxi = -2*sin(eta_val(i,j))*sinh(xi_val(i,j))/(2*cosh(xi_val(i,j)) - 2*cos(eta_val(i,j)))^2;
        temp_deta = cos(eta_val(i,j))/( 2*cosh(xi_val(i,j)) - 2*cos(eta_val(i,j))) ...
            - 2*sin(eta_val(i,j))^2/(2*cosh(xi_val(i,j)) - 2*cos(eta_val(i,j)))^2;
        
        total_dxi(i,j) = total_dxi(i,j) + temp_dxi;
        total_deta(i,j) = total_deta(i,j) + temp_deta;
        
    end
end



%% term 2

t2 = -sin(n*eta)*exp(-n*xi2)*sinh(n*(xi - xi1))/sinh(n*(xi2 - xi1));
t2_dxi = diff(t2,xi);
t2_deta = diff(t2,eta);

for i = 1:size(xi_val,1)
    for j = 1:size(xi_val,2)
        
        temp_dxi = 0;
        temp_deta = 0;
        for k = 1:num_n
            
            %             temp_dxi = temp_dxi + eval(subs(t2_dxi,[xi,eta,n,xi1,xi2],[xi_val(i,j),eta_val(i,j),k,xi1_val,xi2_val]));
            %             temp_deta = temp_deta + eval(subs(t2_deta,[xi,eta,n,xi1,xi2],[xi_val(i,j),eta_val(i,j),k,xi1_val,xi2_val]));
            %
            
            temp_dxi = temp_dxi - cosh((xi_val(i,j) - xi1_val)*k)*k*exp(-xi2_val*k)*sin(eta_val(i,j)*k)/sinh((xi2_val - xi1_val)*k);
            temp_deta = temp_deta - sinh((xi_val(i,j) - xi1_val)*k)*exp(-xi2_val*k)*cos(eta_val(i,j)*k)*k/sinh((xi2_val - xi1_val)*k);
            
        end
        
        total_dxi(i,j) = total_dxi(i,j) + temp_dxi;
        total_deta(i,j) = total_deta(i,j) + temp_deta;
        
    end
end

temp_dxi
temp_deta

%% term 3

t3 = sin(n*eta)*exp(n*xi1)*sinh(n*(xi - xi2))/sinh(n*(xi2 - xi1));
t3_dxi = diff(t3,xi);
t3_deta = diff(t3,eta);

for i = 1:size(xi_val,1)
    for j = 1:size(xi_val,2)
        
        temp_dxi = 0;
        temp_deta = 0;
        for k = 1:num_n
            
            %             temp_dxi = temp_dxi + eval(subs(t3_dxi,[xi,eta,n,xi1,xi2],[xi_val(i,j),eta_val(i,j),k,xi1_val,xi2_val]));
            %             temp_deta = temp_deta + eval(subs(t3_deta,[xi,eta,n,xi1,xi2],[xi_val(i,j),eta_val(i,j),k,xi1_val,xi2_val]));
            %
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

dxdxi = diff(x,xi);
dydxi = diff(y,xi);
dxdeta = diff(x,eta);
dydeta = diff(y,eta);

for i = 1:size(xi_val,1)
    for j = 1:size(xi_val,2)
        
        %         temp_dxdxi = eval(subs(dxdxi,[xi,eta,n,xi1,xi2],[xi_val(i,j),eta_val(i,j),k,xi1_val,xi2_val]));
        %         temp_dydxi = eval(subs(dydxi,[xi,eta,n,xi1,xi2],[xi_val(i,j),eta_val(i,j),k,xi1_val,xi2_val]));
        %         temp_dxdeta = eval(subs(dxdeta,[xi,eta,n,xi1,xi2],[xi_val(i,j),eta_val(i,j),k,xi1_val,xi2_val]));
        %         temp_dydeta = eval(subs(dydeta,[xi,eta,n,xi1,xi2],[xi_val(i,j),eta_val(i,j),k,xi1_val,xi2_val]));
        
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
        %         x_store(i,j) = eval(subs(x,[xi,eta],[xi_val(i,j),eta_val(i,j)]));
        %         y_store(i,j) = eval(subs(y,[xi,eta],[xi_val(i,j),eta_val(i,j)]));
        x_store(i,j) = a*sinh(xi_val(i,j))/(cosh(xi_val(i,j))-cos(eta_val(i,j)));
        y_store(i,j) = a*sin(eta_val(i,j))/(cosh(xi_val(i,j))-cos(eta_val(i,j)));
        
    end
end

%% plot stuff

figure(1)
plot(x_store,y_store,'kx')
axis square

figure(2)
quiver(x_store,y_store,u_store,v_store)
axis square

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

% % add in boundary values
% ni = 50;
% 
% x_min = meshdomain(1);
% x_max = meshdomain(2);
% y_min = meshdomain(3);
% y_max = meshdomain(4);
% 
% x1 = linspace(x_min,x_max,ni*2);
% x2 = x_max*ones(1,ni*2);
% x3 = linspace(x_min,x_max,ni*2);
% x4 = x_min*ones(1,ni*2);
% 
% y1 = y_min*ones(1,ni*2);
% y2 = linspace(y_min,y_max,ni*2);
% y3 = y_max*ones(1,ni*2);
% y4 = linspace(y_min,y_max,ni*2);
% 
% x = [x1, x2(2:end), x3(1:end-1), x4(2:end-1)];
% y = [y1, y2(2:end), y3(1:end-1), y4(2:end-1)];
% 
% x_all(last_count:last_count + length(x) - 1) = x;
% y_all(last_count:last_count + length(y) - 1) = y;
% u_all(last_count:last_count + length(x) - 1) = U;
% v_all(last_count:last_count + length(x) - 1) = 0;
% 
% u_interp = TriScatteredInterp(x_all',y_all',u_all');
% v_interp = TriScatteredInterp(x_all',y_all',v_all'); 
% 
% eval(['save flow_interp_U', num2str(U),'_d',num2str(100*d),' u_interp v_interp'])