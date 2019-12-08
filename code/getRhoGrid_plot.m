function [rho_store, oob_flag, pos_flag]  = getRhoGrid_plot( xi, yi, d, r_e, r_c,...
    mu, epsi, rho0, meshdomain, ind_fig,U,uxinterp,uyinterp)

global plot_debug do_flow

% initialize storage variables
rho_store = 0;
last_count = 1;
oob_flag = 0;
pos_flag = 0;

num_steps = 20000;
x_storen = zeros(length(xi),num_steps);
y_storen = zeros(length(xi),num_steps);
rho_storen = zeros(length(xi),num_steps);
count_storen = zeros(1,length(xi));

x_storen2 = zeros(length(xi),num_steps);
y_storen2 = zeros(length(xi),num_steps);
rho_storen2 = zeros(length(xi),num_steps);
count_storen2 = zeros(1,length(xi));

if length(xi) == 1
    endind = length(xi);
else
    endind = length(xi)-1;
end

for i = 1:endind+1
    %% go backwards to emitter
    % initialize spatial location
    x_temp = xi(i);
    y_temp = yi;
    
    % initialize time
    time = 0;
    
    % intialize distance from emitter
    dist_e = sqrt((x_temp - 0)^2 + (y_temp - (d + r_e))^2);
    
    % initialize index variable
    count = 1;
    
    % initialize storage variables
    x_store = zeros(1,num_steps);
    y_store = zeros(1,num_steps);
    t_store = zeros(1,num_steps);
    
    out_of_bounds = 0;
    
    while dist_e > r_e
        
        [Ex_temp, Ey_temp] = tri2gridDEX(x_temp,y_temp,uxinterp,uyinterp);
        
        % get potential flow contribution
        if do_flow
            [ux, uy] = getPF( x_temp, y_temp, r_c, U );

            if isnan(ux)
                ux = 0;
            end
            if isnan(uy)
                uy  = -U;
            end
        else
            ux = 0;
            uy = 0;
        end
        
        % determine step size based on local E-field strength
        E_temp = sqrt(Ex_temp^2 + Ey_temp^2);
        u_temp = sqrt(ux^2 + uy^2);
        if dist_e < 0.02
            dt = 1e-4/(mu*E_temp + u_temp);
        else
            dt = 1e-3/(mu*E_temp + u_temp);
        end
        
        % determine change in characteristic line
        dxdt = mu*Ex_temp + ux;
        dydt = mu*Ey_temp + uy;
        
        x_store(count) = x_temp;
        y_store(count) = y_temp;
        
        % step backwards
        x_temp = x_temp - dxdt*dt;
        y_temp = y_temp - dydt*dt;
        
        if x_temp > 0
            pos_flag = 1;
        end
        
        % determine change in time
        t_store(count) = time;
        
        time = time + dt;
        
        % increment index
        count = count + 1;
        
        % recompute distances from emitter
        dist_e = sqrt((x_temp - 0)^2 + (y_temp - (d + r_e))^2);
        
    end
    
    if isnan(Ex_temp) || isnan(Ey_temp)
        
        count = count - 1;
        rho = 0;
        out_of_bounds = 1;
        oob_flag = 1;
        
    else
        
        x_temp = x_temp + dxdt*dt;
        y_temp = y_temp + dydt*dt;
        time = time - dt;
        
        % take final step to emitter surface
        dt_temp = linspace(0,dt,100);
        x_dt = x_temp - dxdt*dt_temp;
        y_dt = y_temp - dydt*dt_temp;
        dist_temp = sqrt((x_dt - 0).^2 + (y_dt - (d + r_e)).^2);
        [~, ind] = min(abs(dist_temp - r_e));
        
        dt = dt_temp(ind);
        
        x_temp = x_temp - dxdt*dt;
        y_temp = y_temp - dydt*dt;
        time = time + dt;
        x_store(count) = x_temp;
        y_store(count) = y_temp;
        t_store(count) = time;
        
        t_store = abs(t_store - time);
        
        % solve for space charge along current characteristic line
        rho = (1/rho0 + (mu/epsi)*t_store).^(-1);
        
    end
    
    if plot_debug
        figure(ind_fig)
        hold on
        plot(x_store(1:count),y_store(1:count))
        plot(-x_store(1:count),y_store(1:count))
        axis(meshdomain)
    end
    
    x_storen(i,:) = x_store;
    y_storen(i,:) = y_store;
    rho_storen(i,:) = rho;
    count_storen(i) = count;
    
    %         % assign space charge to point and place in master storage variable
    %         rho_store(1,last_count:last_count + count - 1) = x_store(1:count);
    %         rho_store(2,last_count:last_count + count - 1) = y_store(1:count);
    %         rho_store(3,last_count:last_count + count - 1) = rho(1:count);
    %         last_count = last_count + count;
    %         rho_store(1,last_count:last_count + count - 1) = -x_store(1:count);
    %         rho_store(2,last_count:last_count + count - 1) = y_store(1:count);
    %         rho_store(3,last_count:last_count + count - 1) = rho(1:count);
    %         last_count = last_count + count;
    
    t_add = t_store(1);
    
    %% go forward to collector
    
    % re-initialize spatial location
    x_temp = xi(i);
    y_temp = yi;
    
    % initialize time
    time = t_add;
    
    % initialize index variable
    count = 1;
    
    % initialize storage variables
    x_store = zeros(1,num_steps);
    y_store = zeros(1,num_steps);
    t_store = zeros(1,num_steps);
    
    % intialize distance from collector
    dist_c = sqrt((x_temp - 0)^2 + (y_temp + r_c)^2);
    
    first_step = 1;
    
    while dist_c > r_c
        
        [Ex_temp, Ey_temp] = tri2gridDEX(x_temp,y_temp,uxinterp,uyinterp);
        
        % get potential flow contribution
        if do_flow
            [ux, uy] = getPF( x_temp, y_temp, r_c, U );
            
            if isnan(ux)
                ux = 0;
            end
            if isnan(uy)
                uy  = -U;
            end
        else
            ux = 0;
            uy = 0;
        end
        
        if isnan(Ex_temp)
            break
        end
        
        % determine step size based on local E-field strength
        E_temp = sqrt(Ex_temp^2 + Ey_temp^2);
        if dist_c < 0.02
            dt = 1e-4/(mu*E_temp + u_temp);
        else
            dt = 1e-3/(mu*E_temp + u_temp);
        end
        
        % determine change in characteristic line
        dxdt = mu*Ex_temp + ux;
        dydt = mu*Ey_temp + uy;
        
        if first_step == 0
            
            x_store(count) = x_temp;
            y_store(count) = y_temp;
            
        end
        
        % step forward
        x_temp = x_temp + dxdt*dt;
        y_temp = y_temp + dydt*dt;
        
        if x_temp > 0
            pos_flag = 1;
        end
        
        % determine change in time
        if first_step == 0
            
            t_store(count) = time;
            
            % increment index
            count = count + 1;
            
        end
        
        time = time + dt;
        
        % recompute distances from collector
        dist_c = sqrt((x_temp - 0)^2 + (y_temp + r_c)^2);
        
        if first_step
            first_step = 0;
        end
        
    end
    
    if ~isnan(Ex_temp) % only take final step to collector if still in solution domain
        
        x_temp = x_temp - dxdt*dt;
        y_temp = y_temp - dydt*dt;
        time = time - dt;
        
        % take final step to emitter surface
        dt_temp = linspace(0,dt,100);
        x_dt = x_temp + dxdt*dt_temp;
        y_dt = y_temp + dydt*dt_temp;
        dist_temp = sqrt((x_dt - 0).^2 + (y_dt + r_c).^2);
        [~, ind] = min(abs(dist_temp - r_c));
        dt = dt_temp(ind);
        
        x_temp = x_temp + dxdt*dt;
        y_temp = y_temp + dydt*dt;
        time = time + dt;
        x_store(count) = x_temp;
        y_store(count) = y_temp;
        t_store(count) = time;
        
    else
        
        count = count - 1;
        
    end
    
    % solve for space charge along current characteristic line
    if out_of_bounds
        rho = 0;
    else
        rho = (1/rho0 + (mu/epsi)*t_store).^(-1);
    end
    
    x_storen2(i,:) = x_store;
    y_storen2(i,:) = y_store;
    rho_storen2(i,:) = rho;
    count_storen2(i) = count;
    
    %         % assign space charge to point and place in master storage variable
    %         rho_store(1,last_count:last_count + count - 1) = x_store(1:count);
    %         rho_store(2,last_count:last_count + count - 1) = y_store(1:count);
    %         rho_store(3,last_count:last_count + count - 1) = rho(1:count);
    %         last_count = last_count + count;
    %         rho_store(1,last_count:last_count + count - 1) = -x_store(1:count);
    %         rho_store(2,last_count:last_count + count - 1) = y_store(1:count);
    %         rho_store(3,last_count:last_count + count - 1) = rho(1:count);
    %         last_count = last_count + count;
    
    if plot_debug
        figure(ind_fig)
        hold on
        plot([xi(i),x_store(1:count)],[yi,y_store(1:count)])
        plot([-xi(i),-x_store(1:count)],[yi,y_store(1:count)])
        axis(meshdomain)
    end
    
end

if plot_debug
    hold off
end

% put storage matrix together
for i = 1:length(xi)
    rho_store(1,last_count:last_count + count_storen(i)-1) = x_storen(i,1:count_storen(i));
    rho_store(2,last_count:last_count + count_storen(i)-1) = y_storen(i,1:count_storen(i));
    rho_store(3,last_count:last_count + count_storen(i)-1) = rho_storen(i,1:count_storen(i));
    last_count = last_count + count_storen(i);
    rho_store(1,last_count:last_count + count_storen(i)-1) = -x_storen(i,1:count_storen(i));
    rho_store(2,last_count:last_count + count_storen(i)-1) = y_storen(i,1:count_storen(i));
    rho_store(3,last_count:last_count + count_storen(i)-1) = rho_storen(i,1:count_storen(i));
    last_count = last_count + count_storen(i);
    rho_store(1,last_count:last_count + count_storen2(i)-1) = x_storen2(i,1:count_storen2(i));
    rho_store(2,last_count:last_count + count_storen2(i)-1) = y_storen2(i,1:count_storen2(i));
    rho_store(3,last_count:last_count + count_storen2(i)-1) = rho_storen2(i,1:count_storen2(i));
    last_count = last_count + count_storen2(i);
    rho_store(1,last_count:last_count + count_storen2(i)-1) = -x_storen2(i,1:count_storen2(i));
    rho_store(2,last_count:last_count + count_storen2(i)-1) = y_storen2(i,1:count_storen2(i));
    rho_store(3,last_count:last_count + count_storen2(i)-1) = rho_storen2(i,1:count_storen2(i));
    last_count = last_count + count_storen2(i);
end

%% add in zero on boundaries

%outer domain
ni = length(xi);

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

rho_store(1,last_count:last_count + length(x) - 1) = x;
rho_store(2,last_count:last_count + length(y) - 1) = y;
rho_store(3,last_count:last_count + length(x) - 1) = 0;

last_count = last_count + length(x);

% % inner domain
% x_min = meshdomain_in(1);
% x_max = meshdomain_in(2);
% y_min = meshdomain_in(3);
% y_max = meshdomain_in(4);
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
% rho_store(1,last_count:last_count + length(x) - 1) = x;
% rho_store(2,last_count:last_count + length(y) - 1) = y;
% rho_store(3,last_count:last_count + length(x) - 1) = 0;

end







