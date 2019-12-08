function [rho_store, oob_flag, inf_flag]  = getRhoGrid( xi, yi, d, r_e, r_c,...
    mu, epsi, rho0, meshdomain,ind_fig,U,uxinterp,uyinterp)

global plot_debug do_flow

% initialize storage variables
rho_store = 0;
last_count = 1;
oob_flag = 0;
inf_flag = 0;

num_steps = 20000;
x_storen = zeros(length(yi),num_steps);
y_storen = zeros(length(yi),num_steps);
rho_storen = zeros(length(yi),num_steps);
count_storen = zeros(1,length(yi));

x_storen2 = zeros(length(yi),num_steps);
y_storen2 = zeros(length(yi),num_steps);
rho_storen2 = zeros(length(yi),num_steps);
count_storen2 = zeros(1,length(yi));

if length(xi) == 1 && length(yi) == 1
    endind = 1;
elseif length(xi) > 1 && length(yi) == 1
    endind = length(xi)-1;
elseif  length(xi) == 1 && length(yi) > 1
    endind = length(yi)-1;
end

for i = 1:endind
    %% go backwards to emitter
    % initialize spatial location
    
    if length(xi) == 1 && length(yi) == 1
        x_temp = xi;
        y_temp = yi;
    elseif length(xi) > 1 && length(yi) == 1
        x_temp = xi(i);
        y_temp = yi;
    elseif  length(xi) == 1 && length(yi) > 1
        x_temp = xi;
        y_temp = yi(i);
    end
    
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
        
        nan_switch = 0;
        
        if isnan(Ex_temp)
            nan_switch = 1;
            break
        end
        
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
        if ~do_flow
            if dist_e < 0.005
                dt = 1e-4/(mu*E_temp + u_temp);
            elseif dist_e < 0.03
                dt = 1e-3/(mu*E_temp + u_temp);
            else
                dt = 1e-2/(mu*E_temp + u_temp);
            end
        else
            if dist_e < 0.005
                dt = 1e-4/(mu*E_temp + u_temp);
            elseif dist_e < 0.1
                dt = 1e-3/(mu*E_temp + u_temp);
            else
                dt = 1e-2/(mu*E_temp + u_temp);
            end
        end

        x_store(count) = x_temp;
        y_store(count) = y_temp;
        t_store(count) = time;
     
        
        %[ x_temp,y_temp,time,dxdt,dydt ] = Eulerf( x_temp,y_temp,-Ex_temp,-Ey_temp,time,dt,mu,-ux,-uy );
        [ x_temp,y_temp,time,dxdt,dydt,x_old,y_old ] = ...
            RK4( x_temp,y_temp,Ex_temp,Ey_temp,ux,uy,time,dt,mu,uxinterp,uyinterp,r_c,U,0 );
        
        % increment index
        count = count + 1;
        
        if isnan(x_temp)
            nan_switch = 1;
            x_temp = x_old;
            y_temp = y_old;
            break
        end
        
        % recompute distances from emitter
        dist_e = sqrt((x_temp - 0)^2 + (y_temp - (d + r_e))^2);
        
              
        if count > 20000
            break
        end
        
%         figure(1)
%         hold on
%         plot(x_temp,y_temp,'.')
%         drawnow
        
    end    
   
    if nan_switch
        
        count = count - 1;
        rho = 0;
        out_of_bounds = 1;
        oob_flag = 1;
        
    elseif count > 20000
        
        count = count - 1;
        inf_flag = 1;
        rho = 0;
        
    else
        
        x_temp = x_temp - dxdt*dt;
        y_temp = y_temp - dydt*dt;
        time = time - dt;
        
        % take final step to emitter surface
        dt_temp = linspace(0,dt,25);
        x_last = zeros(1,length(dt_temp));
        y_last = zeros(1,length(dt_temp));
        time_last = zeros(1,length(dt_temp));
        for mmm = 1:length(dt_temp)
            [ x_last(mmm),y_last(mmm),time_last(mmm) ] = ...
                RK4( x_temp,y_temp,Ex_temp,Ey_temp,ux,uy,time,dt_temp(mmm),mu,uxinterp,uyinterp,r_c,U,0 );
        end

        dist_temp = sqrt((x_last - 0).^2 + (y_last - (d + r_e)).^2);
        [~, ind] = min(abs(dist_temp - r_e));
        dt = dt_temp(ind);
        
        x_temp = x_last(ind);
        y_temp = y_last(ind);
        time = time_last(ind);
        
        x_store(count) = x_temp;
        y_store(count) = y_temp;
        t_store(count) = time;
        
        t_store = abs(t_store - time);
        
        % solve for space charge along current characteristic line
        rho = (1/rho0 + (mu/epsi)*t_store).^(-1);
        
    end
    
    if plot_debug
        figure(1)
        hold on;
        plot(x_store(1:count),y_store(1:count),'r')
        axis(meshdomain)
        drawnow
    end
    
    x_storen(i,:) = x_store;
    y_storen(i,:) = y_store;
    rho_storen(i,:) = rho;
    count_storen(i) = count;
        
    t_add = t_store(1);
    
    if plot_debug
        figure(2)
        hold on
        scatter(x_store,y_store,[],rho)
    end
    %% go forward to collector
    
    % re-initialize spatial location
    if length(xi) == 1 && length(yi) == 1
        x_temp = xi;
        y_temp = yi;
    elseif length(xi) > 1 && length(yi) == 1
        x_temp = xi(i);
        y_temp = yi;
    elseif  length(xi) == 1 && length(yi) > 1
        x_temp = xi;
        y_temp = yi(i);
    end
    
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
        
        nan_switch = 0;
        
        if isnan(Ex_temp)
            nan_switch = 1;
            break
        end
        
        % determine step size based on local E-field strength
        E_temp = sqrt(Ex_temp^2 + Ey_temp^2);
        if ~do_flow
            if dist_c < 0.03
                dt = 1e-3/(mu*E_temp + u_temp);
            else
                dt = 1e-2/(mu*E_temp + u_temp);
            end
        else
            if dist_c < 0.1
                dt = 1e-3/(mu*E_temp + u_temp);
            else
                dt = 1e-2/(mu*E_temp + u_temp);
            end
        end
        
        if first_step == 0            
            x_store(count) = x_temp;
            y_store(count) = y_temp;
            t_store(count) = time;
            
            % increment index
            count = count + 1;            
        end
        
        %[ x_temp,y_temp,time,dxdt,dydt ] = Eulerf( x_temp,y_temp,Ex_temp,Ey_temp,time,dt,mu,ux,uy );
        [ x_temp,y_temp,time,dxdt,dydt,x_old,y_old ] = ...
            RK4( x_temp,y_temp,Ex_temp,Ey_temp,ux,uy,time,dt,mu,uxinterp,uyinterp,r_c,U,1 );
        
        if isnan(x_temp)
            nan_switch = 1;
            x_temp = x_old;
            y_temp = y_old;
            break
        end
        
        % recompute distances from collector
        dist_c = sqrt((x_temp - 0)^2 + (y_temp + r_c)^2);
        
        if first_step
            first_step = 0;
        end
        
        %break if outside domain
        if ((x_temp<meshdomain(1)) || (x_temp > meshdomain(2)) || (y_temp < meshdomain(3)) || y_temp > meshdomain(4))
            %nan_switch = 1;
            break
        end
        
        %figure(1)
        %hold on
        %plot(x_temp,y_temp,'.')
        %drawnow
        
    end
    
    if ~nan_switch  % only take final step to collector if still in solution domain
        
        x_temp = x_temp - dxdt*dt;
        y_temp = y_temp - dydt*dt;
        time = time - dt;
        
        % take final step to collector surface
        dt_temp = linspace(0,dt,25);
        x_last = zeros(1,length(dt_temp));
        y_last = zeros(1,length(dt_temp));
        time_last = zeros(1,length(dt_temp));
        for mmm = 1:length(dt_temp)
            [ x_last(mmm),y_last(mmm),time_last(mmm) ] = ...
                RK4( x_temp,y_temp,Ex_temp,Ey_temp,ux,uy,time,dt_temp(mmm),mu,uxinterp,uyinterp,r_c,U,1 );
        end

        dist_temp = sqrt((x_last - 0).^2 + (y_last + r_c).^2);
        [~, ind] = min(abs(dist_temp - r_c));
        dt = dt_temp(ind);

        x_temp = x_last(ind);
        y_temp = y_last(ind);
        time = time_last(ind);
        
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
    
    if plot_debug
        figure(ind_fig)
        line(x_store(1:count),y_store(1:count))
        axis(meshdomain)
    end
    
    figure(2)
    hold on
    scatter(x_store,y_store,[],rho)
    
end

% if plot_debug
%     hold off
% end

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
if length(xi)>1
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
end

end







