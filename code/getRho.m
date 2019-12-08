function [ rho_store, rhoc_store ] = getRho( xe, ye, d, r_e, r_c,  ...
    mu, epsi, rho0, meshdomain,uxinterp,uyinterp)

% initialize storage variables
rho_store = 0;
last_count = 1;

num_steps = 20000;
x_storen = zeros(length(xe),num_steps);
y_storen = zeros(length(xe),num_steps);
rho_storen = zeros(length(xe),num_steps);
count_storen = zeros(1,length(xe));
no_nan = ones(1,length(xe));

for i = 1:length(xe)
    
    % initialize spatial location
    x_temp = xe(i);
    y_temp = ye(i);
    
    % initialize time
    time = 0;
    
    % intialize distance from collector
    dist_c = sqrt((x_temp - 0)^2 + (y_temp + r_c)^2);
    dist_e = sqrt((x_temp - 0)^2 + (y_temp - (d + r_e))^2);
    
    % initialize index variable
    count = 1;
    
    % initialize storage variables
    x_store = zeros(1,num_steps);
    y_store = zeros(1,num_steps);
    t_store = zeros(1,num_steps);
    
    % generate function for RKF
    
   
        prime = @(t,x)mu.*Eeval(t,x,uxinterp,uyinterp);
    
%     htest = figure;
    
    while dist_c > r_c
        
        % get electric field at current point
        [Ex_temp, Ey_temp] = tri2gridDEX(x_temp,y_temp,uxinterp,uyinterp);
        
        
            ux = 0;
            uy = 0;
       
        
        if isnan(Ex_temp)
            no_nan(i) = 0;
            break
        end
        
        % store variables of interest
        x_store(count) = x_temp;
        y_store(count) = y_temp;
        t_store(count) = time;
        
        % determine step size
        E_temp = sqrt(Ex_temp^2 + Ey_temp^2);
        u_temp = sqrt(ux^2 + uy^2);
        if dist_e < 0.03 || dist_c < 0.03
            dt = 1e-4/(mu*E_temp + u_temp);
        else
            dt = 1e-3/(mu*E_temp + u_temp);
        end
        
        if time == 0
            h = dt;
        end
        
        [x, ~, h, h_old, dxdt] = RK23( [x_temp; y_temp],time,h,prime,1e-8 );
        x_temp = x(1,:);
        y_temp = x(2,:);
        dydt = dxdt(2,:);
        dxdt = dxdt(1,:);
        dt = h_old;
        time = time+h_old;
        
%         figure(htest)
%         line(x_temp,y_temp,'Marker','o','LineStyle','none')
        
        if isnan(x_temp)
            no_nan(i) = 0;
            break
        end
        
        % increment index
        count = count + 1;
        
        % recompute distances from collector
        dist_c = sqrt((x_temp - 0)^2 + (y_temp + r_c)^2);
        dist_e = sqrt((x_temp - 0)^2 + (y_temp - (d + r_e))^2);
        
        %break if outside domain
        if ((x_temp<meshdomain(1)) || (x_temp > meshdomain(2)) || (y_temp < meshdomain(3)) || y_temp > meshdomain(4))
            %no_nan(i) = 0;
            break
        end
        
        %[dist_c,dist_e]
        %figure(1)
        %hold on
        %plot(x_temp,y_temp,'.')
        %drawnow
        
    end
    
    if ~isnan(Ex_temp) || ~isnan(x_temp) % only take final step to collector if still in solution domain
        
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
                RK4( x_temp,y_temp,Ex_temp,Ey_temp,time,dt_temp(mmm),mu,uxinterp,uyinterp,r_c,1 );
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
        
        %rhoc_store(:,i) = [x_store(count);y_store(count)];
        
    else
        count = count - 1;
    end
    
    % solve for space charge along current characteristic line
    rho = (1/rho0 + (mu/epsi)*t_store).^(-1);
    
    
    if ~isnan(Ex_temp) || ~isnan(x_temp)
        x_storen(i,:) = x_store;
        y_storen(i,:) = y_store;
        rho_storen(i,:) = rho;
        count_storen(i) = count;
    end
    
    if 0 %useful for extensive debugging
        figure;
        subplot(312)
        line(x_store(1:count),y_store(1:count))
        axis(meshdomain)
        drawnow
    end
    
end

for i = 1:length(xe)
    if no_nan(i) == 1
        rho_store(1,last_count:last_count + count_storen(i)-1) = x_storen(i,1:count_storen(i));
        rho_store(2,last_count:last_count + count_storen(i)-1) = y_storen(i,1:count_storen(i));
        rho_store(3,last_count:last_count + count_storen(i)-1) = rho_storen(i,1:count_storen(i));
        last_count = last_count + count_storen(i);
        rho_store(1,last_count:last_count + count_storen(i)-1) = -x_storen(i,1:count_storen(i));
        rho_store(2,last_count:last_count + count_storen(i)-1) = y_storen(i,1:count_storen(i));
        rho_store(3,last_count:last_count + count_storen(i)-1) = rho_storen(i,1:count_storen(i));
        last_count = last_count + count_storen(i);
    end
end

% remove (0,0) from rhoc_store
%keep = ones(1,length(rhoc_store));
%for i = 1:length(rhoc_store)
%    if rhoc_store(1,i) == 0 && rhoc_store(2,i) == 0
%        keep(i) = 0;
%    end
%end

rhoc_store = rho_store;%(:,keep==1);

end













