function [ rho_store ] = getRho_alt( xe, ye, d, r_e, r_c,  ...
    mu, epsi, rho0, meshdomain, ind_fig,u_interp,v_interp,U,uxinterp,uyinterp)

% initialize storage variables
rho_store = 0;
last_count = 1;

x_storen = zeros(length(xe),5000);
y_storen = zeros(length(xe),5000);
rho_storen = zeros(length(xe),5000);

parfor j = 1:length(xe)
            
    % initialize spatial location
    x_temp = xe(j);
    y_temp = ye(j);
    
    % initialize time
    time = 0;
    
    % intialize distance from collector
    dist_c = sqrt((x_temp - 0)^2 + (y_temp + r_c)^2);
    dist_e = sqrt((x_temp - 0)^2 + (y_temp - (d + r_e))^2);
    
    % initialize index variable
    count = 1;
    
    % initialize storage variables
    x_store = zeros(1,5000);
    y_store = zeros(1,5000);
    t_store = zeros(1,5000);
    
    while dist_c > r_c
        
        % get electric field at current point
        %[Ex_temp, Ey_temp] = tri2gridCKG(p,t,temp,x_temp,y_temp);
        [Ex_temp, Ey_temp] = tri2gridDEX(x_temp,y_temp,uxinterp,uyinterp);
        
        if isnan(Ex_temp)
            break
        end
        
        % determine step size based on local E-field strength
        E_temp = sqrt(Ex_temp^2 + Ey_temp^2);
        
        if dist_e < 0.01 || dist_c < 0.01
            dt = 1e-5/(mu*E_temp);
        else
            dt = 1e-4/(mu*E_temp);
        end
        
        % determine change in characteristic line
        dxdt = mu*Ex_temp ;
        dydt = mu*Ey_temp ;
        
        x_store(count) = x_temp;
        y_store(count) = y_temp;
        
        x_temp = x_temp + dxdt*dt;
        y_temp = y_temp + dydt*dt;
        
        % determine change in time
        t_store(count) = time;
        
        time = time + dt;
        
        % increment index
        count = count + 1;
        
        % recompute distances from collector
        dist_c = sqrt((x_temp - 0)^2 + (y_temp + r_c)^2);
        dist_e = sqrt((x_temp - 0)^2 + (y_temp - (d + r_e))^2);
        
    end
    
    if ~isnan(Ex_temp) % only take final step to collector if still in solution domain
        
        x_temp = x_temp - dxdt*dt;
        y_temp = y_temp - dydt*dt;
        time = time - dt;
        
        % take final step to collector surface
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
    rho = (1/rho0 + (mu/epsi)*t_store).^(-1);
    
    x_storen(j,:) = x_store;
    y_storen(j,:) = y_store;
    rho_storen(j,:) = rho;
    
end

% put storage matrix together
for i = 1:length(xe)
    for j = 1:5000
        if x_storen(i,j)==0&&y_storen(i,j)==0            
            rho_store(1,last_count:last_count + j - 2) = x_storen(i,1:j-1);
            rho_store(2,last_count:last_count + j - 2) = y_storen(i,1:j-1);
            rho_store(3,last_count:last_count + j - 2) = rho_storen(i,1:j-1);
            last_count = last_count + j - 1;
            rho_store(1,last_count:last_count + j - 2) = -x_storen(i,1:j-1);
            rho_store(2,last_count:last_count + j - 2) = y_storen(i,1:j-1);
            rho_store(3,last_count:last_count + j - 2) = rho_storen(i,1:j-1);
            last_count = last_count + j - 1;
            break
        end
    end    
end

end













