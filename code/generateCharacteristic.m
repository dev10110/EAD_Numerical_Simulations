function [rho_store, endpos] = generateCharacteristic(theta, d, r_e, r_c, mu, epsi, rho0, meshdomain,uxinterp,uyinterp)

num_steps = 20000;

%find starting point
x_start = r_e*cos(theta);
y_start = (d + r_e) + r_e*sin(theta);


x_temp = x_start;
y_temp = y_start;

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


while (dist_c > r_c) && (count< num_steps)
    
    % get electric field at current point
    [Ex_temp, Ey_temp] = tri2gridDEX(x_temp,y_temp,uxinterp,uyinterp);
    
    if isnan(Ex_temp)
        no_nan(i) = 0;
        break
    end
    
    % store variables of interest
    x_store(count) = x_temp;
    y_store(count) = y_temp;
    t_store(count) = time;
    
    
if 0 %useful for extensive debugging

    plot(x_store(end),y_store(end),'.')
    axis(meshdomain)
    drawnow
    axis equal
end
    
    % determine step size
    E_temp = norm(Ex_temp,Ey_temp);
    
    
    if dist_e < 0.03 || dist_c < 0.03
        dt = 1e-3/(mu*E_temp );
    else
        dt = 1e-2/(mu*E_temp);
    end
    
    if time == 0
        h = dt;
    end
    
    [x, ~, h, h_old, dxdt] = RK23( [x_temp; y_temp],time,h,prime,1e-7 );
    x_temp = x(1,:);
    y_temp = x(2,:);
    dydt = dxdt(2,:);
    dxdt = dxdt(1,:);
    dt = h_old;
    time = time+h_old;
    
    if isnan(x_temp)
        no_nan(i) = 0;
        break
    end
    
    
    % increment index
    count = count + 1;
    
    % recompute distances from collector
    dist_c = sqrt((x_temp - 0)^2 + (y_temp + r_c)^2);
    dist_e = sqrt((x_temp - 0)^2 + (y_temp - (d + r_e))^2);
    
    %if dist_e < r_e
    %    rho_store = [];
    %    endpos = 0;
    %    return
    %end
    
    if dist_c <=  r_c 
        endpos = 4;
    end
    
    %break if outside domain
    if (x_temp<meshdomain(1)) %gone to the left
        endpos = 2;
        break;
    elseif (y_temp < meshdomain(3))%gone to the bottom
        endpos = 3;
        break
    elseif (y_temp > meshdomain(4)) %gone to the top
        endpos = 1;
        break
    end
    
end
if 0
    plot(x_store,y_store,'.')
    axis(meshdomain)
    drawnow
    axis equal
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

rho_store(1,:) = x_store;
rho_store(2,:) = y_store;
rho_store(3,:) = rho;



if 0 %useful for extensive debugging
    figure;
    subplot(312)
    line(x_store(1:count),y_store(1:count))
    axis(meshdomain)
    drawnow
end

end

