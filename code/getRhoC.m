function [ rho3_store ] = getRhoC( d, r_e, r_c,  ...
    mu, epsi, rho0, meshdomain, ind_fig,U,rhoc_store,uxinterp,uyinterp)

global plot_debug do_flow

%% find 'last' point at which we have charge data for on collector
rho_temp = rhoc_store(1,:);
rho_temp = rho_temp(rhoc_store(2,:) < -r_c);
rho_temp2 = rhoc_store(2,:);
rho_temp2 = rho_temp2(rhoc_store(2,:) < -r_c);
[x_max, max_ind] = max(rho_temp);
y_max = rho_temp2(max_ind);

theta = atan((y_max + r_c)/x_max);

% determine points I want to loop over
theta_pts = linspace(theta,pi()/2,10);

% initialize storage variables
rho3_store = 0;
last_count = 1;
oob = 0;

for j = 2:length(theta_pts)-1
            
    %% go backwards to emitter
        % initialize spatial location
        x_temp = -r_c*cos(theta_pts(j));
        y_temp = -r_c*sin(theta_pts(j)) - r_c;
        
        % initialize time
        time = 0;
        
        % intialize distance from emitter
        dist_e = sqrt((x_temp - 0)^2 + (y_temp - (d + r_e))^2);
        dist_c = sqrt((x_temp - 0)^2 + (y_temp  + r_c)^2);
        
        % initialize index variable
        count = 1;
        
        % initialize storage variables
        x_store = 0;
        y_store = 0;
        t_store = 0;
        
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
            if dist_e < 0.005
                dt = 1e-4/(mu*E_temp + u_temp);
            elseif dist_e < 0.03 || dist_c < 0.03
                dt = 1e-3/(mu*E_temp + u_temp);
            else
                dt = 1e-2/(mu*E_temp + u_temp);
            end            
           
            x_store(count) = x_temp;
            y_store(count) = y_temp;          
            t_store(count) = time;
            
            %[ x_temp,y_temp,time,dxdt,dydt ] = EulerfBack( x_temp,y_temp,Ex_temp,Ey_temp,time,dt,mu,ux,uy );
            [ x_temp,y_temp,time,dxdt,dydt ] = ...
                RK4( x_temp,y_temp,Ex_temp,Ey_temp,ux,uy,time,dt,mu,uxinterp,uyinterp,r_c,U,0 );
            
            % increment index
            count = count + 1;
            
            % recompute distances from emitter
            dist_e = sqrt((x_temp - 0)^2 + (y_temp - (d + r_e))^2);
            
        end
        
        if isnan(Ex_temp) || isnan(x_temp)
            
            count = count - 1;
            rho = 0;
            oob = 1;
            
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
            figure(ind_fig)
            hold on
            plot(x_store,y_store)
            axis(meshdomain)
        end
        
        %% store values
        
        rho3_store(1,last_count:last_count + count - 1) = x_store;
        rho3_store(2,last_count:last_count + count - 1) = y_store;
        rho3_store(3,last_count:last_count + count - 1) = rho;
        last_count = last_count + count;
        rho3_store(1,last_count:last_count + count - 1) = -x_store;
        rho3_store(2,last_count:last_count + count - 1) = y_store;
        rho3_store(3,last_count:last_count + count - 1) = rho;
        last_count = last_count + count;
        
end

%% routine for no out-of-bounds characteristic

if oob == 0
    'doing oob routine'
    if U < 10
        pert = 1e-3;
    else
        pert = 5e-4;
    end
    % find outermost point
    [x_max, ind] = max(abs(rho3_store(1,:)));
    y_max = rho3_store(2,ind);
    x = -x_max - pert;
    y = y_max;
    while oob == 0
    %while pos_flag
        % use getrhogrid to find new characteristic line
        [rho_vec, oob, inf_flag]  = getRhoGridOOB( x,y,d,r_e,r_c,...
            mu,epsi,rho0,meshdomain,ind_fig,U,uxinterp,uyinterp);
        x = x - pert;
        if ~inf_flag 
            rho3_store = [rho3_store,rho_vec];
        end
    end
    %rho3_store = [rho3_store,rho_vec];
end

hold off

end













