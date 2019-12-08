function diff = find_emitter_eField(rho0, Ecrit, p,e,t, geometry, constants, plot_debug)

%disp(rho0)


r_c = geometry.r_c;
r_e = geometry.r_e;
d = geometry.d;
meshdomain = geometry.meshdomain;


mu = constants.mu;
epsi = constants.epsi;
num_emit_default = 40;


%get the emitter and collector points
[pts_e, pts_c, e_inds, c_inds] = getBPts(p,e,t,d,r_e,r_c);
[xc,yc,dAc,thetac,c_inds] = getMOCEnd( num_emit_default,r_c,pts_c,c_inds);
[xe_alt,ye_alt,dAe,theta,e_inds] = getMOCStart_alt( num_emit_default, d, r_e,pts_e,e_inds); % for current calculation




%% solve for the phi

b = @pdebound_1d;
a = 0;
c = 1;

f = 0; %assume 0 space charge at first

[phi, ~] = pdenonlin(b,p,e,t,c,a,f);


% find electric field
[Ex, Ey] = pdegrad(p,t,phi);
grad_phi = [Ex; Ey];

% initialize potential field storage variables
phi_last = phi;



%start it high
exit_test = 1;

while exit_test > 0.005
    
    if plot_debug
        figure;
        subplot(221)
        pdeplot(p,e,t,'XYData',phi_last);
        %scatter(p(1,:),p(2,:),[],phi_last,'.');
        colorbar
        title('Potential')
        axis equal
        drawnow
    end
    
    % find interpolation surface for electric field
    [ uxinterp,uyinterp ] = getEInterp( p,t,grad_phi,r_c );
    
    if plot_debug
        subplot(222)
        quiver(uxinterp.Points(:,1),uxinterp.Points(:,2),(uxinterp.Values)./sqrt(uxinterp.Values.^2+uyinterp.Values.^2),uyinterp.Values./sqrt(uxinterp.Values.^2+uyinterp.Values.^2),'AutoScaleFactor',1,'MaxHeadSize',0.01)
        axis equal
        title('Unit E Field Vector')
        drawnow
    end
    
    % find all starting points from emitter
    num_emit = num_emit_default;
    [xe,ye,theta_full] = getMOCStartDev( num_emit, d, r_e, []);
    
    %get space charge dist
    [rho_store, ~] = ...
        getRho(xe, ye, d, r_e, r_c, mu, ...
        epsi, rho0, meshdomain,uxinterp,uyinterp);
    
    touches_all_bounds = 0;
    
    if plot_debug
            subplot(223)
            scatter(rho_store(1,:),rho_store(2,:),[],rho_store(3,:),'.');
            colorbar
            title('Space Charge')
            axis equal
            drawnow
    end
    
    while (touches_all_bounds == 0)
        
        %check if the space charge lines touch all bounds:
        %check that the min to max of x is within 99% of the full domain:
        if ((max(rho_store(1,:))-min(rho_store(1,:))) > 0.99*(meshdomain(2)-meshdomain(1))) && ((max(rho_store(2,:))-min(rho_store(2,:))) > 0.99*(meshdomain(4)-meshdomain(3)))
            touches_all_bounds = 1;
        else
            %if it isnt, double the number of exit points  and check
            %again
            num_emit = round(num_emit*2);
            fprintf(', %i', num_emit);
            [xe_new,ye_new,theta_full] = getMOCStartDev( num_emit, d, r_e, theta_full);
            
            
            [rho_store_temp, ~] = ...
                getRho(xe_new, ye_new, d, r_e, r_c, mu, ...
                epsi, rho0, meshdomain,uxinterp,uyinterp);
            
            rho_store_temp2 = [rho_store, rho_store_temp];
            rho_store = rho_store_temp2;
            %rhoc_store = [rhoc_store, rhoc_store_temp];
            
        end
        
        if plot_debug
            subplot(223)
            scatter(rho_store(1,:),rho_store(2,:),[],rho_store(3,:),'.');
            colorbar
            title('Space Charge')
            axis equal
            drawnow
        end
        
    end
    
%     if touches_all_bounds == 0 %ie. it still didnt touch all the boundaries, set the boundary charge to 0
%         %bottom edge
%         startx = meshdomain(1):0.01:meshdomain(2);
%         starty = meshdomain(3)*ones(size(startx));
%         rho_store_add = [startx;starty;zeros(size(startx))];
%         rho_store = [rho_store, rho_store_add];
%         
%         %top edge:
%         startx = meshdomain(1):0.01:meshdomain(2);
%         starty = meshdomain(4)*ones(size(startx));
%         rho_store_add = [startx;starty;zeros(size(startx))];
%         rho_store = [rho_store, rho_store_add];
%         
%         %left edge
%         starty = meshdomain(3):0.01:meshdomain(4);
%         startx = meshdomain(1)*ones(size(starty));
%         rho_store_add = [startx;starty;zeros(size(startx))];
%         rho_store = [rho_store, rho_store_add];
%         
%         %right edge
%         starty = meshdomain(3):0.01:meshdomain(4);
%         startx = meshdomain(2)*ones(size(starty));
%         rho_store_add = [startx;starty;zeros(size(startx))];
%         rho_store = [rho_store, rho_store_add];
%     end
%     
    
    
    %rho_store = getRho_efficient(d, r_e, r_c, mu, epsi, rho0, meshdomain,uxinterp,uyinterp);
    
    if plot_debug
        subplot(223)
        scatter(rho_store(1,:),rho_store(2,:),[],rho_store(3,:),'.');
        colorbar
        title('Space Charge')
        axis equal
        drawnow
    end
    % interpolate and fill FEM domain with space charge
    [f_new,finterp] = interpf( p,t,rho_store,epsi);
    
    if plot_debug
        %figure(ind_fig)
        subplot(224)
        % visualize space charge
        pdeplot(p,e,t,'xydata',finterp(p(1,:),p(2,:)),'colormap','jet')
        axis equal
        title('Space Charge Interp')
        drawnow
    end
    
    %% get the E field
    
    % emitter
    [ Exe,Eye ] = getE_alt( e_inds, grad_phi );
    
    % compute emitter e efield
    
    Ee = sum(sqrt(Exe.^2 + Eye.^2).*dAe)/sum(dAe); % 'average' electric field
    
    
    %% re-solve for potential field: Davis and Hoburg, Step 4
    
    [phi_new, ~] = pdenonlin(b,p,e,t,c,a,f_new);
    
    %% determine exit criteria, exit if met: Davis and Hoburg, Step 5
    
    exit_test = norm(phi_new-phi_last)/norm(phi_last);
    
    fprintf('\n rho0: %.3g, Ee: %.3g, exit: %.1f perc',rho0, Ee, exit_test*100);
    
    %% reset necessary variables
    omega = .5;
    phi_last = phi_new;
    phi = phi_new*omega + phi*(1-omega);
    
    [Ex, Ey] = pdegrad(p,t,phi);
    grad_phi = [Ex; Ey];
    
end


diff = Ee - Ecrit;
disp('diff')
disp(diff)
disp(' ')
