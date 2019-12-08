function [eta_1d, eta_alt, eta, Ee, Ie, Ic, Jd, Je, Jc, res_norm, grad_phi, phi_last, f_new, phi_lap] = ...
    EHDsolverf(p,e,t,r_e,r_c,d,kV,U,meshdomain,rho0,max_iter,first_it,phi,Sext,pbdry)

%% declare global variables

global do_flow plot_debug no_flux

%% set other parameters
mu = 2*10^-4;
epsi = 8.85418782e-12;
ng = 50;

%% initialize solution (i.e. Laplace): Davis and Hoburg, Step 1

% get boundary conditions through custom function
b = @pdebound_1d;
a = 0;
c = 1;

if first_it
    f = 0;
    % solve for potential field
    if no_flux
        [phi, ~] = pdenonlin(b,p,e,t,c,a,f);
    else
        [K,M,F,Q,G,H,R]=assempde(b,p,e,t,c,a,f);
        K(pbdry,pbdry) = K(pbdry,pbdry)+Sext(1:length(pbdry),1:length(pbdry));
        phi = assempde(K,M,F,Q,G,H,R);
    end
    
    phi_lap = phi;
    
end

% find electric field
[Ex, Ey] = pdegrad(p,t,phi);
grad_phi = [Ex; Ey];

% find all starting points from emitter
num_emit_default = 20;
num_emit = num_emit_default;
%[xe,ye] = getMOCStart( num_emit, d, r_e);
[pts_e, pts_c, e_inds, c_inds] = getBPts(p,e,t,d,r_e,r_c);
[xc,yc,dAc,thetac,c_inds] = getMOCEnd( num_emit,r_c,pts_c,c_inds);
[xe_alt,ye_alt,dAe,theta,e_inds] = getMOCStart_alt( num_emit, d, r_e,pts_e,e_inds); % for current calculation


[xe,ye,theta_full] = getMOCStartDev( num_emit, d, r_e, []);

xe = [0, xe];
ye = [(d + r_e), ye];

% get area of each FEM element
dAFEM = getdAFEM( p,t );

% initialize potential field storage variables
phi_last = phi;

dummy_ind_fig = 0;


%% begin iteration

% initialize iteration variables
i = 1;
exit_test = 1;


while exit_test > 0.01
    
    % show iteration number
    i
    
    figure;
    subplot(311)
    scatter(p(1,:),p(2,:),[],phi_last,'.');
    colorbar
    title('Potential')
    axis equal
    drawnow
    
    %% determine space charge distribution using MOC: Davis and Hoburg, Step 3
    
    % find interpolation surface for electric field
    [ uxinterp,uyinterp ] = getEInterp( p,t,grad_phi,r_c );
    [ uxinterpf,uyinterpf ] = getEInterpFull( p,t,grad_phi );
     
    %to plot the normalized e field strength
    
    %quiver(uxinterp.Points(:,1),uxinterp.Points(:,2),(uxinterp.Values)./sqrt(uxinterp.Values.^2+uyinterp.Values.^2),uyinterp.Values./sqrt(uxinterp.Values.^2+uyinterp.Values.^2),'AutoScaleFactor',1,'MaxHeadSize',0.01)
    touches_all_bounds = 0;
    % determine characteristic lines per starting point and return space
    % charge along lines
    %do this once
    [rho_store, rhoc_store] = ...
        getRho(xe, ye, d, r_e, r_c, mu, ...
        epsi, rho0, meshdomain, dummy_ind_fig,U,uxinterp,uyinterp);
     
    while (touches_all_bounds == 0) && (num_emit < 161)

    %check if the space charge lines touch all bounds:
        %check that the min to max of x is within 99% of the full domain:
        if ((max(rho_store(1,:))-min(rho_store(1,:))) > 0.99*(meshdomain(2)-meshdomain(1))) && ((max(rho_store(2,:))-min(rho_store(2,:))) > 0.99*(meshdomain(4)-meshdomain(3)))
            touches_all_bounds = 1;
        else
            %if it isnt, double the number of exit points  and check
            %again
            num_emit = round(num_emit*2)
            
            [xe_new,ye_new,theta_full] = getMOCStartDev( num_emit, d, r_e, theta_full);
            
            
            [rho_store_temp, rhoc_store_temp] = ...
                getRho(xe_new, ye_new, d, r_e, r_c, mu, ...
                epsi, rho0, meshdomain, dummy_ind_fig,U,uxinterp,uyinterp);
            
            rho_store_temp2 = [rho_store rho_store_temp];
            rho_store = rho_store_temp2;
            %rhoc_store = [rhoc_store, rhoc_store_temp];
            
        end
    
    
    subplot(312)
    scatter(rho_store(1,:),rho_store(2,:),[],rho_store(3,:),'.');
    colorbar
    title('Space Charge')
    axis equal
    drawnow
            
    end
    
    if touches_all_bounds == 0 %ie. it still didnt touch all the boundaries, set the boundary charge to 0
        %bottom edge
        startx = meshdomain(1):0.01:meshdomain(2);
        starty = meshdomain(3)*ones(size(startx));
        rho_store_add = [startx;starty;zeros(size(startx))];
        rho_store = [rho_store, rho_store_add];
        
        %top edge:
        startx = meshdomain(1):0.01:meshdomain(2);
        starty = meshdomain(4)*ones(size(startx));
        rho_store_add = [startx;starty;zeros(size(startx))];
        rho_store = [rho_store, rho_store_add];
        
        %left edge
        starty = meshdomain(3):0.01:meshdomain(4);
        startx = meshdomain(1)*ones(size(starty));
        rho_store_add = [startx;starty;zeros(size(startx))];
        rho_store = [rho_store, rho_store_add];
        
        %right edge
        starty = meshdomain(3):0.01:meshdomain(4);
        startx = meshdomain(2)*ones(size(starty));
        rho_store_add = [startx;starty;zeros(size(startx))];
        rho_store = [rho_store, rho_store_add];
    end
    
    %reset the number of emitter lines:
    num_emit = num_emit_default;
    [xe,ye,theta_full] = getMOCStartDev( num_emit, d, r_e, []);

    
    % calculate information required for getRhoGrid
    x_max = max(abs(rho_store(1,:)));
    ni = 6; % keep even to avoid infinite field lines
    %[xi,yi] = findInterps(ni,d,meshdomain,x_max);
    [xi,yi] = findInterpsy(ni,d,meshdomain,x_max);
    
    if 1 %prevent the get rho grid script from running
        rho2_store = [];
    else
    rho2_store = ...
        getRhoGrid(xi, yi,d, r_e, r_c, mu, epsi, rho0, ...
        meshdomain, ind_fig,U,uxinterp,uyinterp);
    end
    
    %size(rho2_store)
    
    if ~no_flux && ~do_flow
        ni = 4;
        [xi2,yi2] = findInterpsVert(ni,meshdomain,xi(1),yi(1));
        rho_temp = ...
            getRhoGrid(xi2,yi2,d,r_e, r_c, mu, epsi, rho0, ...
            meshdomain, ind_fig,U,uxinterp,uyinterp);
        rho2_store = [rho2_store,rho_temp];
        [xi2,yi2] = findInterpsTop( ni,xi(1),yi2(end) );
        rho_temp = ...
            getRhoGrid(xi2,yi2,d,r_e, r_c, mu, epsi, rho0, ...
            meshdomain, ind_fig,U,uxinterp,uyinterp);
        rho2_store = [rho2_store,rho_temp];
        [xi2,yi2] = findInterpsTop( 1,xi(1)/10,-yi2(end) );
        rho_temp = ...
            getRhoGrid(xi2,yi2,d,r_e, r_c, mu, epsi, rho0, ...
            meshdomain, ind_fig,U,uxinterp,uyinterp);
        rho2_store = [rho2_store,rho_temp];
    end
    
    % work back from collector
    if do_flow
        rho3_store = ...
            getRhoC(d, r_e, r_c, mu, epsi, rho0, ...
            meshdomain, ind_fig,U,rhoc_store,uxinterp,uyinterp);
    else
        rho3_store = [];
    end
    
    % put together space charge interpolation vector
    rho_store = [rho_store,rho2_store,rho3_store];
    
    % interpolate and fill FEM domain with space charge
    [f_new,finterp] = interpf( p,t,rho_store,epsi,meshdomain,dummy_ind_fig );
    
    if plot_debug
        %figure(ind_fig)
        subplot(313)
        % visualize space charge
        pdeplot(p,e,t,'xydata',finterp(p(1,:),p(2,:)),'colormap','jet')
        axis equal
        title('Space Charge Interp')
        drawnow
    end
    %ind_fig = ind_fig + 1;
    
    %% find electric field at emitter and collector
    
    % emitter
    [ Exe,Eye ] = getE_alt( e_inds, grad_phi );
    
    % compute outer loop exit metric
    Ee(i) = sum(sqrt(Exe.^2 + Eye.^2).*dAe)/sum(dAe); % 'average' electric field
    
    % collector
    [ Exc,Eyc ] = getE_alt( c_inds, grad_phi );
    Ec = sqrt(Exc.^2 + Eyc.^2);
    
    % emitter, space charge
    rhoe = rho0; % assumed constant at emitter (may change this later)
    
    % collector, space charge
    rhoc  = getRhoEnd( xc,yc,rho_store,finterp );
    
    % determine currents at emitter and collector
    Ie(i) = sum((Exe.^2 + Eye.^2).^(0.5)*mu.*rhoe.*dAe);
    %Ie(i) = rhoe*Ee(i)*mu*2*pi()*r_e;
    Ic(i) = sum((Exc.^2 + Eyc.^2).^(0.5)*mu.*rhoc.*dAc);
    
    %% alternate method to determine charge conservation
    
    [Jin, Jout, Jnet, Jdiff] = findJConv(r_c,finterp,U,mu,1000,[0,d/2],d/2,uxinterpf,uyinterpf);
    Jd(i) = Jdiff;
    
    %% alternate method to determine emitter and collector currents
    
    [~, Jout, ~, ~] = findJConv(r_c,finterp,U,mu,1000,[0,d+r_e],r_e*5,uxinterpf,uyinterpf);
    Je(i) = Jout;
    
    [Jin, ~, ~, ~] = findJConv(r_c,finterp,U,mu,1000,[0,-r_c],r_c*5,uxinterpf,uyinterpf);
    Jc(i) = Jin;
    
    %% re-solve for potential field: Davis and Hoburg, Step 4
    
    if no_flux
        [phi_new, ~] = pdenonlin(b,p,e,t,c,a,f_new);
    else
        [K,M,F,Q,G,H,R]=assempde(b,p,e,t,c,a,f_new);
        K(pbdry,pbdry) = K(pbdry,pbdry)+Sext(1:length(pbdry),1:length(pbdry));
        phi_new = assempde(K,M,F,Q,G,H,R);
    end
    
    %% determine exit criteria, exit if met: Davis and Hoburg, Step 5
    
    exit_test = norm(phi_new-phi_last)/norm(phi_last);
    
    res_norm(i) = exit_test*100;
    
    res_norm(end);
    
    % reset necessary variables
    omega = .5;
    phi_last = phi_new;
    phi = phi_new*omega + phi*(1-omega);
    [Ex_new, Ey_new] = pdegrad(p,t,phi);
    grad_phi = [Ex_new; Ey_new];
    
    %% determine EHD thrust
    
    % compute thrust
    T = sum(f_new*epsi.*Ey_new.*dAFEM);
    TA = T/(2*pi()*r_e);
    
    % find T/P
    P = (kV(1)-kV(2))*1e3*Je(i);
    eta(i) = T/P*1e3;
    eta_1d = 1e3/(mu*(kV(1)-kV(2))*1e3/d + U);
    
    % space charge distribution
    x = linspace(meshdomain(1),meshdomain(2),101);
    y = linspace(meshdomain(3),meshdomain(4),101);
    [X,Y] = meshgrid(x,y);
    f_plot = finterp(X,Y);
    
    % calculate more accurate version of 1D T/P
    eta_alt = getTP(d,mu,kV,Ec,rhoc,yc,dAc,r_e);
    
    %print stuff
    fprintf("Thrust: %.3f \n", T)
    fprintf("Power: %.3f \n", P)
    fprintf("Thrust-to-Power Ratio: %.3f \n", eta(i))
    fprintf("Emmitter E Field:  %.3e \n", Ee(i))
    fprintf("Emitter Current:   %.3g,  %.3g \n", Ie(i), Je(i))
    fprintf("Collector Current: %.3g,  %.3g \n", Ic(i), Jc(i))
    fprintf("Rho0: %.3e \n",rho0)
    fprintf("Exit Test: %.2f perc \n",100*exit_test)
    % increment iteration variables
    i = i + 1;
    
    % set maximum number of iterations
    if i > max_iter
        break
    end
    
end
    

end









