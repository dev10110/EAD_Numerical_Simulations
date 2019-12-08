function results = get_results_using_rho(rho0, p,e,t, geometry, constants, plot_debug)
%% almost the same, but now gives the full results instead.

global kV

r_c = geometry.r_c;
r_e = geometry.r_e;
d = geometry.d;
meshdomain = geometry.meshdomain;

collector_top_fun = @(x)  sqrt(r_c^2 - x.^2) - r_c;
collector_bot_fun = @(x) -sqrt(r_c^2 - x.^2) - r_c;
emitter_top_fun = @(x)  sqrt(r_e^2 - x.^2) + d + r_e;
emitter_bot_fun = @(x) -sqrt(r_c^2 - x.^2) + d + r_e;



% get area of each FEM element
dAFEM = getdAFEM( p,t );

mu = constants.mu;
epsi = constants.epsi;

b = @pdebound_1d;
a = 0;
c = 1;

num_emit_default = 20;

%get the emitter and collector points
[pts_e, pts_c, e_inds, c_inds] = getBPts(p,e,t,d,r_e,r_c);
[xc,yc,dAc,c_thetas,c_inds] = getMOCEnd( num_emit_default,r_c,pts_c,c_inds);
[xe_alt,ye_alt,dAe,e_thetas,e_inds] = getMOCStart_alt( num_emit_default, d, r_e,pts_e,e_inds); % for current calculation




%% solve for the phi
f = 0; %assume 0 space charge at first

[phi, ~] = pdenonlin(b,p,e,t,c,a,f);

% get area of each FEM element
dAFEM = getdAFEM( p,t );

% find electric field
[Ex, Ey] = pdegrad(p,t,phi);
grad_phi = [Ex; Ey];

% initialize potential field storage variables
phi_last = phi;


[ Exe,Eye ] = getE_alt( e_inds, grad_phi );
[ Exc,Eyc ] = getE_alt( c_inds, grad_phi );

%start it high
exit_test = 1;

while exit_test > 0.005
    
    if plot_debug
        figure;
        subplot(221)
        pdeplot(p,e,t,'XYData',phi_last)
        %scatter(p(1,:),p(2,:),[],phi_last,'.');
        colorbar
        title('Potential')
        axis equal
        drawnow
    end
    
    % find interpolation surface for electric field
    [ uxinterp,uyinterp ] = getEInterp( p,t,grad_phi,'All' );    
    
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
    while (touches_all_bounds == 0)
        
        %check if the space charge lines touch all bounds:
        %check that the min to max of x is within 99% of the full domain:
        if ((max(rho_store(1,:))-min(rho_store(1,:))) > 0.99*(meshdomain(2)-meshdomain(1))) && ((max(rho_store(2,:))-min(rho_store(2,:))) > 0.99*(meshdomain(4)-meshdomain(3)))
            touches_all_bounds = 1;
        else
            %if it isnt, double the number of exit points  and check
            %again
            num_emit = round(num_emit*2);
            fprintf(', %i',num_emit);
            [xe_new,ye_new,theta_full] = getMOCStartDev( num_emit, d, r_e, theta_full);
            
            
            [rho_store_temp, ~] = ...
                getRho(xe_new, ye_new, d, r_e, r_c, mu, ...
                epsi, rho0, meshdomain,uxinterp,uyinterp);
            
            rho_store_temp2 = [rho_store rho_store_temp];
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
    
    results.num_emit = num_emit
    
    % interpolate and fill FEM domain with space charge
    [f_new,finterp] = interpf( p,t,rho_store,epsi);
    %note, f_new(point) = rho(point)/epsi
    %finterp(x,y) = rho(x,y)
    
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
    results.Ee = Ee;
    
    
    % collector
    [ Exc,Eyc ] = getE_alt( c_inds, grad_phi );
    Ec = sqrt(Exc.^2 + Eyc.^2);
    results.Ec = Ec;
    
    % emitter, space charge
    rhoe = rho0; % assumed constant at emitter (may change this later)
    
    % collector, space charge
    rhoc  = getRhoEnd( xc,yc,rho_store,finterp );
    
    % determine currents at emitter and collector
    Ie_old = sum((Exe.^2 + Eye.^2).^(0.5)*mu.*rhoe.*dAe);
    results.Ie_old= Ie_old;
    
    dIe = (Exe.*cos(e_thetas) + (Eye.*sin(e_thetas))).*mu.*rhoe.*r_e;
    Ie  = trapz(e_thetas,dIe);
    results.Ie = Ie;
    
    Ic_old = sum((Exc.^2 + Eyc.^2).^(0.5).*mu.*rhoc.*dAc);
    results.Ic_old = Ic_old;
    
    dIc = ((Exc.*cos(c_thetas)) + (Eyc.*sin(c_thetas))).*mu.*rhoc.*r_c;
    Ic  = trapz(c_thetas,dIc);
    results.Ic = Ic;
    
    %% alternate method to determine emitter and collector currents
    
    if 0 %if you want to calculate Je
        [ uxinterpf,uyinterpf ] = getEInterpFull( p,t,grad_phi );
        [~, Jout, ~, ~] = findJConv(r_c,finterp,mu,1000,[0,d+r_e],r_e*5,uxinterpf,uyinterpf);
        Je = Jout;
        results.Je = Je;

        [Jin, ~, ~, ~] = findJConv(r_c,finterp,mu,1000,[0,-r_c],r_c*5,uxinterpf,uyinterpf);
        Jc = Jin;
        results.Jc = Jc;
    else
        results.Je = 0;
        results.Jc = 0;
    end
    
    
    
    
    %% re-solve for potential field: Davis and Hoburg, Step 4
    
    [phi_new, ~] = pdenonlin(b,p,e,t,c,a,f_new);
    
    %% determine exit criteria, exit if met: Davis and Hoburg, Step 5
    
    exit_test = norm(phi_new-phi_last)/norm(phi_last);
    
    fprintf('\n rho0: %.3g, Ee: %.3g, exit: %.1f perc',rho0, Ee, exit_test*100);
    
    %% reset necessary variables
    omega = .5;
    phi_last = phi_new;
    phi = phi_new*omega + phi*(1-omega);
    
    [Ex_new, Ey_new] = pdegrad(p,t,phi);
    grad_phi = [Ex_new; Ey_new];
    
    %% thrust and power
    T_old = sum(f_new*epsi.*Ey_new.*dAFEM)
    results.T_old = T_old;
    
    
    T = calculateThrust_using_interpolants(finterp,uyinterp,r_e,r_c,meshdomain,collector_bot_fun,collector_top_fun,emitter_bot_fun,emitter_top_fun)
    
    results.T = T;
    
    if plot_debug
        figure;
        subplot(121)
        pdeplot(p,e,t,'XYData',f_new*epsi.*Ey_new.*dAFEM)
        colorbar
        title('thrust per unit area (old)')
        axis equal

        subplot(122)
        [Xq,Yq] = meshgrid(meshdomain(1):0.1:meshdomain(2),meshdomain(3):0.1:meshdomain(4));
        Vq =finterp(Xq,Yq).*uyinterp(Xq,Yq);
        [~,h]=contourf(Xq,Yq,Vq,100);
        set(h,'LineColor','none')
        title('thrust per unit area (interpolant)');       
        colorbar       
        axis equal
        
        drawnow;
        
    end
    
    
    P = (kV(1)-kV(2))*1e3*Ie;
    results.P = P;
    
    results.rho0 = rho0;
    
    results.eta = T/P*1e3;
    
    
    
    results.p = p;
    results.e = e;
    results.t = t;
    results.phi = phi;
    results.f_new = f_new;
    results.finterp=finterp;
    results.Ey_new = Ey_new;
    results.dAFEM = dAFEM;
    results.uxinterp = uxinterp;
    results.uyinterp = uyinterp;
    
    
end
