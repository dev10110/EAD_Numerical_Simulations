function results = convergent_solver(max_rho0_guess, geometry, constants, plot_debug)

global kV

r_c = geometry.r_c;
r_e = geometry.r_e;
d = geometry.d;
meshdomain = geometry.meshdomain;

%% generate geometry
num_ref = 3;
[p, e, t] = getGeoTD_1D(d, 1, r_c, r_e, meshdomain,num_ref);

%% determine critical e field using peeks criteria:
pr = 76; % cm Hg
T = 298; % K
delt = 3.92*pr/T;
E_crit = 30*delt*(1 + 0.301/sqrt(delt*r_e*100))*1e5; % V/m


%% use fzero to find the appropriate rho0
if 1
    diff_func = @(rho0) find_emitter_eField(rho0, E_crit, p,e,t, geometry, constants, plot_debug);
    
    
    if 1
        options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval},'TolX',1e-7); % show iterations
    else
        options = optimset('Display','off','TolX',1e-7);
    end
    
    if length(max_rho0_guess) == 1
        rho0 = fzero(diff_func, [0, max_rho0_guess], options);
    elseif length(max_rho0_guess) == 2
        rho0 = fzero(diff_func, [max_rho0_guess(1), max_rho0_guess(2)], options);
    end
end


%% solve the problem using this rho0

results = get_results_using_rho(rho0, p,e,t, geometry, constants, plot_debug);

results.mesh.p = p;
results.mesh.e = e;
results.mesh.t = t;

results.E_crit = E_crit;
results.geometry = geometry;

results.kV_emitter = kV(1);
results.kV_collector = kV(2);
end



