function results = solver_given_rho0(rho0, geometry, constants, plot_debug)

global kV

r_c = geometry.r_c;
r_e = geometry.r_e;
d = geometry.d;
meshdomain = geometry.meshdomain;

%% generate geometry
num_ref = 2;
[p, e, t] = getGeoTD_1D(d, 1, r_c, r_e, meshdomain,num_ref);




%% determine critical e field using peeks criteria:
pr = 76; % cm Hg
T = 298; % K
delt = 3.92*pr/T;
E_crit = 30*delt*(1 + 0.301/sqrt(delt*r_e*100))*1e5; % V/m

%% solve the problem using this rho0

results = get_results_using_rho(rho0, p,e,t, geometry, constants, plot_debug);

results.E_crit = E_crit;
results.geometry = geometry;

results.kV_emitter = kV(1);
results.kV_collector = kV(2);
end



