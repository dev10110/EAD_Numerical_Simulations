function T = calculateThrust_using_interpolants(rho_interp,Ey_interp,r_e,r_c,meshdomain,collector_bot_fun,collector_top_fun,emitter_bot_fun,emitter_top_fun)
%% calculates the thrust by computing the double integral from the numerical result

%% get integrand

fun = @(x,y) rho_interp(x,y).*Ey_interp(x,y);

%% first integrate over the full domain


Tfull = integral2(fun,meshdomain(1),meshdomain(2),meshdomain(3),meshdomain(4));

%% integrate over the collector


Tcollector = integral2(fun,-r_c, r_c, collector_bot_fun,collector_top_fun);



%% integrate over the emitter


Temitter = integral2(fun,-r_e, r_e, emitter_bot_fun,emitter_top_fun);



%% compute results

T = Tfull - Tcollector - Temitter;
