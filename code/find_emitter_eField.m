function diff = find_emitter_eField(rho0, Ecrit, p,e,t, geometry, constants, plot_debug)


%solve the problem using the guessed rho0
results = get_results_using_rho(rho0, p,e,t, geometry, constants, plot_debug);

diff = results.Ee - Ecrit;
disp('diff')
disp(diff)
disp(' ')
