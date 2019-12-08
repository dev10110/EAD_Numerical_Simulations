clear all; close all;clc;


ref = load('solution_70kv_25cm_mesh_3_tolx_1e7_efield_0R5P.mat')
ref(2) = load('solution_70kv_25cm_mesh_4_tolx_5e8_efield_0R25P.mat')

%%


for i = 1:2
    num_of_elements(i) = size(ref(i).full_results.t,2);
    num_of_emit(i) = ref(i).full_results.num_emit;
    Ee(i) = ref(i).full_results.Ee;
    Ie(i) = ref(i).full_results.Ie;
    Ic(i) = ref(i).full_results.Ic;
    T(i)  = ref(i).full_results.T;
    rho0(i) = ref(i).full_results.rho0;
    
end

%%
numdiff = 100*(num_of_elements(2)-num_of_elements(1))/num_of_elements(1)
Eediff = 100*(Ee(2)-Ee(1))/Ee(1)
Iediff = 100*(Ie(2)-Ie(1))/Ie(1)
Icdiff = 100*(Ic(2)-Ic(1))/Ic(1)
 Tdiff = 100*(T(2)-T(1))/T(1)
rho0diff = 100*(rho0(2)-rho0(1))/rho0(1)