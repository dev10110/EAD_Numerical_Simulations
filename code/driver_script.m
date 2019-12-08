clear all;
clc;
close all;

warning('off','all');

%% This is the standard driver template

%% Add in path for correct code directory

%addpath('../Version 14 copy - mesh convergence')

%% Set global variables

global kV

%% Guess at an initial charge density concentration (for now assume constant on emitter wire surface
% Davis and Hoburg, Step 2

% charge density
%max_rho0_guess_values = [3.40E-05	2.06E-05	1.38E-05	1.10E-05
%    6.12E-05	3.32E-05	2.09E-05	1.68E-05
%    8.37E-05	5.06E-05	3.24E-05	2.45E-05
%    1.19E-04	6.87E-05	4.17E-05	3.14E-05];%[2e-4 3e-4];

max_rho0_guess_values = [5.40E-05];
max_rho0_guess_values = [2e-05, 4e-5]

%%%%%%%%%%% define collector and emitter voltagesCase: 
kV_emitter_vals = [70]; %,50, 60,70,80]; % kV_emitter_vals specifies array of run cases 
kV_collector_vals = 0; 
%%%%%%%%%%%

plot_debug_flag = 0; % Produces extra plots for debugging

%% get list of solutions already found

files = dir('/Users/Devansh/OneDrive - Imperial College London/MIT/UROP/Electro-Static-Solver/Version 14 copy - mesh convergence/results_new');
%% initialize geometry and parameters

% electrode geometry [m]
r_c_vals = 0.0381;% %2.5*0.00635/2;
r_e_vals = 0.000254;%0.202/2*10^-3;
d_vals = 0.150:0.05:0.300;
d_vals = 0.250;

% solution domain
meshdomain = [-2.5/2,2.5/2,-2.5/2,2.5/2];


%create geometry vals
m=1;
for i = 1:length(r_c_vals)
    for j=1:length(r_e_vals)
        for k=1:length(d_vals)
            geometry_vals(m).r_c = r_c_vals(i);
            geometry_vals(m).r_e = r_e_vals(j);
            geometry_vals(m).d   = d_vals(k);
            geometry_vals(m).meshdomain = meshdomain;
            m=m+1;
        end
    end
end



% set other parameters
mu = 2.409638554e-4;%2*10^-4;
epsi = 8.85418782e-12;

constants.mu = mu;
constants.epsi = epsi;


m=1;

for j=1:length(kV_collector_vals)
    for k=1:length(geometry_vals)
        for i=1:length(kV_emitter_vals)
            close all force
            
            kV_emitter = kV_emitter_vals(i);
            
            kV_collector = kV_collector_vals(j);
            
            geometry = geometry_vals(k);
            
            kV = [kV_emitter, kV_collector];
            
            %check that the solution hasnt already been found.
            if ~isempty(find(strcmp({files.name}, sprintf('solution_%ikVe_%icm.mat',kV(1),round(geometry.d*100)))==1))
                continue
            end
            if kV_emitter == 60 && geometry.d == 0.250;
                continue
            end
            
            disp('***************')
            disp('Now solving case:')
            fprintf('kv_e = %i, d = %.3f \n',kV_emitter,geometry.d)
            disp(datetime('now'))
            
            try
                %full_results = convergent_solver([0.5*max_rho0_guess_values(i,k),1.3*max_rho0_guess_values(i,k)], geometry, constants, plot_debug_flag);
                full_results = convergent_solver(max_rho0_guess_values, geometry, constants, plot_debug_flag);
                filename = sprintf('results_new/solution_%ikVe_%icm.mat',kV(1),round(geometry.d*100));
                save(filename,'full_results');
            catch %allows the script to carry on to the next test case if the previous one failed (without throwing errors and stopping all the executions)
            end
            
            %full_results(m).kV_emitter = kV_emitter;
            %full_results(m).kV_collector = kV_collector;
            
            disp('DONE')
            disp(datetime('now'))
        end
    end
end






%results_summary = [kV_emmiter, kV_collector, d, T, P, T/P, Ee1, Je(end), Jc(end), rho0];













