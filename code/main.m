%% Version description

% v12.1: 
% 1. RKF23 adaptive time integration implemented in getRho

%% Do outer iteration loop

%% setup exit criteria for outer iteration loop and some outer loop and
% inner loop options
if ~do_newton & ~do_basic_rho_update 
    max_out_it = 1;
end

% setup parameters for newton step
pert = 0.01;

first_it = 1;
temp = 0;
phi = 0;

% setup storage variables
Ee_store = zeros(max_out_it,max_iter);
eta_store = zeros(max_out_it,max_iter);
Ie_store = zeros(max_out_it,max_iter);
Ic_store = zeros(max_out_it,max_iter);
res_store = zeros(max_out_it,max_iter);

%% Do loop

if do_parallel
    matlabpool open 2
end

for iii = 1:max_out_it
    
    iii
    
    %% Do inner iteration loop    
    rho0 = rho0_now;
    tic
    [eta_1d, eta_alt, eta, Ee, Ie, Ic, Jd, Je, Jc, res_norm, temp_new, phi_new, f_plot, phi_lap] = ...
        EHDsolverf(p,e,t,r_e,r_c,d,kV,U,meshdomain,rho0,max_iter,first_it,phi,...
        Sext,pbdry);
    toc
    Ee1 = Ee(end);
    
    % store variables
    ns = length(Ee);
    Ee_store(iii,1:ns) = Ee;
    eta_store(iii,1:ns) = eta;
    Ie_store(iii,1:ns) = Ie;
    Ic_store(iii,1:ns) = Ic;
    res_store(iii,1:ns) = res_norm;
    
    %% check to see if Peek's Law criteria has been met: Davis and Hoburg, Step 6
    
    exit_out = abs(Ee1 - E_crit)/E_crit;
    if exit_out < out_thresh
        break
    end
    
    drawnow
    
    %% Do newton step if not sufficiently converged
    
    if do_newton
        
        rho0 = rho0_now*(1+pert);
        [~, ~, ~, Ee] = ...
            EHDsolverf(p,e,t,r_e,r_c,d,kV,U,meshdomain,rho0,max_iter,first_it,phi,...
            Sext,pbdry);
        Ee2 = Ee(end);
        
        dEdr = (Ee2 - Ee1)/(pert*rho0_now);
        %%% slow down factor:
        newtonsFactor = 0.6;
        %%%
        rho0_now = rho0_now - newtonsFactor*(Ee1 - E_crit)/dEdr;
        
        rho0_store(iii+1) = rho0_now;
        
        % reset starting point
        %         first_it = 0;
        %         phi = phi_new;
        
        if plot_debug
            close all;
        end
        
    end
    
    if do_basic_rho_update
        
        if Ee1 > E_crit %rho needs to increase
            rho0_now = rho0_now*1.2
        else
            rho0_now = rho0_now*0.8;
        end
        rho0_store(iii+1) = rho0_now;
        
        % reset starting point
        %         first_it = 0;
        %         phi = phi_new;
        
    end
        
end

if do_parallel
    matlabpool close
end






