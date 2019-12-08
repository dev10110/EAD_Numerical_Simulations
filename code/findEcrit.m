function [ E_crit ] = findEcrit(p,e,t,d,r_e,r_c)

% get boundary conditions through custom function
b = @pdebound_incept;
a = 0;
c = 1;

f = 0;
% solve for potential field
[phi, ~] = pdenonlin(b,p,e,t,c,a,f);

% find electric field at emitter
[Ex, Ey] = pdegrad(p,t,phi);
grad_phi = [Ex; Ey];

[pts_e, ~, e_inds, ~] = getBPts(p,e,t,d,r_e,r_c);

[~,~,dAe_alt] = getMOCStart_alt( 25, d, r_e, pts_e, e_inds);
[ Exe,Eye ] = getE_alt( e_inds, grad_phi );
E_crit = sum(sqrt(Exe.^2 + Eye.^2).*dAe_alt)/sum(dAe_alt);

end

