function [xe_new,ye_new,theta_full] = getMOCStartDev( num_emit, d, r_e, theta_full)

theta_new = [linspace(pi/2,0.5005*pi,num_emit+1), linspace(0.5005*pi,3*pi/2,num_emit+1)];

%removes the ones at pi/2 and 3pi/2
N = length(theta_new);

theta_new = theta_new(1:(N-1));

theta_new = setdiff(theta_new,theta_full);



xe_new = r_e*cos(theta_new);
ye_new = (d + r_e) + r_e*sin(theta_new);


if length(theta_full) == 0
    theta_full = theta_new;
    
else
    theta_full = union(theta_new,theta_full);
end



end
