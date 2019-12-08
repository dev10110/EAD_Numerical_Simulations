function rho_store = getRho_efficient( d, r_e, r_c,  ...
    mu, epsi, rho0, meshdomain,uxinterp,uyinterp)


num_emit = 40;

%create angles to test
theta = linspace(pi/2,3*pi/2,num_emit+2);
theta = theta(2:(end-1));

endpos = zeros(size(theta));

rho_store = [];

%generate the basic characteristic lines
for i=1:length(theta)
    [rho_store_temp, endpos(i)] = generateCharacteristic(theta(i), d, r_e, r_c, mu, epsi, rho0, meshdomain,uxinterp,uyinterp);
    rho_store=[rho_store,rho_store_temp];
end
%note: end pos can have 4 values:
%1 if touching top boundary
%2 if touching left boundary
%3 if touching bottom boundary
%4 if touching collector

%remove duplicates from rho_store
rho_store = (unique(rho_store','rows'))';

%iterate starting points if boundary is not touched
while ismember(1, endpos) == 0 %not touching top boundary
    %need to increase resoution before first theta
    theta = [theta(1)/2, theta];
    
    [rho_store_temp, endpos_temp] = generateCharacteristic(theta(1), d, r_e, r_c, mu, epsi, rho0, meshdomain,uxinterp,uyinterp);
    
    rho_store=[rho_store,rho_store_temp];
    
    endpos=[endpos_temp,endpos];
    
end

%check left boundary
while ismember(2,endpos) == 0
    
    %find the last position with endpos = 1
    theta_min_ind = find(endpos < 2,1,'last');
    theta_min = theta(theta_min_ind);
    theta_max = theta(theta_min_ind+1);
    
    theta_new = (theta_max + theta_min)/2;
    
    [rho_store_temp, endpos_temp] = generateCharacteristic(theta_new, d, r_e, r_c, mu, epsi, rho0, meshdomain,uxinterp,uyinterp);
    
    %restore the variables
    endpos = [endpos(1:theta_min_ind), endpos_temp, endpos((theta_min_ind+1):end)];
    theta = [theta(1:theta_min_ind), theta_new, theta((theta_min_ind+1):end)];
    rho_store=[rho_store,rho_store_temp];
    
    
end

%check bottom boundary
while ismember(3,endpos) == 0
    
    %find the last position with endpos = 1
    theta_min_ind = find(endpos < 3,1,'last');
    theta_min = theta(theta_min_ind);
    theta_max = theta(theta_min_ind+1);
    
    theta_new = (theta_max + theta_min)/2;
    
    [rho_store_temp, endpos_temp] = generateCharacteristic(theta_new, d, r_e, r_c, mu, epsi, rho0, meshdomain,uxinterp,uyinterp);
    
    %restore the variables
    endpos = [endpos(1:theta_min_ind), endpos_temp, endpos((theta_min_ind+1):end)];
    theta = [theta(1:theta_min_ind), theta_new, theta((theta_min_ind+1):end)];
    rho_store=[rho_store,rho_store_temp];
    
    
end

rho_store = (unique(rho_store','rows'))';

%remove the all zero col:
zeroxind = find(rho_store(1,:) == 0);
zeroyind = find(rho_store(2,:) == 0);

zeroxyind = intersect(zeroxind,zeroyind);

rho_store(:,zeroxyind) = [];

%now it should be touching all four sides

%reflect the rho_store

rho_store2(1,:) = -rho_store(1,:);
rho_store2([2,3],:) = rho_store([2,3],:);

rho_store = [rho_store, rho_store2];

end

