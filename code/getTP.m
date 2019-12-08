function eta_alt = getTP(d,mu,kV,Ec,rhoc,yc,dAc,r_e)

d_temp = d - yc;
eta_temp = d_temp./(mu*(kV(1)-kV(2)));
T = Ec.*rhoc.*dAc.*d_temp;

eta_alt = sum(T.*eta_temp)/sum(T);
%eta_alt = sum(rhoc.*d_temp.*eta_temp)/sum(rhoc.*d_temp);
%eta_alt = sum(d_temp.*eta_temp)/sum(d_temp);

end

