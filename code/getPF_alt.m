function [ u ] = getPF_alt( x, y, r_c, U )

%% transform x, y to appropriate coordinates

y_new = x;
x_new = -y + 0.023557386828733;

y = y_new;
x = x_new;

[u,v] = getPotFlowPt_alt(r_c,U,x,y);

%% change back to coordinates in FEM frame

u_new = v;
v_new = -u;

u = u_new;
v = v_new;

u = [u;v];

end




















