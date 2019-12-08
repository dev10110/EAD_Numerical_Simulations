function [p, e, t, b] = getGeoTD_1D(d, V, r_c, r_e, meshdomain,num_ref)
% Single Stage E2 thruster. Opens PDE Solver GUI.
% Inputs: d: gap length (cm)
%         V: emitter voltage (kV)

%% Initialize PDE Solver
[pde_fig, ax] = pdeinit;
pdetool('appl_cb',5);                           % Electrostatics
pdetool('snapon','on');
set(ax,'DataAspectRatio',[1 0.75 1]);
set(ax,'PlotBoxAspectRatio',[300 200 1]);
set(ax,'XLim',[meshdomain(1) meshdomain(2)]);
set(ax,'YLim',[meshdomain(3) meshdomain(4)]);
set(ax,'XTick',meshdomain(1):50:meshdomain(2));
set(ax,'YTick',meshdomain(3):50:meshdomain(4));

%% Set Geometry
pderect(meshdomain, 'R1');
pdecirc(0, (d + r_e), r_e, 'C1');
pdecirc(0, -r_c, r_c, 'C2');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1-C1-C2');

%% Set Boundary Conditions
pdetool('changemode',0)

%% Solve in PDE Toolbox
% Generate Mesh
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
pdetool('initmesh')
for i = 1:num_ref
    pdetool('refine')
end
pdetool('jiggle')

h = findobj(get(pde_fig,'Children'),'flat','Tag','PDEMeshMenu');

hp = findobj(get(h,'Children'),'flat','Tag','PDEInitMesh');
he = findobj(get(h,'Children'),'flat','Tag','PDERefine');
ht = findobj(get(h,'Children'),'flat','Tag','PDEMeshParam');

p = get(hp,'UserData');
e = get(he,'UserData');
t = get(ht,'UserData');

hbound = findobj(get(pde_fig,'Children'),'flat','Tag','PDEBoundMenu');
hb_cld = get(hbound, 'children');
b = get(hb_cld(7), 'userdata');

end














