function [ Sext, pbdry ] = Ballooningv2( meshstruct,bdry, alpha, m )
%BALLOONING calculates stiffness matrix of boundary to ballooned elements
%   Input:
%   meshstruct - data structure containing p,e,t PDE toolbox mesh data
%   bdry - array containing list of boundary edge labels
%   alpha - ballooning geometric ratio
%   m - number of ballooned layers
%   Output:
%   Sext - ballooned stiffness matrix
%   pbdry - index numbers of points on outer boundary

%%

%constants
epsi = 8.85418782e-12;

%pull out struct data
p = meshstruct.p;
e = meshstruct.e;
t = meshstruct.t;

%record info on the original mesh
Np = size(p,2);
Ne = size(e,2);
Nt = size(t,2);

%% Generate first layer of ballooned elements

% reduce edge matrix to outer boundary edges
% 1 and 2 refer to edge end points (references for p matrix)
% 5 is subdomain reference
etemp = [e([1 2 5],:) ; 1:Ne];
ered = [];

% pull out edge coordinates for those edges that lie on the outer boundary
for ii=1:length(bdry)
    ered = [ered etemp(:,etemp(3,:)==bdry(ii))];
end
% number of edges on outer boundary
Nered = size(ered,2);
clear etemp

% calculate ballooned coordinates of edge nodes

% eliminate redundant boundary points
pbdry = unique([ered(1,:) ered(2,:)]);
% pull out x/y coordinates of boundary nodes
pnew = [p(:,pbdry) alpha.*p(:,pbdry); pbdry pbdry];
% get number of boundary nodes
Npbdry = length(pbdry);

% generate triangle connectivity
tnew = zeros(4,Nered*2);
ptemp = 1:size(pnew,2);

% each edge on the boundary results in two triangles
for ii=1:Nered
    % find node number of point on outer domain
    p2 = ptemp(pbdry == ered(1,ii));
    % find node number of point on indder domain
    p1 = ptemp(pbdry == ered(2,ii));    
    
    t1 = [p1 p2+Npbdry p1+Npbdry 1];
    t2 = [p1 p2        p2+Npbdry 1];
    
    % build new triangle point matrix
    tnew(:,(2*ii-1):2*ii) = [t1' t2'];
end
Ntnew = size(tnew,2);
clear p1 p2 t1 t2

%% Calculate stiffness matrix of balloned layer

% do some pre-calcs
x1 = pnew(1,tnew(1,:));
x2 = pnew(1,tnew(2,:));
x3 = pnew(1,tnew(3,:));

y1 = pnew(2,tnew(1,:));
y2 = pnew(2,tnew(2,:));
y3 = pnew(2,tnew(3,:));

% populate stiffness matrix
Sext = zeros(2*Npbdry,2*Npbdry);

for ii=1:Ntnew
    P = [1 x1(ii) y1(ii);1 x2(ii) y2(ii);1 x3(ii) y3(ii)];
    C = inv(P);
    A = abs(det(P))/2;
    grad = C(2:3,:);
    s = A*(grad)'*(grad);
    indx = tnew(1:3,ii);
    Sext(indx,indx) = Sext(indx,indx)+s;
end
Sext = sparse(Sext);

%% Compute stiffness matrix for m ballooned layers

Sext = SextCompute(Sext,m,Npbdry);

end

function Sm = SextCompute(Sext1,m,Npbdry)

%initialize current stiffness matrix
Scur = Sext1;

%for each annuli, compute reduced matrix
for ii = 2:m
    Scur = ScurCompute(Scur,Npbdry);
end

Sm = Scur;
end

function Scur = ScurCompute(S,Npbdry)
%SCURCOMPUTE calculates next reduced stiffness matrix given current
%reduced stiffness matrix and current annuli stiffness matrix
%   Input:
%   S - current reduced stiffness matrix
%   T - current annuli stiffness matrix
%   Npbdry - number of nodes on domain outer boundary
%   Output:
%   Scur - new current reduced stiffness matrix

%%

%extract quadrants of stiffness matrix
S11 = S(1:Npbdry,1:Npbdry);
S12 = S(1:Npbdry,(1+Npbdry):2*Npbdry);
S21 = S((1+Npbdry):2*Npbdry,1:Npbdry);
S22 = S((1+Npbdry):2*Npbdry,(1+Npbdry):2*Npbdry);
    
Scur = [S11-S12*((S22+S11)\S21)  -S12*((S22+S11)\S12); ...
    -S21*((S22+S11)\S21)     S22-S21*((S22+S11)\S12)];

end
















