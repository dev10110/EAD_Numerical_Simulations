function [qmatrix,gmatrix,hmatrix,rmatrix] = pdebound_incept(p,e,u,time)

global kVi

ne = size(e,2); % number of edges
qmatrix = zeros(1,ne);
gmatrix = qmatrix;
hmatrix = zeros(1,2*ne);
rmatrix = zeros(1,2*ne);

V = kVi*1000;

for k = 1:ne
    switch e(5,k)
        case {5,6,7,8} 
            hmatrix(k) = 1;
            hmatrix(k+ne) = 1;
            rmatrix(k) = V;
            rmatrix(k+ne) = V;
        case {9,10,11,12,13} 
            hmatrix(k) = 1;
            hmatrix(k+ne) = 1;
            rmatrix(k) = 0;
            rmatrix(k+ne) = 0;
    end
end