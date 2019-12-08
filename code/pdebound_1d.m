function [qmatrix,gmatrix,hmatrix,rmatrix] = pdebound_1d(p,e,u,time)

global kV no_flux

ne = size(e,2); % number of edges
qmatrix = zeros(1,ne);
gmatrix = qmatrix;
hmatrix = zeros(1,2*ne);
rmatrix = zeros(1,2*ne);

Ve = kV(1)*1000; %emitter voltage
Vc = kV(2)*1000; %collector voltage

for k = 1:ne
    switch e(5,k)
        case {1,2,3,4}
%             if ~no_flux
%                 hmatrix(k) = 1;
%                 hmatrix(k+ne) = 1;
%                 rmatrix(k) = 0;
%                 rmatrix(k+ne) = 0;
%             end
            hmatrix(k) = 1;
            hmatrix(k+ne) = 1;
            rmatrix(k) = 0; %define the voltage of the bounding box to be 0
            rmatrix(k+ne) = 0; %define the voltage of the bounding box to be 0
        case {5,6,7,8}
            hmatrix(k) = 1;
            hmatrix(k+ne) = 1;
            rmatrix(k) = Ve;
            rmatrix(k+ne) = Ve;
        case {9,10,11,12,13}
            hmatrix(k) = 1;
            hmatrix(k+ne) = 1;
            rmatrix(k) = Vc;
            rmatrix(k+ne) = Vc;
    end
end