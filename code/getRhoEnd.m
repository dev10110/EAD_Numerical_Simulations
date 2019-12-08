function [ rhoc ] = getRhoEnd( xc,yc,rho_store,F )

%rhoc = zeros(1,length(xc));

%rho = TriScatteredInterp(rho_store(1,:)',rho_store(2,:)',rho_store(3,:)');

% for i = 1:length(xc)
%     
%     rhoc(i) = rho(xc(i),yc(i));
%     
%        
% end

rhoc = F(xc,yc);

end

