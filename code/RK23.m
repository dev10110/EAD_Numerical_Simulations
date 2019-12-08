function [ x3, terr, hopt, h_old, prime] = RK23( x0,t0,h,f,err_max )
%Adaptive Runge-Kutta Fehldberg
%  x0 - initial value
%  t0 - initial time
%  h - timestep
%  f - derivative function of form x' = f(t,x)
%  x1 - solution value at t+h
%  terr - truncation error estimate

%populate coefficients

h_old = h;

a = [0 .25 27/40 1];
b = [0       0       0; ...
    1/4      0       0; ...
    -189/800 729/800 0; ...
    214/891  1/33    650/891];

c = [214/891 533/2106; ...
    1/33     0; ...
    650/891  800/1053; ...
    0        -1/78];

K = zeros(4,size(x0,2),size(x0,1));

K(1,:,:) = h.*f(t0,x0); 

K(2,:,:) = h.*f(t0+a(2).*h, x0 + reshape(b(2,1).*K(1,:,:),size(x0))); 

K(3,:,:) = h.*f(t0+a(3).*h, x0 + reshape(b(3,1).*K(1,:,:)+b(3,2).*K(2,:,:),size(x0))); 

K(4,:,:) = h.*f(t0+a(4).*h, x0 + reshape(b(4,1).*K(1,:,:)+b(4,2).*K(2,:,:)+b(4,3).*K(3,:,:),size(x0))); 

prime = [sum(repmat(c(:,2),1,size(x0,2)).*K(:,:,1),1); ...
    sum(repmat(c(:,2),1,size(x0,2)).*K(:,:,2),1)]./h;

x2 = x0+[sum(repmat(c(1:3,1),1,size(x0,2)).*K(1:3,:,1),1); ...
    sum(repmat(c(1:3,1),1,size(x0,2)).*K(1:3,:,2))];
x3 = x0+h.*prime;

%adjust step size for error
terr = max(abs(x3-x2),[],1);
beta = .9;
hopt = h;


ind = terr>err_max;
hopt(ind) = beta.*hopt(ind).*(err_max./terr(ind))^(1/3);

if max(ind)==1
    %disp('re-evaluating for better accuracy')
    [ x3(:,ind), terr(ind), hopt(ind), h_old(ind), prime(:,ind)] = RK23( x0(:,ind),t0(ind),hopt(ind),f,err_max );
end

ind = terr < 1*err_max;
% if max(ind) == 1;
%     disp('increasing step size')
%     
% end
hopt(ind) = beta.*hopt(ind).*(err_max./terr(ind))^(1/4);

end


