function [xe,ye,dA,theta,inds] = getMOCStart_alt( num_emit, d, r_e, pts,inds)

% offset = 0;
% 
% x_e = 0;
% y_e = d + r_e;
% 
% r = r_e*(1 + offset);
% theta = linspace(0,2*pi(),num_emit*4) + pi()/2;
% 
% xe = r*cos(theta) + x_e;
% ye = r*sin(theta) + y_e;

xe = pts(1,:);
ye = pts(2,:);

for i = 1:length(xe)
    
    x = xe(i);
    y = ye(i) - (d + r_e);
    
    theta(i) = atan2(y,x);
    
%     if y >= 0 && x >= 0
%         theta(i) = temp;
%     elseif y <= 0 && x >= 0
%         theta(i) = 2*pi() + temp;
%     elseif y >= 0 && x <= 0
%         theta(i) = pi() + temp;
%     else
%         theta(i) = pi() + temp;
%     end   
    
end

[theta, ind] = sort(theta);
xe = xe(ind);
ye = ye(ind);
inds = inds(ind);

% find area element

for i = 1:length(theta)
    
    if i == 1
        ind1 = length(theta);
        ind2 = i + 1;
    elseif i == length(theta)
        ind1 = i - 1;
        ind2 = 1;
    else
        ind1 = i-1;
        ind2 = i+1;
    end
    
    theta1 = theta(ind1);
    theta2 = theta(ind2);    
    thetam = theta(i);
    
    if i == 1
        theta1 = theta1 - 2*pi();
    elseif i == length(theta)
        theta2 = theta2 + 2*pi();
    end
    
    temp1 = (thetam + theta1)/2;
    temp2 = (thetam + theta2)/2;
    
    dA(i) = (temp2-temp1)*r_e;     
    
end

end

