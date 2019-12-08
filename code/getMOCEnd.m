function [xc,yc,dA,theta,inds] = getMOCEnd( num_emit, r_c, pts, inds)

% offset = 0;
% 
% x_c = 0;
% y_c = -r_c;
% 
% r = r_c*(1 + offset);
% theta = linspace(0,2*pi(),num_emit*4)+pi()/2;
% 
% xc = r*cos(theta) + x_c;
% yc = r*sin(theta) + y_c;

xc = pts(1,:);
yc = pts(2,:);

for i = 1:length(xc)
    
    x = xc(i);
    y = yc(i) + r_c;
    
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
xc = xc(ind);
yc = yc(ind);
inds = inds(ind);

% find area element

for i = 1:length(theta)
    
    if i == 1
        ind1 = length(xc);
        ind2 = i + 1;
    elseif i == length(xc)
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
    
    dA(i) = (temp2-temp1)*r_c;
    
end

end

