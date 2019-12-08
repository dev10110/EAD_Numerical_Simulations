function [xe,ye,theta] = getMOCStart( num_emit, d, r_e)

offset = 0;

x_e = 0;
y_e = d + r_e;

space = pi()/4/(num_emit);
%space = pi()/(num_emit);

r = r_e*(1 + offset);
theta1 = linspace(space,pi()/4,floor(1.5*num_emit)+1)+pi()/2;
theta2 = linspace(pi()/4+space,pi(),num_emit)+pi()/2;
%theta2 = linspace(pi()/2+space,2*pi()-space,num_emit)+pi()/6;
theta = [theta1, theta2];
%theta = theta2;

xe = r*cos(theta) + x_e;
ye = r*sin(theta) + y_e;

end

