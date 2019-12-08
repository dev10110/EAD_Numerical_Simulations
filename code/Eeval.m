function [ prime ] = Eeval(~, x, Ex, Ey )
%Eeval calculates the electric field strength at a give point x
%   t - current time
%   x - position [x; y], multiple columns allowed
%   Ex - interpolating function for Ex in form f(x,y)
%   Ey - interpolating function for Ey in form f(x,y)


prime = [Ex(x(1,:),x(2,:));...
         Ey(x(1,:),x(2,:))];

end

