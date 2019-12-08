function [ flow_switch ] = checkFlow( U, d )

% check to see if interpolation function exists already for potential flow
% configuration

list = dir;
file = ['flow_interp_U',num2str(U),'_d',num2str(100*d),'.mat'];

flow_switch = 0;

for i = 1:length(list)
    if strcmp(file,list(i).name)
        flow_switch = 1;
    end
end

end

