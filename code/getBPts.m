function [ pts_e, pts_c, e_inds, c_inds ] = getBPts( p, e, t, d, r_e, r_c )

% find triangles on boundaries
pinds = e(1,:);
count = 1;

for i = 1:length(t)
    count_pts = 0;
    for j = 1:3
        for k = 1:length(pinds)
            if t(j,i) == pinds(k)
                count_pts = count_pts + 1;
                break
            end
        end
    end
    if count_pts == 2
        t_store(count) = i;
        count = count + 1;
    end
end

% eliminate duplicates
t_store = unique(t_store);

% Triangle point indices
it1=t(1,t_store);
it2=t(2,t_store);
it3=t(3,t_store);

% Find centroids of triangles
xpts=(p(1,it1)+p(1,it2)+p(1,it3))/3;
ypts=(p(2,it1)+p(2,it2)+p(2,it3))/3;

de = sqrt(xpts.^2 + (ypts - (d + r_e)).^2);
dc = sqrt(xpts.^2 + (ypts - (-r_c)).^2);

count_e = 1;
count_c = 1;

for i = 1:length(de)
    
    if de(i) < 2*r_e + 0.001
        
        pts_e(1,count_e) = xpts(i);
        pts_e(2,count_e) = ypts(i); 
        e_inds(count_e) = t_store(i);
        count_e = count_e + 1;
        
    elseif dc(i) < 2*r_c
        
        pts_c(1,count_c) = xpts(i);
        pts_c(2,count_c) = ypts(i);   
        c_inds(count_c) = t_store(i);
        count_c = count_c + 1;        
        
    end
    
end

end

