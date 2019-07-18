function [res] = loc_compare(loc1, loc2)
% Compare two loc's

    n = size(loc1, 1);
    m = size(loc2, 1);
    d = 0;
    
    for i = 1:n
        f = false;
        for j = 1:m
            
            if loc1(i, :) == loc2(j, :)
                f = true;
            end
        end
        
        if ~f
            d = d + 1;
        end
    end
    
    res = d;
    fprintf('%d diff in %d and %d.\n', d, n, m);
end
