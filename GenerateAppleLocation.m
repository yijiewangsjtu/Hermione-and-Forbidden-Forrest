% Generate Apple Location
function [res] = GenerateAppleLocation(loc, x, y)

    n = length(loc);
    m = length(x);
    
    count = 0;
    res = zeros(100, 2);
    for i = 1:n
        for j = 1:m
        
            appleLoc = [loc(i, 1), loc(i, 2)];
            sensorLoc = [x(j), y(j)];
            if norm(appleLoc - sensorLoc) < 1
            
                count = count + 1;
                res(count, :) = loc(i, :);
                break;
            end
        end
    end
end
