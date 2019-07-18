function [loc_] = auto_run(x__, y__, c1__, loc_)
% temp

    for i = 1:10
        
        loc_ = annealing([x__, y__], c1__, loc_);
    end
end
