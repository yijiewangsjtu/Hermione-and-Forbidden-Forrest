function [loc_, costRecord] = annealing(mloc, c)
% Simulated Annealing

    %% Initialization
    % Configurations of SA
    maxIterations = 400000;
    numCandPerIteration = 100;
    decTempCount = 300;
    sigma = 23;  % temperature, 
        % sigma = 1/alpha = 20 ... 0.05
    record = zeros(numCandPerIteration, 4);  % [index, x, y, expCost]
    decCount = 0;  % Only add 1 when a better move is found
    randCount = 0;
    stayCount = 0;
    cycle = 100;  % print and record cycle
    costRecord = zeros(maxIterations / cycle, 2);
    % Initialization for locations
    N = 2500;  % Total number of points
    l = 100;  % Length (both edge)
    loc_ = rand(N, 2) * l;
    C = [c, zeros(length(c), 1)];  % First col true observations, second col current obs
    % figure
    f = figure('Name', 'Annealing');
    set(f, 'position', [600, 1000, 800, 700]);
    f1 = subplot(2, 2, 1);
    title('True sensor distribution');
    filter = c > 0;
    scatter3(mloc(filter, 1), mloc(filter, 2), c(filter), 8, c(filter));
    colorbar
    view(2);
    f2 = subplot(2, 2, 2);
    f3 = subplot(2, 2, 3);
    f4 = subplot(2, 2, 4);

    %% Simulated Annealing
    [grid_mloc, grid_num] = SplitMloc(mloc, l);
    [cost_, C] = CalculateCost(loc_, mloc, grid_mloc, grid_num, C, l);
    UpdateFigure(loc_, mloc, C, f2, f3, f4, sigma);  % First plot
    for i = 1:maxIterations

        % print
        if mod(i, cycle) == 0
            
            fprintf('\nIter:%d  cost:%d  dec:%d/%d  fail:%d/%d', i, cost_, decCount, decTempCount, randCount, cycle);
            costRecord(i/cycle, 1) = cost_;
            costRecord(i/cycle, 2) = sigma;
            randCount = 0;
            stayCount = stayCount + 1;
            if stayCount > 7 && decCount/decTempCount < 0.33
            
                decCount = decTempCount;
                fprintf('\nSkip this sigma.');
            end
        end
        % Whether to cool down
        if (decCount == decTempCount)

            if (stayCount > 4)
                if sigma > 14.1, sigma = sigma - 0.35;end
                if sigma > 3.1 && sigma <=14.1, sigma = sigma - 1;end
                if sigma > 0.91 && sigma <=3.1, sigma = sigma - 0.3;end
                if sigma > 0.51 && sigma <=0.91
                    decTempCount = 50;
                    sigma = sigma - 0.01;
                end
                if sigma > 0.011 && sigma <=0.51, sigma = sigma - 0.01;end
                if sigma > 0.0011 && sigma <=0.011, sigma = sigma - 0.001;end
                if sigma <= 0.001, sigma = 0.0011;end
            end
            decCount = 0;
            stayCount = 0;
            fprintf('\n==== Current sigma: %f ====', sigma);
            UpdateFigure(loc_, mloc, C, f2, f3, f4, sigma);
        end

        flag = false;  % Whether find a better loc_
        for j = 1:numCandPerIteration

            % get a random move and re-evaluate
            k = randi(N);
            while (sigma < 10) && (loc_(k, 1) == 0 || loc_(k, 2)==0 || loc_(k, 1)==l || loc_(k, 2)==l)
            
                k = randi(N);
            end
            record(j, 1) = k;
            delta = normrnd(0, sigma, [1, 2]);  % random move dir and mag
            delta_x = delta(1);
            delta_y = delta(2);
            old_x = loc_(k, 1);
            old_y = loc_(k, 2);
            record(j, 2) = Update(old_x, delta_x, l);
            record(j, 3) = Update(old_y, delta_y, l);
            cost = UpdateCost([old_x, old_y], [record(j, 2), record(j, 3)], mloc, grid_mloc, grid_num, C, l);

            % whether this random move is better
            if (cost < cost_)

                loc_(k, 1) = record(j, 2);
                loc_(k, 2) = record(j, 3);
                [cost_, C] = UpdateCost([old_x, old_y], [loc_(k, 1), loc_(k, 2)], mloc, grid_mloc, grid_num, C, l);
                flag = true;
                decCount = decCount + 1;
                break;
            end

            % record cost
            expCost = exp((cost_ - cost) / sigma);
            record(j, 4) = expCost;
        end
        % If no better move found, randomly choose based on expCost
        if ~flag

            chosen = randsample(numCandPerIteration, 1, true, record(:, 4));
            [cost_, C] = UpdateCost([loc_(record(chosen, 1), 1), loc_(record(chosen, 1), 2)], [record(chosen, 2), record(chosen, 3)], mloc, grid_mloc, grid_num, C, l);
            loc_(record(chosen, 1), 1) = record(chosen, 2);
            loc_(record(chosen, 1), 2) = record(chosen, 3);
            randCount = randCount + 1;
        end
    end
    
    %% output
    
end


function [grid_mloc, grid_num] = SplitMloc(mloc, l)

    grid_mloc = zeros(l, l, 48);
    grid_num = zeros(l, l, 1);
    for i = 1:length(mloc)
    
        xx = ceil(mloc(i, 1));
        yy = ceil(mloc(i, 2));
        
        grid_num(xx, yy) = grid_num(xx, yy) + 1;
        grid_mloc(xx, yy, grid_num(xx, yy)) = i;
    end
end


function [cost, C] = CalculateCost(loc, mloc, grid_mloc, grid_num, C, l)
% Given n points located at 'loc' [n*2]
% Measure m locations located at 'mloc' [m*2]
% Number of points in a 1m radius 'C' [m*2] (both true and observed)

    n = size(loc, 1);

    for i = 1:n

        xx = ceil(loc(i, 1));
        yy = ceil(loc(i, 2));
        for xxi = xx-1:xx+1
            for yyj = yy-1:yy+1
                
                if xxi<1 || yyj<1 || xxi>l || yyj>l

                    continue;
                end
                
                for j = 1:grid_num(xxi, yyj)

                    index = grid_mloc(xxi, yyj, j);
                    if norm([loc(i, 1)-mloc(index, 1), loc(i, 2)-mloc(index, 2)]) < 1

                       C(index, 2) = C(index, 2) + 1; 
                    end
                end
            end
        end
    end
    
    cost = sum(abs(C(:, 1) - C(:, 2)));
end


function [cost, C] = UpdateCost(old_loc, new_loc, mloc, grid_mloc, grid_num, C, l)

    % For old loc
    xx = ceil(old_loc(1));
    yy = ceil(old_loc(2));
    for xxi = xx-1:xx+1
        for yyj = yy-1:yy+1

            if xxi<1 || yyj<1 || xxi>l || yyj>l

                continue;
            end

            for j = 1:grid_num(xxi, yyj)

                index = grid_mloc(xxi, yyj, j);
                if norm([old_loc(1)-mloc(index, 1), old_loc(2)-mloc(index, 2)]) < 1

                   C(index, 2) = C(index, 2) - 1; 
                end
            end
        end
    end
    
    % For new loc
    xx = ceil(new_loc(1));
    yy = ceil(new_loc(2));
    for xxi = xx-1:xx+1
        for yyj = yy-1:yy+1

            if xxi<1 || yyj<1 || xxi>l || yyj>l

                continue;
            end

            for j = 1:grid_num(xxi, yyj)

                index = grid_mloc(xxi, yyj, j);
                if norm([new_loc(1)-mloc(index, 1), new_loc(2)-mloc(index, 2)]) < 1

                   C(index, 2) = C(index, 2) + 1; 
                end
            end
        end
    end

    cost = sum(abs(C(:, 1) - C(:, 2)));
end


function [res] = Update(x, delta, l)

    res = x + delta;
    if res > l, res = l;end
    if res < 0, res = 0;end
end


function [] = UpdateFigure(loc, mloc, C, f2, f3, f4, sigma)

    cla(f2);
    subplot(2, 2, 2);
    filter = C(:, 2) > 0;
    scatter3(mloc(filter, 1), mloc(filter, 2), C(filter, 2), 8, C(filter, 2));
    title(['Currunt sigma: ', num2str(sigma)]);
    colorbar
    xlim([0 100])
    ylim([0 90])
    view(2);
    
    cla(f3);
    subplot(2, 2, 3);
    scatter(loc(:, 1), loc(:, 2), 8, 'filled');
    xlim([0 100])
    ylim([0 90])
    
    cla(f4);
    subplot(2, 2, 4);
    c_ = C(:, 1) - C(:, 2);
    filter = (c_ > 0) | (c_ < 0);
    scatter3(mloc(filter, 1), mloc(filter, 2), c_(filter), 8, c_(filter));
    colorbar
    xlim([0 100])
    ylim([0 90])
    view(2);
    title('Difference |delta|>0');
    
    shg
end
