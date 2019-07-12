% Plot

%f0 = figure('Name', '2D dots');
subplot(2, 3, 1)
hold on
C = [1 0.1 0.1; 0.1 1 0.1; 0.1 0.1 1];
for i = 1:135

    l = nnz(x_(:, i));  % zeros are invalid
    scatter(x_(1:l, i), y_(1:l, i), 4, C(n_(i), :));
end
grid on

%f1 = figure('Name', '2D trajectory');
subplot(2, 3, 4)
hold on
C = [1 0.1 0.1; 0.1 1 0.1; 0.1 0.1 1];
for i = 1:135

    l = nnz(x_(:, i));  % zeros are invalid
    plot(x_(1:l, i), y_(1:l, i), 'Color', C(n_(i), :));
end
grid on

% sample distribution
count = zeros(100, 90);
for i = 1:100
    for j = 1:90
    
        count(i, j) = nnz((i-1 <= x__) & (x__ < i) & (j-1 <= y__) & (y__ < j));
    end
end
[yy, xx] = meshgrid(1:90, 1:100);
filter = count > 0;
xx = xx(:);
yy = yy(:);
subplot(2, 3, 2)
hold on
scatter3(xx(filter), yy(filter), count(filter), 7, count(filter));
colorbar

%f2 = figure('Name', 'c1 height');
subplot(2, 3, 3)
hold on
filter = c1__ > 0;
scatter3(x__(filter), y__(filter), c1__(filter), 8, c1__(filter));
colorbar
grid on

%f3 = figure('Name', 'c3 height');
subplot(2, 3, 6)
hold on
filter = c3__ > 0;
scatter3(x__(filter), y__(filter), c3__(filter), 8, c3__(filter));
colorbar
grid on


