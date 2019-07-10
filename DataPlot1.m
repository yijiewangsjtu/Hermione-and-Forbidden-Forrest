% Plot

f1 = figure('Name', '2D trajectory');
hold on
C = [1 0 0; 0 1 0; 0 0 1];
for i = 1:135

    l = nnz(x_(:, i));  % zeros are invalid
    plot(x_(1:l, i), y_(1:l, i), 'Color', C(n_(i), :));
end

f2 = figure('Name', 'c1 height');
hold on
filter = c1__ > 0;
scatter3(x__(filter), y__(filter), c1__(filter), 8, c1__(filter));
colorbar
grid on

f3 = figure('Name', 'c3 height');
hold on
filter = c3__ > 0;
scatter3(x__(filter), y__(filter), c3__(filter), 8, c3__(filter));
colorbar
grid on
