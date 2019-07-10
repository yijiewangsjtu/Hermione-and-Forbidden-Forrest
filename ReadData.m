% Read in data

D = csvread('data_proj_414.csv', 1, 1);

x__ = D(:, 1);
y__ = D(:, 2);
n__ = D(:, 3) + D(:, 4) * 2 + D(:, 5) * 3;
t__ = D(:, 6);
c1__ = D(:, 7);
c3__ = D(:, 8);

newt = t__ + 50 * n__;
g = findgroups(newt);
x_ = zeros(240, 135);
y_ = zeros(240, 135);
n_ = zeros(1, 135);
c1_ = zeros(240, 135);
c3_ = zeros(240, 135);
for i = 1:135
    
    boolArray = (g == i);
    l = sum(boolArray);
    x_(1:l, i) = x__(boolArray);
    y_(1:l, i) = y__(boolArray);
    n_(i) = floor(mean(newt(boolArray)) / 50);
    c1_(1:l, i) = c1__(boolArray);
    c3_(1:l, i) = c3__(boolArray);
end
