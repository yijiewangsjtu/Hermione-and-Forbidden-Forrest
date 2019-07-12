% Sensor error test

pdf1 = @(x) 0.3 .* (x>-5) .* (x<5) + 10 .* (x>-0.15) .* (x<0.15) ...
        + 7 .* (x>-0.7) .* (x<-0.3) + 7 .* (x>0.3) .* (x<0.7) ...
        + 3 .* (x>-1.2) .* (x<-0.8) + 3 .* (x>0.8) .* (x<1.2);
pdf2 = @(x) normpdf(x, 0, 1) * 10;
pdf3 = @(x) (2-x) .* (x>=0) .* (x<2) + (x+2) .* (x<0) .* (x>-2) + 0.1 .* (x>-5) .* (x<5);
proppdf = @(x, y) normpdf(x, y, 1);
proprnd = @(x) x + rand * 2 - 1;

smpl1 = mhsample(1, 50000, 'pdf', pdf1, 'proppdf', proppdf, 'proprnd', proprnd, 'burnin', 1000);
smpl2 = mhsample(1, 50000, 'pdf', pdf2, 'proppdf', proppdf, 'proprnd', proprnd, 'burnin', 1000);
smpl3 = mhsample(1, 50000, 'pdf', pdf3, 'proppdf', proppdf, 'proprnd', proprnd, 'burnin', 1000);

t = -5:0.5:5;
y1 = zeros(1, length(t));
y2 = zeros(1, length(t));
y3 = zeros(1, length(t));
for i = 1:length(t)

    y1(i) = nnz(smpl1 > t(i)-1 & smpl1 < t(i)+1);
    y2(i) = nnz(smpl2 > t(i)-1 & smpl2 < t(i)+1);
    y3(i) = nnz(smpl3 > t(i)-1 & smpl3 < t(i)+1);
end

figure('Name', '1m radius sensor error');
subplot(331)
plot(-5:0.05:5, pdf1(-5:0.05:5))
subplot(332)
plot(-5:0.05:5, pdf2(-5:0.05:5))
subplot(333)
plot(-5:0.04:5, pdf3(-5:0.04:5))
subplot(334)
hist(smpl1, 100);
subplot(335)
hist(smpl2, 100);
subplot(336)
hist(smpl3, 100);
subplot(337)
plot(t, y1);
subplot(338)
plot(t, y2);
subplot(339)
plot(t, y3);
