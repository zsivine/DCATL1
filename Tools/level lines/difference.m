%% l1-norm, l2-norm
x = -10:.1:10;
[X,Y] = meshgrid(x);
Z1 = abs(X)+abs(Y)-sqrt(X.^2+Y.^2);
Z2 = (abs(X)+abs(Y))./sqrt(X.^2+Y.^2);
figure,
contour(X,Y,Z1), axis square
title('l1-l2') % colorbar
figure,
contour(X,Y,Z2), axis square
title('l1/l2')
% colorbar