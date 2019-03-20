%% level lines for l1 and transformed l1 \rho_a
clear 
close all
clc
%% l1-norm
x = -1:.1:1;
[X,Y] = meshgrid(x);
Z1 = abs(X) + abs(Y);
figure
contour(X,Y,Z1), axis([-1 1 -1 1]);
%title('l1 norm') % colorbar


%% transformed l1
a = 1;
x = -1:.1:1;
[X,Y] = meshgrid(x);
Z1 = (a+1)*abs(X)./(a+abs(X)) + (a+1)*abs(Y)./(a+abs(Y)) ;
figure
contour(X,Y,Z1), axis([-1 1 -1 1]);
%title('transformed l1 with a = 1') % colorbar


a = 100;
x = -1:.1:1;
[X,Y] = meshgrid(x);
Z1 = (a+1)*abs(X)./(a+abs(X)) + (a+1)*abs(Y)./(a+abs(Y)) ;
figure
contour(X,Y,Z1), axis([-1 1 -1 1]);
%title('transformed l1 with a = 100') % colorbar



a = 0.01;
x = -1:.1:1;
[X,Y] = meshgrid(x);
Z1 = (a+1)*abs(X)./(a+abs(X)) + (a+1)*abs(Y)./(a+abs(Y)) ;
figure
contour(X,Y,Z1), axis([-1 1 -1 1]);
%title('transformed l1 with a = 0.01') % colorbar


% a = 0.1;
% x = -1:.1:1;
% [X,Y] = meshgrid(x);
% Z1 = (a+1)*abs(X)./(a+abs(X)) + (a+1)*abs(Y)./(a+abs(Y)) ;
% figure
% contour(X,Y,Z1), axis([-1 1 -1 1]);
% title('transformed l1 with a = 0.1') % colorbar


