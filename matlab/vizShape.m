% viz_shape.m: Visualize contact measurement results
%
% Sudharshan Suresh <suddhu@cmu.edu>, 2020
%
% Copyright: This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License v3 as published by
% the Free Software Foundation. This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of any FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License v3 for more details. You should have received a copy of the GNU General Public License v3
% along with this program; if not, you can access it online at http://www.gnu.org/licenses/gpl-3.0.html.

figure(1); 
%% plot 1: contour 
s1 = subplot(2, 1,1); 
cla(s1);

hold on; axis equal; grid off;

if (~isempty(X_active))
    plot(X_active(:,1), X_active(:,2), 'bo', 'MarkerSize',5, 'LineWidth',2, 'MarkerFaceColor', 'b'); % contact points
    quiver(X_active(:,1), X_active(:,2), Y_active(:,1),Y_active(:,2),0.3,'LineWidth',1.5, 'color', 'b'); % normal directions
end

xlim([-0.09, 0.09])
ylim([-0.09, 0.09])

patch(C(:,1), C(:,2), 'green','FaceAlpha',.3, 'LineWidth',3,  'EdgeColor','green'); % object shape
plot(shape_data(:,1), shape_data(:,2), 'k-.'); % ground truth

title('0 Level-set contour');    
xlabel('x', 'FontSize',12,'FontWeight', 'bold','Color','r'); ylabel('y','FontSize',12,'FontWeight',  'bold','Color','r'); 
legend('contact measurements', 'object contour','ground truth', 'Location','SouthEast');
hold off; 

%% plot 2: GP
subplot(2,1,2);
x_range = -testLim:testRes:testLim; 
y_range = -testLim:testRes:testLim; 
[xt, yt] = meshgrid(x_range,y_range); % scale meshgrid 
s = surf(xt,yt,fMean', 'FaceColor', 'interp'); % plot GP (mex output flips matrix?)
hold on; grid on;  
xlim([-0.09, 0.09])
ylim([-0.09, 0.09])
colormap jet; 
cb = colorbar;
ylabel(cb, 'variance')
s.CData = fVar';
caxis([0, 5*varNoise(1)]);
set(gca, 'ZDir','reverse');

% plot 0 plane
Z = zeros(size(xt));
CO(:,:,1) = zeros(size(xt,1));
CO(:,:,2) = zeros(size(xt,1));
CO(:,:,3) = zeros(size(xt,1)); 
s = surf(xt,yt,Z, CO); alpha(s,.5)
plot(C(:, 1), C(:, 2), 'g-', 'linewidth', 10); % object contour
% fMean
zlim([min(fMean(:)) 0.1])
%view(2);

legend('GP-ISP', 'ISP = 0 plane', 'object contour', 'Location','SouthEast');
title('GP Implicit surface potential'); 
xlabel('x','FontSize',12,'FontWeight',  'bold','Color','r'); ylabel('y','FontSize',12,'FontWeight',  'bold','Color','r'); zlabel('ISP','FontSize',12,'FontWeight',  'bold','Color','r'); 
hold off;    

