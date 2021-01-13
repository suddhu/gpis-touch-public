% GPShape.m: shape modeling of 2D objects from contact and normal measurements
%
% Sudharshan Suresh <suddhu@cmu.edu>, 2020
%
% Copyright: This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License v3 as published by
% the Free Software Foundation. This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of any FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License v3 for more details. You should have received a copy of the GNU General Public License v3
% along with this program; if not, you can access it online at http://www.gnu.org/licenses/gpl-3.0.html.

clc; clear;

% data: x, y, nx, ny
addpath('../mex');


%% load shape data (stored in the data/contacts folder)
shape = 'rect1';
id = '20200115-1026';
% id = '20200131-1608';
% id = '20200419-0028';

% shape = 'ellip2';
% id = '20200115-2208';
% id = '20200213-1445';

% shape = 'butter';
% id = '20200115-2213';

file = strcat('../data/contacts/contacts-', shape,'-', id, '.txt');
data = importdata(file);

data(:,3:4) = normr(data(:,3:4)); % unit normals

%% ground truth shape data
if strcmp(shape,'rect1')
    shape_data = struct2array(load('../data/shapes/rect/rect1.mat','rect1_shape'));
elseif strcmp(shape,'ellip2')
    shape_data = struct2array(load('../data/shapes/ellip/ellip2.mat','ellip2_shape'));
elseif strcmp(shape,'butter')
    shape_data = struct2array(load('../data/shapes/butter/butter.mat','butter_shape'));
end

%% plot contacts and normals (sanity check)
% plot(data(:,1),data(:,2),'o','color',[0 0.392157 0],'MarkerSize',8, 'MarkerFaceColor',[0 0.392157 0]); hold on;
% axis equal;
% quiver(data(:,1),data(:,2),data(:,3),data(:,4), 'color',[0 0.392157 0]);
% plot(shape_data(:,1), shape_data(:,2), 'k-.');
% title(['Contacts and normals']);
% xlabel('x', 'FontSize',12,'FontWeight', 'bold','Color','r'); ylabel('y','FontSize',12,'FontWeight',  'bold','Color','r'); 
% hold off; 

X_all = data(:,1:2);
Y_all = data(:,3:4);
X_active = []; Y_active = []; 

%% Hyperparameters
varNoise = [1e-5, 1e-5, 1e-5]; % measurement noise
priorNoise = [1e-3, 1e-3, 1e-3]; % circular prior noise
testLim = 1.5; testRes =  0.01; localGPs = 1; priorRad = 0.03; % GP parameters

mexGPShape('init',varNoise, priorNoise, testLim, testRes, localGPs, priorRad);  % init GP from cpp
[C, fMean, fVar] = mexGPShape('test'); % initialize with prior
vizShape % plot prior

for i = 1:size(X_all,1)
    X_curr = X_all(i, :); 
    Y_curr = Y_all(i, :);
    mexGPShape('update', X_curr, Y_curr); % update GP
    [C, fMean, fVar] = mexGPShape('test'); % test 2D data
    X_active = [X_active; X_curr]; 
    Y_active = [Y_active; Y_curr]; 
    vizShape % plot results
end

mexGPShape('reset');  

