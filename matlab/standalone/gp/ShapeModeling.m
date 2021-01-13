%% shape modeling of 2D objects from contact and normal measurements

clc; clear; 

% data: x, y, nx, ny
addpath('../utils/')  

%% load shape data 
shape = 'rect1';
id = '20200115-1026';
% id = '20200131-1608';

% shape = 'ellip2';
% id = '20200115-2208';
% id = '20200213-1445';

% shape = 'butter';
% id = '20200115-2213';

file = strcat('../../../data/contacts/contacts-', shape,'-', id, '.txt');
data = importdata(file);

data(:,3:4) = normr(data(:,3:4)); % unit normals

if strcmp(shape,'rect1')
    shape_data = struct2array(load('../../../data/shapes/rect/rect1.mat','rect1_shape'));
elseif strcmp(shape,'ellip2')
    shape_data = struct2array(load('../../../data/shapes/ellip/ellip2.mat','ellip2_shape'));
elseif strcmp(shape,'butter')
    shape_data = struct2array(load('../../../data/shapes/butter/butter.mat','butter_shape'));
end

%% plot contacts and normals
% plotContactsAndNormals(data,shape_data); 

varNoise = [1e-4, 1e-4, 1e-4]; % measurement noise

for i = 2:size(data,1)
    R = 1.1*max(pdist(data(1:i, :),'euclidean')); % kernel parameter (max dist in the set of points)
    [new_data, inputmean, inputscale] = center_and_normalize_data(data(1:i, :)); % preprocess data
    curr_X = new_data(:,1:2);
    curr_Y = reshape([zeros(size(new_data(:,3:4),1), 1)' ; new_data(:,3:4)'], [], 1);
    GPFiltering2D(curr_X, curr_Y, varNoise, R, shape_data, inputscale, inputmean);   
end

