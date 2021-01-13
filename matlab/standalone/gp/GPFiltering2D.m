%% Input:
% data: N x 4;
% varNoise: measurement noise;
% R: kernel parameter
% shape_data: ground truth data 

%% Output:
% data_output: N x 4;
function GPFiltering2D(X, Y, varNoise, R, shape_data, Scale, Mean)
    evaluateAndPlot(X, Y, Scale, Mean, R, varNoise, shape_data, 1.0, 0.05);
end
