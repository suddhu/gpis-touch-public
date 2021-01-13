function [data] = offset_and_scale_data( data, inputscale, inputmean )
    %% re-center and re-scale
    data(:,1) = data(:,1)*inputscale + repmat(inputmean(1),size(data(:,1),1),1);
    data(:,2) = data(:,2)*inputscale + repmat(inputmean(2),size(data(:,2),1),1);
end

