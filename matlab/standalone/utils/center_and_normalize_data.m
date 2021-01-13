function [output_data,inputmean,inputscale] = center_and_normalize_data( dynamics_data )
    %CENTER_AND_NORMALIZE_DATA 
    % output_data: N*D
    dynamics_data=dynamics_data';

    N = size(dynamics_data,1)/2;
    M = size(dynamics_data,2);
    output_data = dynamics_data;
    
    %% compute input mean and apply
%     inputmean = [0; 0];
    inputmean = mean(output_data(1:N,:),2);
    output_data(1:N,:) = output_data(1:N,:) - repmat(inputmean,1,M); % center 

    %% compute input scale and apply
    mx = max(max(output_data(1:N,:)'));
    mn = min(min(output_data(1:N,:)'));
    inputscale=max(abs(mx), abs(mn));
%     inputscale = 1;
    output_data(1:N,:) = output_data(1:N,:)./inputscale; % scale b/w -1 to 1
    
    output_data=output_data';
    inputmean=inputmean';
end

