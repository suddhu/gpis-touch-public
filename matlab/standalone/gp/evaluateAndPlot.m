%% Input:
% X: N x 2 (scaled 1 to -1)
% Y: N x 3
% Scale: data scale
% Mean: data mean
% R: kernel parameter
% KernelTest: 3 x 3
% ground_truth: M x 2
% mesh_lim: max/min 2D test data 
function evaluateAndPlot(X, Y, Scale, Mean, R, v, ground_truth, mesh_lim, mesh_res)
    
    tic; 
    %% generate 2D mesh
    x_range = -mesh_lim:mesh_res:mesh_lim; 
    y_range = -mesh_lim:mesh_res:mesh_lim; 
    [xt, yt] = meshgrid(x_range,y_range); 
    
    test_sz = size(xt);
    x_test = [xt(:),yt(:)];
    
    % scale items
    R = R/Scale; 
    scaled_v = v./Scale^2;

    KernelMat = KernelFun2D(X,X,R) + diag(repmat(scaled_v,1,size(X,1))); % (K(X, X) + sigma^2I)
    fMean = zeros(size(x_test,1),1);
    fVar = zeros(size(x_test,1),1);
    
    % check for PSD (The R used should assure PSD)
    try chol(KernelMat);
%         disp('Matrix is symmetric positive definite.')
    catch ME
        disp('Matrix is not symmetric positive definite')
    end

    %% plot kernel 
%     plotKernel(KernelMat,X);

    %% query GP with test data
    for i=1:size(x_test,1)
        KernelTest = KernelFun2D(x_test(i,1:2), x_test(i,1:2), R); % K(x*, x*)
        KernelTestTrain = KernelFun2D(X, x_test(i,1:2), R); % K*
        
        %% get mean
        if(norm(x_test(i,1:2)) > 1.5)
            fMean(i) = 1;             % outside the object 
        else
            offset = [1,0,0]; 
            FunMean = offset' + KernelTestTrain'*(KernelMat\(Y - repmat(offset',size(Y,1)/3,1)));
            fMean(i) = FunMean(1); % d
        end
        
        %% get variance
        tmp = [1 0 0];
        fVar(i) = tmp*KernelTest*tmp' - (tmp*KernelTestTrain')*(KernelMat\(KernelTestTrain*tmp'));
        fVar(i) = (Scale^2)*fVar(i); % scale back 
    end
    
    disp(['# data points: ', num2str(size(X,1)), ' time: ', num2str(toc)]);
    %% scale and recenter test data
    xt = xt(:)*Scale + repmat(Mean(1),size(xt(:),1),1);
    yt = yt(:)*Scale + repmat(Mean(2),size(yt(:),1),1);
    xt = reshape(xt,test_sz);
    yt = reshape(yt,test_sz);

    figure(1); 
    %% plot 1: contour 
    subplot(2, 1,1); 
    X = offset_and_scale_data( X, Scale, Mean );
    X_on = X(1:end,1:2); % training points input
    Y_on = [Y(2:3:end-1), Y(3:3:end)]; % training points output

    plot(X_on(:,1), X_on(:,2), 'bo', 'MarkerSize',5, 'LineWidth',2, 'MarkerFaceColor', 'b'); % contact points
    hold on; axis equal;
    xlim([-0.08, 0.08]);
    ylim([-0.08, 0.08]);
    
    fMean = reshape(fMean,test_sz); 
    fVar = reshape(fVar,test_sz);    
        
    %% scale mesh to original size
    x_range = x_range*Scale + Mean(1); 
    y_range = y_range*Scale + Mean(2); 
    C = contourc(x_range,y_range,fMean, [0 0]); 
    
    % delete weird outliers
    thresh = 1;
    idx = any(C > thresh,1); 
    m_toDelete = C(:,idx);
    C(:,idx) = [];
    
    patch(C(1,2:end), C(2,2:end), 'green','FaceAlpha',.3, 'LineWidth',3,  'EdgeColor','green'); % shape patch
    plot(ground_truth(:,1), ground_truth(:,2), 'k-.'); % ground truth
    quiver(X_on(:,1), X_on(:,2), Y_on(:,1),Y_on(:,2),0.3,'LineWidth',1.5, 'color', 'b'); % normal directions
    
    title('0 Level-set contour');    
    xlabel('x', 'FontSize',12,'FontWeight', 'bold','Color','r'); ylabel('y','FontSize',12,'FontWeight',  'bold','Color','r'); 
    legend('contact measurements', 'object contour','ground truth', 'Location','SouthEast');
    hold off; 

    %% plot 2: GP
    subplot(2,1,2);
   
    s = surf(xt,yt,fMean, 'FaceColor', 'interp'); % plot GP
    hold on; grid on;  
    xlim([-0.08, 0.08])
    ylim([-0.08, 0.08])
    colormap jet; 
    cb = colorbar;
    ylabel(cb, 'variance')
    s.CData = fVar;
    caxis([0,5*v(1)]); % 5*v(1)
    set(gca, 'ZDir','reverse');
    
    % plot 0 plane
    Z = zeros(size(xt));
    CO(:,:,1) = zeros(size(xt,1));
    CO(:,:,2) = zeros(size(xt,1));
    CO(:,:,3) = zeros(size(xt,1)); 
    s = surf(xt,yt,Z, CO); alpha(s,.7)
    plot(C(1, 2:end), C(2, 2:end), 'g-', 'linewidth', 10); % object contour
    zlim([min(fMean(:)) 0.5]);
%     view(2);

    legend('GP-ISP', 'ISP = 0 plane', 'object contour', 'Location','SouthEast');
    title('GP Implicit surface potential'); 
    xlabel('x','FontSize',12,'FontWeight',  'bold','Color','r'); ylabel('y','FontSize',12,'FontWeight',  'bold','Color','r'); zlabel('ISP','FontSize',12,'FontWeight',  'bold','Color','r'); 
    hold off;    
end
