function plotKernel(KernelMat,X)
    figure(2);
    hold off;
    imagesc(KernelMat); 
    hold on;
    colormap(jet);
    axis equal;
    xlim([1 size(KernelMat,1)]); ylim([1 size(KernelMat,2)]); 
    xlabel('Xi', 'FontSize',12,'FontWeight', 'bold','Color','r'); ylabel('Xj', 'FontSize',12,'FontWeight', 'bold','Color','r');
    title(['Kernel Matrix with ', num2str(size(X,1)), ' datapoints']);     
end

