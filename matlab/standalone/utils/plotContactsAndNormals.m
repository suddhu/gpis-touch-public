function plotContactsAndNormals(data,gt)
    plot(data(:,1),data(:,2),'o','color',[0 0.392157 0],'MarkerSize',8, 'MarkerFaceColor',[0 0.392157 0]); hold on;
    axis equal;
    quiver(data(:,1),data(:,2),data(:,3),data(:,4), 'color',[0 0.392157 0]);
    plot(gt(:,1), gt(:,2), 'k-.');
    title('Contacts and normals');
    xlabel('x', 'FontSize',12,'FontWeight', 'bold','Color','r'); ylabel('y','FontSize',12,'FontWeight',  'bold','Color','r'); 
    hold off; 
end

