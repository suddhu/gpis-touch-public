a = 0.08991/2.0;
b = 0.11258/2.0;
rect2_shape = [[a,b];[-a,b];[-a,-b];[a,-b];[a,b]];
save('rect2.mat', 'rect2_shape');
plot(rect2_shape(:,1), rect2_shape(:,2),'r-')
axis equal;
