a = 0.09/2.0;
b = 0.09/2.0;
rect1_shape = [[a,b];[-a,b];[-a,-b];[a,-b];[a,b]];
save('rect1.mat', 'rect1_shape');
plot(rect1_shape(:,1), rect1_shape(:,2),'r-')
axis equal;