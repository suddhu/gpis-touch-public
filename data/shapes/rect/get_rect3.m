a = 0.13501/2.0;
b = 0.08994/2.0;
rect3_shape = [[a,b];[-a,b];[-a,-b];[a,-b];[a,b]];
save('rect3.mat', 'rect3_shape');
plot(rect3_shape(:,1), rect3_shape(:,2),'r-')
axis equal;
