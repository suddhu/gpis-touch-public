gamma = linspace(0,2*pi);
circle_shape = [0.03, 0.03];
circle_shape = [circle_shape(1)*cos(gamma'), circle_shape(2)*sin(gamma')];
save('circle.mat', 'circle_shape');
plot(circle_shape(:,1), circle_shape(:,2),'r-')
axis equal;