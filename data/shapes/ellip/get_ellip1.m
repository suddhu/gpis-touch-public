gamma = linspace(0,2*pi);
ellip1_shape = [0.105/2, 0.105/2];
ellip1_shape = [ellip1_shape(1)*cos(gamma'), ellip1_shape(2)*sin(gamma')];
save('ellip1.mat', 'ellip1_shape');
plot(ellip1_shape(:,1), ellip1_shape(:,2),'r-')
axis equal;