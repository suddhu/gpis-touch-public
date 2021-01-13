gamma = linspace(0,2*pi);
ellip3_shape = [0.105/2, 0.157/2];
ellip3_shape = [ellip3_shape(1)*cos(gamma'), ellip3_shape(2)*sin(gamma')];
save('ellip3.mat', 'ellip3_shape');
plot(ellip3_shape(:,1), ellip3_shape(:,2),'r-')
axis equal;