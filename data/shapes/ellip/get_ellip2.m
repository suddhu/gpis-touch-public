gamma = linspace(0,2*pi);
ellip2_shape = [0.105/2, 0.13089/2];
ellip2_shape = [ellip2_shape(1)*cos(gamma'), ellip2_shape(2)*sin(gamma')];
save('ellip2.mat', 'ellip2_shape');
plot(ellip2_shape(:,1), ellip2_shape(:,2),'r-')
axis equal;