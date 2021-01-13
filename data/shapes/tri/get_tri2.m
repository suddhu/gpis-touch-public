a = 0.12587;
b = 0.15100;
c = 0.1962;
d = 0.090 / 2.0;  % from the rectangle coordinate system
tri2_shape = [[ d, d]; [d-b,d]; [d,d-a]];
save('tri2.mat', 'tri2_shape');
plot(tri2_shape(:,1), tri2_shape(:,2),'r-')
axis equal;