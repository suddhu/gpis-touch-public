a = 0.12587;
b = 0.12590;
c = 0.178;
d = 0.090 / 2.0;  % from the rectangle coordinate system
tri1_shape = [[d, d]; [d-b,d]; [d,d-a]];
save('tri1.mat', 'tri1_shape');
plot(tri1_shape(:,1), tri1_shape(:,2),'r-')
axis equal;

      