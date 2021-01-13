a = 0.12561;
b = 0.1765;
c = 0.2152;
d = 0.090 / 2.0;  % from the rectangle coordinate system
tri3_shape = [[d, d]; [d-b,d]; [d,d-a]];
save('tri3.mat', 'tri3_shape');
plot(tri3_shape(:,1), tri3_shape(:,2),'r-')
axis equal;