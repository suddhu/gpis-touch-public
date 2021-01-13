a = 0.09895/2.0;
b = 0.09895/2.0;
steel_shape = [[a,b];[-a,b];[-a,-b];[a,-b];[a,b]];
save('steel.mat', 'steel_shape');
plot(steel_shape(:,1), steel_shape(:,2),'r-')
axis equal;