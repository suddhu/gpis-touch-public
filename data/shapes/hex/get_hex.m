side = 0.06050;
n = 6;
hex_shape = [];
for i = 1:n
    theta = (2*pi/n)*i;
    hex_shape = [ hex_shape; [side*cos(theta), side*sin(theta)] ];
end

save('hex.mat', 'hex_shape');
plot(hex_shape(:,1), hex_shape(:,2),'r-')
axis equal;


    