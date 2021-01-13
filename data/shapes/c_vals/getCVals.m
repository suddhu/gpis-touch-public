% compute the eals and get c value by integration 
addpath('./tools');
addpath('./Mesh2dv24');

rc = [0, 0]';
circle = load("/home/suddhu/software/GPIS/data/shapes/circle/circle.mat");
eals = EllipsoidApproxLimitSurface(circle.circle_shape', rc);
fprintf('circle.c: %d.\n',eals.c);

rect1 = load("/home/suddhu/software/GPIS/data/shapes/rect/rect1.mat");
eals = EllipsoidApproxLimitSurface(rect1.rect1_shape', rc);
fprintf('rect1.c: %d.\n',eals.c);

rect2 = load("/home/suddhu/software/GPIS/data/shapes/rect/rect2.mat");
eals = EllipsoidApproxLimitSurface(rect2.rect2_shape', rc);
fprintf('rect2.c: %d.\n',eals.c);

rect3 = load("/home/suddhu/software/GPIS/data/shapes/rect/rect3.mat");
eals = EllipsoidApproxLimitSurface(rect3.rect3_shape', rc);
fprintf('rect3.c: %d.\n',eals.c);

ellip1 = load("/home/suddhu/software/GPIS/data/shapes/ellip/ellip1.mat");
eals = EllipsoidApproxLimitSurface(ellip1.ellip1_shape', rc);
fprintf('ellip1.c: %d.\n',eals.c);

ellip2 = load("/home/suddhu/software/GPIS/data/shapes/ellip/ellip2.mat");
eals = EllipsoidApproxLimitSurface(ellip2.ellip2_shape', rc);
fprintf('ellip2.c: %d.\n',eals.c);

ellip3 = load("/home/suddhu/software/GPIS/data/shapes/ellip/ellip3.mat");
eals = EllipsoidApproxLimitSurface(ellip3.ellip3_shape', rc);
fprintf('ellip3.c: %d.\n',eals.c);

hex = load("/home/suddhu/software/GPIS/data/shapes/hex/hex.mat");
eals = EllipsoidApproxLimitSurface(hex.hex_shape', rc);
fprintf('hex.c: %d.\n',eals.c);

rc = [-0.00541556, 0.00309479]';
tri1 = load("/home/suddhu/software/GPIS/data/shapes/tri/tri1.mat");
eals = EllipsoidApproxLimitSurface(tri1.tri1_shape', rc);
fprintf('tri1.c: %d.\n',eals.c);

rc = [0.00303333, 0.00306333]';
tri2 = load("/home/suddhu/software/GPIS/data/shapes/tri/tri2.mat");
eals = EllipsoidApproxLimitSurface(tri2.tri2_shape', rc);
fprintf('tri2.c: %d.\n',eals.c);

rc = [-0.01387576, 0.0031902]';
tri3 = load("/home/suddhu/software/GPIS/data/shapes/tri/tri3.mat");
eals = EllipsoidApproxLimitSurface(tri3.tri3_shape', rc);
fprintf('tri3.c: %d.\n',eals.c);

% butter = load("/home/suddhu/software/GPIS/data/shapes/butter/butter.mat");
% eals = EllipsoidApproxLimitSurface(butter.butter_shape');
% fprintf('butter.c: %d.\n',eals.c);
