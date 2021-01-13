
% M [ x1 x2 x3 ...
%     y1 y2 y3 ...]
% test if a polygon is convex

function ret = isConvex(M)
    N = size(M,2);
    if (N<4)
        ret = 1;
        return;
    end
        
    x0 = M(1, 1:end);
    x1 = [x0(1,2:end), x0(1,1)];
    x2 = [x0(1,3:end), x0(1,1:2)];
    y0 = M(2, 1:end);
    y1 = [y0(1,2:end), y0(1,1)];
    y2 = [y0(1,3:end), y0(1,1:2)];
    dx1 = x2 - x1;
    dy1 = y2 - y1;
    dx2 = x0 - x1;
    dy2 = y0 - y1;
    zcrossproduct = dx1 .* dy2 - dy1 .* dx2;
    
    t1 = sum(zcrossproduct >= 0);  % allow two consecutive edges to be parallel
    t2 = sum(zcrossproduct <= 0);  
    ret = t1 == N || t2 == N;
    
end