%% KernelFun
function KernelMat = KernelFun2D(x1, x2, R)
    xdatasize = size(x1,1);
    ydatasize = size(x2,1);
    KernelCell = cell(xdatasize,ydatasize);

    for i=1:xdatasize
        for j=1:ydatasize
            KernelCell{i,j}=KernelBlock12D(x1(i,:),x2(j,:),R); % this 
        end
    end
    KernelMat = cell2mat(KernelCell);
end

%% KernelBlock1
% Equation 2 gives us a 3 x 3 matrix 
function KerMat = KernelBlock12D(x1,x2,R)
    N = size(x1, 2); 
    if norm(x1-x2)>1e-4
        KerMat = zeros(N,N);
        KerMat(1,1) = 2*norm(x1-x2)^3 - 3*R*norm(x1-x2)^2 + R^3; % cov(di, dj) 
        KerMat(1,2) = 2*3*norm(x1-x2)*(-x1(1) + x2(1)) - 3*R*2*(-x1(1) + x2(1)); % cov(di, w1j)
        KerMat(1,3)= 2*3*norm(x1-x2)*(-x1(2) + x2(2))-3*R*2*(-x1(2) + x2(2)); % cov(di, w2j)

        KerMat(2,1) = -KerMat(1,2); %2*3*norm(x1-x2)*(x1(1)-x2(1))-3*R*2*(x1(1)-x2(1)); % cov(w1i, dj)
        KerMat(2,2) = 2*3*(-x1(1)+x2(1))/norm(x1-x2)*(x1(1)-x2(1)) + 2*3*norm(x1-x2)*(-1) - 3*R*2*(-1); % cov(w1i, w1j)
        KerMat(2,3) = 2*3*(-x1(2)+x2(2))/norm(x1-x2)*(x1(1)-x2(1)); % cov(w1i, w2j)

        KerMat(3,1) = -KerMat(1,3); %2*3*norm(x1-x2)*(x1(2)-x2(2))-3*R*2*(x1(2)-x2(2)); % cov(w2i, dj)
        KerMat(3,2) = KerMat(2,3); %2*3*(-x1(1)+x2(1))/norm(x1-x2)*(x1(2)-x2(2)); % cov(w2i, w1j)
        KerMat(3,3) = 2*3*(-x1(2)+x2(2))/norm(x1-x2)*(x1(2)-x2(2))+2*3*norm(x1-x2)*(-1)-3*R*2*(-1); % cov(w2i, w2j)
    else
        % edge case? 
        KerMat = zeros(N,N);
        KerMat(1,1) = 2*norm(x1-x2)^3 - 3*R*norm(x1-x2)^2 + R^3;
        KerMat(2,2) = 6*R;
        KerMat(3,3) = 6*R;
    end

end