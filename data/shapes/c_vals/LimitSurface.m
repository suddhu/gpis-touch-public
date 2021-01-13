classdef LimitSurface
  %LimitSurface A class for generating limit surface (LS) of a polygon, and use it in pushing
  %
  properties
    sampleRate = 0.01;
    mus = 0.25; % friction coefficient between the slider and support
    muc = 0.25; % friction coefficient between the pusher and slider
    M;          % shape, colvecs
    rc; % rotation center
    % for sample based LS generation
    ndegQuad = 2; 
    x0 = [0,0]';   %origin for calculating moment
    sampleLimit=[-0.5 0.5 -0.5 0.5]; % [xmin xmax, ymin ymax]
    tess;          % tessalation used by tessquad() 
    k = [0,0,1]';  % rotation type counter clockwise or CCW (default: CCW)
    pressure = 1;  % pressure distrbution? 
%     pressure = 1.0132e+03;  %     0.8374 mass,     0.0081 area
    
    % for ellipsoidal approximation (Lee and Cutosky)
    fmax;
    mmax;
    c;
    
    % for pusher
    radius = 0.00626/2;  %  real probe1
    %radius = 0.004745;  %  real probe2
    %radius = 0.005; % simulated pusher
  end
  
  methods
    
    function obj = LimitSurface()
      % do nothing
    end
    
    function obj = setRotationDirection(obj, str)
      %setRotationDirection set the rotation to be counter clockwise (CCW) or
      %clockwise (CW)
      %   ADDME(str) str = 'CW' or 'CCW'.
      if strcmp(str,'CCW')
        obj.k = [0,0,1];
      else
        obj.k = [0,0,-1];
      end
    end
  
    function [Mf, Ff] = createByRCSampling_tessquad(obj)
      % make LS of a shape by sampling center of rotation
      % using tessquad for integration
      % create limit surface by sampling rotation center
      limit = obj.sampleLimit;
      
      % sample rotation center
      % Xrc = [minx, ... maxx]
      %       [minx, ... maxx]
      % Yrc = [miny, ... miny]
      %       [maxy, ... maxy]
      [Xrc,Yrc] = meshgrid(limit(1):obj.sampleRate:limit(2), ...
                           limit(3):obj.sampleRate:limit(4)); 
      Xrc = Xrc(:);
      Yrc = Yrc(:);
      
      % moments and forces 
      Mf = zeros(1,length(Xrc));
      Ff = zeros(2,length(Xrc));
      for i=1:length(Xrc)
        % numerical integration over the evlaution points (x) with 2nd
        % order gaussian quadrature, and given shape
        Mf(1,i) = tessquad(@(x) mfIntegrand(obj, x, [Xrc(i), Yrc(i)]'), ...
                           obj.ndegQuad, obj.M, obj.tess);
        Ff(:,i) = tessquad(@(x) fIntegrand(obj, x, [Xrc(i), Yrc(i)]'), ...
                           obj.ndegQuad, obj.M, obj.tess)';
      end
    end
 
    %% used by tessquad
    % xrc: rotation center (colvec)
    % X: evaluation points (rowvecs)
    % ret: rowvecs
    function ret = mfIntegrand(obj,X,xrc)
      % Eq (4) in the paper
      ninput = size(X,1);
      ret = zeros(ninput,1);
      for i=1:ninput
        x = X(i,:)';
        va = cross(obj.k,[x-xrc;0]); % k points in vertical axis
        vnorm = va / norm(va);
        ret(i) = obj.mus * cross2d(x-obj.x0, vnorm(1:2));        % where is pressure? 
      end
    end
    
    function ret = fIntegrand(obj,X,xrc)
      % Eq (3) in the paper
      ninput = size(X,1);
      ret = zeros(ninput,2);
      for i=1:ninput
        x = X(i,:)';
        va = cross(obj.k,[x-xrc;0]); % k points in vertical axis
        vnorm = va / norm(va);
        ret(i,:) = obj.mus * vnorm(1:2)';        % where is pressure? 
      end
    end
    
    %% used by matlab's integral2()
    function ret = mfIntegrand_matlab(obj,X,Y,xrc)
      % X, Y are 2D arrays in meshgrid style.
      % xrc: rowvec
      % Eq (4) in the paper
      ninputx = size(X,1);
      ninputy = size(X,2);
      ret = zeros(ninputx,ninputy);
      for i=1:ninputx
      for j=1:ninputy
        x = X(i,j)';
        y = Y(i,j)';
        va = cross(obj.k,[[x,y]'-xrc;0]);
        vnorm = va / norm(va);
        ret(i,j) = obj.mus * cross2d([x,y]'-obj.x0, vnorm(1:2));
      end
      end
    end
        
    function ret = fxIntegrand_matlab(obj,X,Y,xrc)
      % Eq (3) in the paper
      ninputx = size(X,1);
      ninputy = size(X,2);
      ret = zeros(ninputx,ninputy);
      for i=1:ninputx
      for j=1:ninputy
        x = X(i,j)';
        y = Y(i,j)';
        va = cross(obj.k,[[x,y]'-xrc;0]);
        vnorm = va / norm(va);
        ret(i,j) = obj.mus * vnorm(1);
      end
      end
    end
    
    function ret = fyIntegrand_matlab(obj,X,Y,xrc)
      % Eq (3) in the paper
      ninputx = size(X,1);
      ninputy = size(X,2);
      ret = zeros(ninputx,ninputy);
      for i=1:ninputx
      for j=1:ninputy
        x = X(i,j)';
        y = Y(i,j)';
        va = cross(obj.k,[[x,y]'-xrc;0]);
        vnorm = va / norm(va);
        ret(i,j) = obj.mus * vnorm(2);
      end
      end
    end
    
    %% Compute motion cones 
    % should be moved outside
    function [ql, qr, xc, vl, vr, normvec, dist] = getMotionConeWrtObject(obj, xcylinder, xc, normvec)
    % getMotionConeWrtObject get the motion cone in the object coordinate
    % xc: the push point in object frame
    % xcylinder: the position of cylinder in object frame
    % returns motion cone velocities: vr, vl from Eq (7)
      
    % note returned normvec is rowvec.
    % compute the conatact point between cylinder and the object
    
      hasContactInfo = false;
      if nargin == 7; hasContactInfo = true; end
      
      if hasContactInfo
        dist = norm(xc-xcylinder);
      else
        [dist,x_poly,y_poly,~,normvec] = p_poly_dist(xcylinder(1), xcylinder(2), obj.M(1,:), obj.M(2,:));
        normvec = normvec';

        xc = [x_poly, y_poly]';
      end
      
      theta = atan(obj.muc);
      fr = rotate_back_frame2d(normvec, [0,0,-theta]);
      fl = rotate_back_frame2d(normvec, [0,0,theta]);
      
      pr = [fr; cross2d(xc, fr)];
      pl = [fl; cross2d(xc, fl)];

      qr = [obj.c^2 * fr ; pr(3)];
      ql = [obj.c^2 * fl ; pl(3)];
      
      vr_tmp = cross([0;0;qr(3)], [xc; 0]); % whole body rotational velocity
      vl_tmp = cross([0;0;ql(3)], [xc; 0]);
      vr = qr(1:2) + vr_tmp(1:2);  % add translation 
      vl = ql(1:2) + vl_tmp(1:2);
      
      % vr and vl are from Eq (7) 
    end
  end
  
  methods (Static)
  end
end
