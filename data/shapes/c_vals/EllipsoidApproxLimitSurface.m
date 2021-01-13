% This is a more general limit surface

classdef EllipsoidApproxLimitSurface < LimitSurface
  properties    
    vertex     % rowvec
  end
  
  methods
    function obj = EllipsoidApproxLimitSurface(M, rc)
      obj@LimitSurface();
      if  nargin > 1
        obj.rc = rc;
        obj = obj.setM(M);
      end
    end
    function obj = setM(obj, M)
      obj.M = M;
      if(~isConvex(M))
        fprintf('Warning: the shape is not convex, will use delauny triangles\n');
        [obj.vertex, obj.tess] = divideObjectMesh2DData(M);
      else
        [obj.vertex, obj.tess] = divideObjectData(M);
      end
      obj = obj.createByEllipsoidApprox();      
    end
    function [obj, x, y, z, mmax, fmax] = createByEllipsoidApprox(obj)
      M = obj.M;
      
      %if isa(M, 'TaylorVar')
        mmax = obj.mus * obj.pressure * mfIntegral_matlab(obj.rc, obj.vertex, obj.tess);
      %else
      %  mmax = obj.mus * obj.pressure * mfIntegral(obj.rc, obj.vertex, obj.tess);
      %end
      fmax = obj.mus * obj.pressure * polyarea(M(1,:)', M(2,:)');
      
      if nargout > 1
          [x, y, z] = ellipsoid(0,0,0,fmax,fmax,mmax,30);
      end
      
      obj.mmax = mmax;
      obj.fmax = fmax;
      obj.c = mmax / fmax;
    end
    
    % pushit is the f function in the paper
    %
    % x: current object location in global frame, x,y,theta
    % M: shape described in object frame
    % push_pt: the point where contact happens z_t=[x,y] in paper (in global frame)
    % push_end: z_t+1 =[x,y] (in global frame)

    % x_new represented in global frame

    function [x_new] = pushit(obj,x,push_pt,push_end)
      c = obj.c;
      vp = rotate_to_frame2d(push_end - push_pt,x);
      xc = transform_to_frame2d(push_pt,x);
      
      denom = (c^2 + xc(1)^2 + xc(2)^2);
      % vx, vy, omega are in object frame
      vx = ((c^2+xc(1)^2) * vp(1) + xc(1)*xc(2)*vp(2)) / denom;
      vy = (xc(1)*xc(2)*vp(1) + (c^2+xc(2)) * vp(2)) / denom;
      omega = (xc(1)*vy - xc(2)*vx) / (c^2);
      
      %fprintf('x(%f %f %f), xc(%f %f) vp(%f %f), vx, vy, omega denom  %f %f %f %f push_pt %f %f push_end %f %f\n', x, xc, vp, vx, vy, omega, denom, push_pt, push_end)
      v_global = rotate_back_frame2d([vx, vy]',x);
      x_new = x + [v_global; omega];
      
    end
    
    
    function [x_new, xc, slip, normvec_now, incontact] = pushit_slip_cylinder(obj,x,xcylinder,xcylinder_end, toviz, xc, normvec)
    % pushit is the f function in the paper
    %
    % x: current object location in global frame, x,y,theta
    % M: shape described in object frame
    % push_pt: the point where contact happens z_t=[x,y] in paper (in global frame)
    % push_end: z_t+1 =[x,y] (in global frame)

    % x_new represented in global frame

    % push_pt in cylinder
    % push_end
    
    % contact threshold
      threshold = 1e-1;

      hasContactInfo = false;
      if nargin == 7; hasContactInfo = true; end
      
      c = obj.c;
      xcylinder_obj = transform_to_frame2d(xcylinder,x);
      xcylinder_end_obj = transform_to_frame2d(xcylinder_end,x);
      % normvec, xc should be in object frame too
      if hasContactInfo
        [ql, qr, xc, vl, vr, normvec, d] = obj.getMotionConeWrtObject(xcylinder_obj, xc, normvec);
      else
        [ql, qr, xc, vl, vr, normvec, d] = obj.getMotionConeWrtObject(xcylinder_obj);
      end
      slip = true;
      incontact = true;
      
%       should remove this when doing optimization
%       too far away then not in contact
      if(~hasContactInfo && abs(d - obj.radius) > threshold)
        x_new = x;
        incontact = false;
        normvec_now = normvec;  % this will not be useful % should return [Inf Inf]
        return
      end

%       pusher cannot pull
%       should remove this when doing optimization
      opposite = dot(normvec, xcylinder_end - xcylinder) < 0;
      opposite = all(opposite(:));
      if opposite
        x_new = x;
        incontact = false;
        normvec_now = normvec;
        return
      end
      
      %vp  push velocity
      vp = rotate_to_frame2d(xcylinder_end - xcylinder,x);
      
      %% todo move visualization outside
      if toviz
        scale_normal = 0.02;
        scale_cone = 0.01;
        vr_global = rotate_back_frame2d(vr,x);
        vl_global = rotate_back_frame2d(vl,x);
        vr_global = vr_global / norm(vr_global) * scale_cone;
        vl_global = vl_global / norm(vl_global) * scale_cone;
        normvec_global = rotate_back_frame2d(normvec,x) * scale_normal; 
        xc_global = transform_back_frame2d(xc,x); 
        hold on
        % plot motion cone
        plot([0 vr_global(1)]+xc_global(1), [0 vr_global(2)]+xc_global(2));
        plot([0 vl_global(1)]+xc_global(1), [0 vl_global(2)]+xc_global(2));
        vp_global
        % plot normal vector
        plot([0 normvec_global(1)]+xc_global(1), [0 normvec_global(2)]+xc_global(2), 'g');
        % plot object pose
        plot(xc_global(1), xc_global(2), 'b*');
      end
      %%
      
      if cross2d(vr, vp) >= 0 && cross2d(vl, vp) <= 0 % inside motion cone
        slip = false;
        denom = (c^2 + xc(1)^2 + xc(2)^2);
        % vx, vy, omega are in object frame
        vx = ((c^2+xc(1)^2) * vp(1) + xc(1)*xc(2)*vp(2)) / denom;
        vy = (xc(1)*xc(2)*vp(1) + (c^2+xc(2)^2) * vp(2)) / denom;
        omega = (xc(1)*vy - xc(2)*vx) / (c^2);
        %fprintf('Inside motioncone x(%f %f %f), xc(%f %f) vp(%f %f), vx, vy, omega denom  %f %f %f %f xcylinder_pt %f %f xcylinder_end %f %f\n', x, xc, vp, vx, vy, omega, denom, xcylinder, xcylinder_end)
      else
        slip = true;
        if cross2d(vp, vr) > 0
          vb = vr;
          qb = qr;
        else
          vb = vl;
          qb = ql;
        end
        
        kappa = (vp(1:2)' * normvec) / (vb(1:2)' * normvec);
        assert(~any(isnan(kappa)), 'nan happened here');
        q = kappa(1) * qb;
        vx = q(1);
        vy = q(2);
        omega = q(3);
      end
      
      nextM = rotate_to_frame2d(obj.M,[vx, vy, omega]');
      
      if ~hasContactInfo
        [dist,~,~,contact_type,normvec_end] = ...
          p_poly_dist(xcylinder_end_obj(1), xcylinder_end_obj(2), nextM(1,:), nextM(2,:));
        if contact_type == 2 % edge
          corrected_v = (obj.radius - dist) * normvec_end' + [vx, vy]';
          vx = corrected_v(1);
          vy = corrected_v(2);
        end
        normvec_now = normvec;
      else
        normvec_now = normvec;
      end
      
      v_global = rotate_back_frame2d([vx, vy]',x);
      
      x_new = x + [v_global; omega];
      
    end
  end
end
