function ret = integrand(x,ric)
  %pressure  =.01;  % assume uniform distribution
  x_ = x*ric(1)*0;
  x_ = x_ + x;
  x_(:,1) = x(:,1)-ric(1);
  x_(:,2) = x(:,2)-ric(2);
  normm = sqrt(x_(:,1).^2 + x_(:,2).^2);
  
  ret = dot(x, x_, 2) ./ normm;
  
end