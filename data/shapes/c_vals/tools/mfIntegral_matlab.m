
% pushpt is 2*1 in object frame
% pushvect is 2*1 in object frame

function ret = mfIntegral_matlab(xic, p, t)
  % ric in push frame
  %ric = [xic, 0]'; 
  %ret = tessquad(@(x)integrand4(x, ric), 2 , p, t);
  ret = tessquadhack([], 2, p, t, xic);
end


