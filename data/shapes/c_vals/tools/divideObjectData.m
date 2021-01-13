
function [p,t] = divideObjectData(M)
  nvert = size(M,2);
  %Mtri = zeros(2,(nvert-2)*3);
  t = zeros(nvert-2,3);
  for i=1:nvert-2
    t(i,1) = 1;
    t(i,2) = i+1;
    t(i,3) = i+2;
  end
  p = M';
end
