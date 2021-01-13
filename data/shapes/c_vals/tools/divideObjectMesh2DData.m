
function [p,t] = divideObjectMesh2DData(M)
  options.output = false;
  [p,t] = mesh2d(M',[],[], options);
end