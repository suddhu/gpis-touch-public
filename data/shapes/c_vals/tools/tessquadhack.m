function res = tessquadhack(fun,order,data,tess,ric)
% tessquad: numerical integration of a function over an n-d tessellated domain 
% usage 1: res = tessquad(fun,order,data)
% usage 2: res = tessquad(fun,order,data,tess)
%
% REQUIRES the function simplexquad as found here:
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=9435&objectType=file
%
% arguments: (input)
%  fun - function (anonymous, inline, nested function, or m-file)
%        fun must operate on rows of a matrix as points in the
%        n-d data space.
%
%  order   - Order of the Gaussian quadrature to use
%
%  data - nxp array of data (p is the dimension of the space)
%
%  tess - (OPTIONAL) tessellation as delaunayn might produce.
%        If tess is not supplied, then delaunayn will be used
%        to build this array. If delaunayn is used, then the
%        tessellation will cover the convex hull of the data.
%
% arguments: (output)
%  res - integral of fun over the domain defined by data and
%       (optionally) tess.
%
% Example usage 1: Integrate the function x1*x2, over the unit
% square. The result should be 1/4.
%
%  data = [0 0;0 1;1 0;1 1];
%  tess = [1 2 4;1 3 4];
%  res = tessquad(@(x) x(:,1).*x(:,2) ,2,data,tess)
%
% res =
%       0.25
%
% Example usage 2: Integrate the function x1*x2*x3, over the
% unit cube as sampled by 1000 random points. The result should
% approach 1/8 as an asymptote as the number of points grows
% large.
%
%  data = rand(1000,3);
%  res = tessquad(@(x) prod(x,2),2,data)
%
% res =
%      0.11723

% check for the existence of simplexquad  % very expensive
%if ~exist('simplexquad','file')
%  error('Please download simplexquad.m from the file exchange. URL is in the help.');
%end

% number of points, dimension of the space 
%[n,p] = size(data);

% check for errors in order
%if (nargin<2) || isempty(order) || ~isnumeric(order) || ...
%   (length(order)>1) || (fix(order)~=order)
%  error 'order must be a scalar integer'
%end

% do we need to build the tessellation?
%if (nargin<4) || isempty(tess)
%  tess = delaunayn(data);
%end
nt = size(tess,1);

% build nodes and weights for a single "standard"
% simplex with simplexquad
%if order ==2 & p==2
    X0 = [0.280019915499074 0.644948974278318;0.666390246014702 0.155051025721682;0.0750311102226081 0.644948974278318;0.178558728263616 0.155051025721682];
    W0 = [0.0909793091280115;0.159020690871989;0.0909793091280115;0.159020690871989];
    m = 4;
%else
%    V0 = eye(p+1,p);
%    [X0,W0] = simplexquad(order,V0);
%    m = length(W0);
%end

res = 0;

for i = 1:nt
  % i'th simplex
  vert = data(tess(i,:),:);
  v1 = vert(1,:);
  vert = vert(2:end,:) - [v1;v1]; %repmat(v1,p,1);
  
  % transform X0 and W0 from the unit simplex to
  % the current simplex
  W = W0*abs(vert(1,1)*vert(2,2) - vert(1,2)*vert(2,1));
  X = X0*vert + [v1;v1;v1;v1]; %repmat(v1,m,1);
  
  %f = feval(fun,X);
  f = integrand(X,ric);
  res = res + W'*f(:);
end

