function k = kemenyfullestimate(P,samples,varargin)
%KEMENYFULLESTIMATE stochastic estimation of Kemeny constant using the
%Hutchplusplus algorithm.
%   INPUT: P stochastic matrix
%          samples matrix-vector product budget for the Hutch++ application
%          "ITERATIVE" or "DIRECT" to select the solution strategy for the
%          computation of the matrix-vector product oracle.
%   OUTPUT: k a (randomized) estimate of Kemeny's constant.

n = size(P,1);
I = speye(n);
e = ones(n,1);

if nargin == 3
    type = varargin{1};
else
    type = "ITERATIVE";
end

switch upper(type)
    case "DIRECT"
        oracle = @(x) (I - P + e*e'/n)\x;
    case "ITERATIVE"
        A = @(x) x - P*x + e*(e'*x)/n;
        [L,U] = ilu(I - P);
        oracle = @(B) colapply(@(x) Ainv(A,x,L,U),B);
end
k = hutchplusplus(oracle,samples,n);
k = k - 1;
end

function y = Ainv(A,x,L,U)
%AINV Solve linear systems with GMRES
    [y,~] = gmres(A,x,[],1e-6,size(L,1),L,U);
end

function V = colapply(A,B)
%COLAPPLY Applies A to all columns of B

V = zeros(size(B));
for i=1:size(B,2)
    V(:,i) = A(B(:,i));
end

end
