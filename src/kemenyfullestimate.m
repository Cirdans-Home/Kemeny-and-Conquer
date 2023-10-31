function k = kemenyfullestimate(P,samples,varargin)
%KEMENYFULLESTIMATE stochastic estimation of Kemeny constant using the
%Hutchplusplus algorithm.
%   INPUT: P stochastic matrix
%          samples matrix-vector product budget for the Hutch++ application
%          "ITERATIVE" ("ITERATIVE2"/"ITERATIVE3") or "DIRECT" to select 
%          the solution strategy for the computation of the matrix-vector 
%          product oracle.
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
        % Incomplete LU discarding rank-1 update
        A = @(x) x - P*x + e*(e'*x)/n;
        [L,U] = ilu(I - P);
        oracle = @(B) colapply(@(x) Ainv(A,x,L,U),B);
    case "ITERATIVE2"
        % Incomplete LU with rank-1 update
        A = @(x) x - P*x + e*(e'*x)/n;
        [L,U] = ilu(I - P);
        Le = L\(e/sqrt(n));
        Ue = ((e/sqrt(n))'/U)';
        d = Ue'*Le;
        Uinv = @(x) U\(x - (Le*(Ue'*x))/(1+d));
        oracle = @(B) colapply(@(x) Ainv(A,x,L,Uinv),B);
    case "ITERATIVE3"
        % Neumann-Series
        [pi,~] = eigs(P',1,"largestabs","MaxIterations",10000);
        pi = pi/sum(pi);
        A = @(x) x - P*x + e*(pi'*x);
        Prec = @(x) x + P*(P*(P*x)) - 3*e*(pi'*x);
        oracle = @(B) colapply(@(x) Ainv(A,x,Prec,[]),B);
end
k = hutchplusplus(oracle,samples,n);
k = k - 1;
end

function y = Ainv(A,x,L,U)
%AINV Solve linear systems with BiCGstab
    [y,flag,relres] = bicgstab(A,x,1e-3,size(x,1),L,U);
    if flag ~= 0
        warning("BiCGstab flag %d relres on exit is %e",flag,relres);
    end
end

function V = colapply(A,B)
%COLAPPLY Applies A to all columns of B

V = zeros(size(B));
for i=1:size(B,2)
    V(:,i) = A(B(:,i));
end

end
