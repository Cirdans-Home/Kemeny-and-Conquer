function k = recursivekemenyhodlr(P,pi)
%KEMENYRANKSTRUCTURE uses rank-structured matrices to perform the divide
%and conquer strategy.
%   INPUT:  P stochastic matrix, if the matrix is not in hodlr format it
%           get converted on input
%           pi stationary vector of the matrix P
%   OUTPUT: k Kemeny constant (approximation within hodlr format accuracy)

if ~isa(P,"hodlr")
    P = hodlr(P);
end

m = size(P.A11,1);
n = size(P,1);

if n <= 500
    k = kemenydirect(P);
    return
else
    pi1 = pi(1:m)/sum(pi(1:m));
    pi2 = pi(m+1:n)/sum(pi(m+1:n));
    I2 = hodlr(speye(m,m));
    I1 = hodlr(speye(n-m,n-m));
    % Build censored chains
    [L,U] = lu(I1-P.A22);
    P1 = P.A11 + hodlr('low-rank',P.U12,(P.V12.'*hodlr('low-rank',U\(L\P.U21),P.V21)).');
    P2 = P.A22 + hodlr('low-rank',P.U21,(P.V21.'*hodlr('low-rank',(I2-P.A11)\P.U12,P.V12)).');
    % Compute theta
    x = U\(L\ones(n-m,1));
    y = (I2 - P1 - hodlr('low-rank',ones(size(pi1)),pi1))\(P.U12*(P.V12'*x));
    y = x+U\(L\((P.U21*(P.V21'*y))));
    theta = pi2'*(x+y);
    gamma = sum(pi(1:m))*theta-sum(pi(m+1:n));
    k1 = recursivekemenyhodlr(P1,pi1);
    k2 = recursivekemenyhodlr(P2,pi2);
    % Kemeny constant by update
    k = k1+k2+gamma;
end



end