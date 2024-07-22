function [PopObj,g] = SMOP8(X,theta)
[N,D] = size(X);
M = 2;
K = ceil(theta*(D-M+1));
g = sum(g3(X(:,M:M+K-1),mod(X(:,M+1:M+K)+pi,2)),2) + sum(g3(X(:,M+K:end-1),X(:,M+K+1:end)*0.9),2);
PopObj = repmat(1+g/(D-M+1),1,M).*fliplr(cumprod([ones(N,1),cos(X(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(X(:,M-1:-1:1)*pi/2)];end
function g = g3(x,t)
    g = 4-(x-t)-4./exp(100*(x-t).^2);
end