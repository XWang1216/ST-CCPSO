function [PopObj,g] = SMOP7(X,theta)
[N,D] = size(X);
M = 2;
K = ceil(theta*(D-M+1));
g = sum(g2(X(:,M:M+K-1),pi/3),2) + sum(g2(X(:,M+K:end),X(:,[M+K+1:end,M+K])*0.9),2);
PopObj = repmat(1+g/(D-M+1),1,M).*fliplr(cumprod([ones(N,1),cos(X(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(X(:,M-1:-1:1)*pi/2)];
end
function g = g2(x,t)
    g = 2*(x-t).^2 + sin(2*pi*(x-t)).^2;
end