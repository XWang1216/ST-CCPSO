function [PopObj,g] = SMOP1(X,theta)
[N,D] = size(X);
M = 2;
K = ceil(theta*(D-M+1));%1
g = sum(g1(X(:,M:M+K-1),pi/3),2) + sum(g2(X(:,M+K:end),0),2);
PopObj = repmat(1+g/(D-M+1),1,M).*fliplr(cumprod([ones(N,1),X(:,1:M-1)],2)).*[ones(N,1),1-X(:,M-1:-1:1)];
%cumprod «¿€≥À
end
function g = g1(x,t)
g=(x-t).^2;
end
function g = g2(x,t)
    g = 2*(x-t).^2 + sin(2*pi*(x-t)).^2;
end