function [PopObj,g] = SMOP5(X,theta)
[N,D] = size(X);
M = 2;
K = ceil(theta*(D-M+1));
g = sum(g1(X(:,M:end),pi/3).*g2(X(:,M:end),0),2) + abs(K-sum(X(:,M:end)~=0,2));
PopObj = repmat(1+g/(D-M+1),1,M).*fliplr(cumprod([ones(N,1),1-cos(X(:,1:M-1)*pi/2)],2)).*[ones(N,1),1-sin(X(:,M-1:-1:1)*pi/2)];
end
function g = g1(x,t)
    g = (x-t).^2;
end
function g = g2(x,t)
    g = 2*(x-t).^2 + sin(2*pi*(x-t)).^2;
end