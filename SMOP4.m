function [PopObj,g,sortg] = SMOP4(X,theta)
[N,D] = size(X);
M = 2;
K = ceil(theta*(D-M+1));
[g,sortg] = sort(g3(X(:,M:end),0),2);
g = sum(g(:,1:D-M-K+1),2);
PopObj = repmat(1+g/(D-M+1),1,M).*fliplr(cumprod([ones(N,1),1-cos(X(:,1:M-1)*pi/2)],2)).*[ones(N,1),1-sin(X(:,M-1:-1:1)*pi/2)];
end
function g = g3(x,t)
    g = 4-(x-t)-4./exp(100*(x-t).^2);
end