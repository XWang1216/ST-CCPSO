function [PopObj,g] = SMOP3(X,theta)
[N,D] = size(X);
M = 2;
K = ceil(theta*(D-M+1));
g = sum(g1(X(:,M:M+K-1),pi/3),2);
for i = 1 : ceil((D-M-K+1)/10)
    temp = 50 - sum(g1(X(:,M+K+(i-1)*10:min(M+K+i*10-1,end)),0),2);
    g(temp<50) = g(temp<50) + temp(temp<50);
end
PopObj = repmat(1+g/(D-M+1),1,M).*fliplr(cumprod([ones(N,1),X(:,1:M-1)],2)).*[ones(N,1),1-X(:,M-1:-1:1)];
end
function g = g1(x,t)
    g = (x-t).^2;
end