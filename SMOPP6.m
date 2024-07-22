function [temp,g] = SMOPP6(X,theta)
[N,D] = size(X);
M = 2;
K = ceil(theta*(D-M+1));
g = g4(X(:,M:end),repmat(linspace(0,1,D-M+1),N,1));
[g,rank] = sort(g,2);%由低到高排序.rank是g的行索引
temp = false(size(rank));%开始全部都是F
for i = 1 : size(rank,1)%每一行中x=0的temp=true，对每一个粒子有
    temp(i,X(i,M-1+rank(i,:))==0) = true;
end
temp(:,1:K) = false;%前K小的temp重新为false
g(temp) = 0;%temp是true的g置0
temp = [false,temp];
g = sum(g,2);
%PopObj = repmat(1+g/(D-M+1),1,M).*fliplr(cumprod([ones(N,1),1-cos(X(:,1:M-1)*pi/2)],2)).*[ones(N,1),1-sin(X(:,M-1:-1:1)*pi/2)];
end
function g = g4(x,t)
g = (x-pi/3).^2 + t.*sin(6*pi*(x-pi/3)).^2;
end