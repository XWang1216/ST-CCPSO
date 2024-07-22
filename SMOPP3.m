function [Temp] = SMOPP3(X,theta)
[N,D] = size(X);
M = 2;
K = ceil(theta*(D-M+1));
%g = sum(g1(X(:,M:M+K-1),pi/3),2);
for i = 1 : ceil((D-M-K+1)/10)
    temp(:,i) = 50 - sum(g1(X(:,M+K+(i-1)*10:min(M+K+i*10-1,end)),0),2);
    Temp(:,i) = temp(:,i)==50;
end
end
function g = g1(x,t)
    g = (x-t).^2;
end