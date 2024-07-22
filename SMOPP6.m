function [temp,g] = SMOPP6(X,theta)
[N,D] = size(X);
M = 2;
K = ceil(theta*(D-M+1));
g = g4(X(:,M:end),repmat(linspace(0,1,D-M+1),N,1));
[g,rank] = sort(g,2);%�ɵ͵�������.rank��g��������
temp = false(size(rank));%��ʼȫ������F
for i = 1 : size(rank,1)%ÿһ����x=0��temp=true����ÿһ��������
    temp(i,X(i,M-1+rank(i,:))==0) = true;
end
temp(:,1:K) = false;%ǰKС��temp����Ϊfalse
g(temp) = 0;%temp��true��g��0
temp = [false,temp];
g = sum(g,2);
%PopObj = repmat(1+g/(D-M+1),1,M).*fliplr(cumprod([ones(N,1),1-cos(X(:,1:M-1)*pi/2)],2)).*[ones(N,1),1-sin(X(:,M-1:-1:1)*pi/2)];
end
function g = g4(x,t)
g = (x-pi/3).^2 + t.*sin(6*pi*(x-pi/3)).^2;
end