function [gBest] = gBest_get(A,value_A)
% 找到gBest 
% 方法：构建一个新的函数f_gBest，将value_A带入，其中最小的是gBest
%       在已经更新的档案中寻找gBest所以应该顺便输出更新档案后的value_A的数值
%获得gBest
[M,~] = size(A);
r = [];
%r_matr = [];
f_gBest = [];
gBest = [];
% case{'ZDT3'}
%         [func_dim,~] = info_ZDT3();
%         r = rand(1,2);%列向量
%         r_matr = repmat(r,M,1);
%         [value_A,~,~] = ZDT3(A); %获得A的值
%         %计算函数f(列向量）
%         f_gBest = sum(r_matr.*value_A,2)/sum(r);
%         %找到f中最小的值的索引
%         [~,index] = min(f_gBest);
%         %赋值给gBest
%         gBest = A(index,:);
r = rand(1,size(value_A,2));
r_matr = repmat(r,M,1);
f_gBest = sum(r_matr.*value_A,2)/sum(r);
[~,index] = min(f_gBest);
gBest = A(index,:);
end