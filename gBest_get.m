function [gBest] = gBest_get(A,value_A)
% �ҵ�gBest 
% ����������һ���µĺ���f_gBest����value_A���룬������С����gBest
%       ���Ѿ����µĵ�����Ѱ��gBest����Ӧ��˳��������µ������value_A����ֵ
%���gBest
[M,~] = size(A);
r = [];
%r_matr = [];
f_gBest = [];
gBest = [];
% case{'ZDT3'}
%         [func_dim,~] = info_ZDT3();
%         r = rand(1,2);%������
%         r_matr = repmat(r,M,1);
%         [value_A,~,~] = ZDT3(A); %���A��ֵ
%         %���㺯��f(��������
%         f_gBest = sum(r_matr.*value_A,2)/sum(r);
%         %�ҵ�f����С��ֵ������
%         [~,index] = min(f_gBest);
%         %��ֵ��gBest
%         gBest = A(index,:);
r = rand(1,size(value_A,2));
r_matr = repmat(r,M,1);
f_gBest = sum(r_matr.*value_A,2)/sum(r);
[~,index] = min(f_gBest);
gBest = A(index,:);
end