function [A,value_A] = A_get(pBest,value_pBest,A,value_A,NA,xrange,funcname,theta)
%% ����
% ����A��pBest�Ƚϣ��ҵ�FrontNo==1�����ӣ�����µ�A
% �� ���������p����֮ǰ����A��ֵ����ô���ø��º�ĵ�����pֵ
% �� ���a��pֵû��֧���ϵ����ֱ�ӽ�p���뵵��

%% ����A
%�ڴ�
dist_index = [];
%�ϲ��Ƚ�
Compare_Obj = [value_A;value_pBest];%ǰNA����A��������pBest
%�������Ƚ�
[FrontNo,~] = NDSort(Compare_Obj,Inf);
%[FrontNo,MM] = NDSort(Compare_Obj,100);
compare_x = [A;pBest];
%�ҳ�front
FrontNo_index = find(FrontNo == 1);
%FrontNo_index = find(FrontNo <= MM);
%�����ɵ�A��value
A = compare_x(FrontNo_index,:);
value_A = Compare_Obj(FrontNo_index,:);
%value_A = Compare_Obj(FrontNo_index,:);
[A,index] = unique(A,'rows');
value_A = value_A(index,:);

%% ����A
%ʹ��ELS��A�����Ŷ�
ELS_A = A;
for i=1:size(A,1)%��ȫ��������
    
    j = ceil(rand*size(A,2));
    ELS_A(i,j) = ELS_A(i,j)+ (xrange(j,2)-xrange(j,1))*normrnd(0,1);
    if ELS_A(i,j)>xrange(j,2)
        ELS_A(i,j)=xrange(j,2);
    elseif ELS_A(i,j)<xrange(j,1)
        ELS_A(i,j)=xrange(j,1);
    end

ELS_A(find(abs(ELS_A)<0.005))=0;

end
%���㺯��ֵ
switch funcname
    case {'SMOP1'}
        [Obj_fitness,~] = SMOP1(ELS_A,theta);
    case {'SMOP2'}
        [Obj_fitness,~] = SMOP2(ELS_A,theta);
    case {'SMOP3'}
        [Obj_fitness,~] = SMOP3(ELS_A,theta);
    case {'SMOP4'}
        [Obj_fitness,~,~] = SMOP4(ELS_A,theta);
    case {'SMOP5'}
        [Obj_fitness,~] = SMOP5(ELS_A,theta);
        case {'SMOP6'}
        [Obj_fitness,~] = SMOP6(ELS_A,theta);
    case {'SMOP7'}
        [Obj_fitness,~] = SMOP7(ELS_A,theta);
        case {'SMOP8'}
        [Obj_fitness,~] = SMOP8(ELS_A,theta);
end
A = [A;ELS_A];
value_A = [value_A;Obj_fitness];
%% ��֦
row_A = size(A,1);
if row_A>NA
    %% NGSA-II
%     [Front_rank,max_rank] = NDSort(value_A,Inf);
%     sum_size = 0;
%     Pruned_A = [];
%     Pruned_vA = [];
%     % rank�Ӻ�<NA rank
%     for i = 1:max_rank
%         rank_size(i) = size(find(Front_rank==i),2);
%         sum_size = sum_size + rank_size(i);
%         if i > 1
%             Pruned_A = [Pruned_A;A(find(Front_rank==(i-1)),:)];
%             Pruned_vA = [Pruned_vA;value_A(find(Front_rank==(i-1)),:)];
%         end
%         if sum_size > NA
%             selected_rank = i;
%             break
%         end
%     end
%     %���������rank = i��ֵ�����ϳ�ά�ȣ���ȥ�п��࣬����Ҫ��ѡһ���ּ���
%     %���ԣ���v2��ѡ���߱�����0�ĸ�����ģ�
%     %���ԣ� NGSA II
%     selected_A = A(find(Front_rank==selected_rank),:);
%     selected_vA = value_A(find(Front_rank==selected_rank),:);
%     CrowdDis = CrowdingDistance(selected_vA,ones(1,rank_size(selected_rank)));
%     % CrowdDisԽ��Խ���ױ�ѡ��
%     [~,index_selected] = sort(CrowdDis,'descend');
%     sizeOfSelected_Ai = NA - size(Pruned_A,1);
%     Pruned_A = [Pruned_A;selected_A(index_selected(1:sizeOfSelected_Ai),:)];
%     Pruned_vA = [Pruned_vA;selected_vA(index_selected(1:sizeOfSelected_Ai),:)];
%     A = Pruned_A;
%     value_A = Pruned_vA;
%% SPEA2
%SPEA��ͨ��������ķ�ʽ��A���гͷ���
%����ÿ���㶼��һ��
%����ƽ��������С�����ಢ��һ��
%�ﵽ��Ŀ�󣬽�ÿһ����ƽ��������С��һ����Ϊ���������ĵ�
%��������������ⲿ�ļ���˵�ǲ��Ѻõģ�ȱ����չ��
% �� ��� SPEA2
% SPEA2��֦�Ļ���˼����NSGAII���ƣ��������ҳ�ǰ�����ϵ���
% �����rank=max���ǲ������ӣ�Ҳ��ѡ��ϡ���
% NSGA-II:ѡ����ϡ��ģ�ǰ���������Ӹ�ά��֮��ļӺ�
% SPEA2��ȥ�����ܵģ���������С�ģ������������֮���ŷʽ����
[FrontNo,MaxFNo] = NDSort(value_A,100);
    Next = false(1,length(FrontNo));%ȫ����false����100֮ǰ��MaxFNo-1ȫ����true��֮�����һ������ѡ���Ե���
    %��NGSAII��˼·��ͬ
    Next(FrontNo<MaxFNo) = true;%�ҳ��������һ���ȼ���֮ǰ���еȼ�����
    
    PopObj = value_A;
    fmax   = max(PopObj(FrontNo==1,:),[],1);
    fmin   = min(PopObj(FrontNo==1,:),[],1);
    PopObj = (PopObj-repmat(fmin,size(PopObj,1),1))./repmat(fmax-fmin,size(PopObj,1),1);%Ϊ��ѡ������һ������һ

    %% Select the solutions in the last front
    Last = find(FrontNo==MaxFNo);
    del  = Truncation(PopObj(Last,:),length(Last)-100+sum(Next));
    Next(Last(~del)) = true;%del==0�ı�ѡ��
    % Population for next generation
    %Population = Population(Next);
    A = A(Next,:);
    value_A= value_A(Next,:);
end

end
function Del = Truncation(PopObj,K)
%�ض�
%ÿ���ҳ�����������֮�������С��һ��ȥ��
% Select part of the solutions by truncation

    N = size(PopObj,1);
    
    %% Truncation
    Distance = pdist2(PopObj,PopObj);%�Գƾ���
    Distance(logical(eye(length(Distance)))) = inf;%�Լ����Լ��ľ�����inf
    Del = false(1,N);%ȫ������false
    while sum(Del) < K %K��Ҫȥ���ĸ���
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);%����һ�����򣬶Գƾ���OK
        %Temp���������ڣ���ÿһ������С��Ԫ���ҳ���������i������Ԫ�ص���С����
        
        [~,Rank] = sortrows(Temp);%��һ���о�����С�ı��true
        %Rank���������ڣ�������Ԫ���е���С�����е���С�����ҳ�����ô��һ��Ԫ�ؾ���Ҫ���Ƴ��ģ���True��
        
        Del(Remain(Rank(1))) = true;
    end
end
