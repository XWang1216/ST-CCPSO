function [A,value_A] = A_get(pBest,value_pBest,A,value_A,NA,xrange,funcname,theta)
%% 解释
% 档案A与pBest比较，找到FrontNo==1的粒子，组成新的A
% ？ 如果出现有p优于之前档案A的值，那么就用更新后的档案和p值
% ？ 如果a与p值没有支配关系，就直接将p归入档案

%% 生成A
%内存
dist_index = [];
%合并比较
Compare_Obj = [value_A;value_pBest];%前NA行是A，后面是pBest
%将其带入比较
[FrontNo,~] = NDSort(Compare_Obj,Inf);
%[FrontNo,MM] = NDSort(Compare_Obj,100);
compare_x = [A;pBest];
%找出front
FrontNo_index = find(FrontNo == 1);
%FrontNo_index = find(FrontNo <= MM);
%新生成的A与value
A = compare_x(FrontNo_index,:);
value_A = Compare_Obj(FrontNo_index,:);
%value_A = Compare_Obj(FrontNo_index,:);
[A,index] = unique(A,'rows');
value_A = value_A(index,:);

%% 更新A
%使用ELS对A进行扰动
ELS_A = A;
for i=1:size(A,1)%对全部的粒子
    
    j = ceil(rand*size(A,2));
    ELS_A(i,j) = ELS_A(i,j)+ (xrange(j,2)-xrange(j,1))*normrnd(0,1);
    if ELS_A(i,j)>xrange(j,2)
        ELS_A(i,j)=xrange(j,2);
    elseif ELS_A(i,j)<xrange(j,1)
        ELS_A(i,j)=xrange(j,1);
    end

ELS_A(find(abs(ELS_A)<0.005))=0;

end
%计算函数值
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
%% 剪枝
row_A = size(A,1);
if row_A>NA
    %% NGSA-II
%     [Front_rank,max_rank] = NDSort(value_A,Inf);
%     sum_size = 0;
%     Pruned_A = [];
%     Pruned_vA = [];
%     % rank加和<NA rank
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
%     %多余出来的rank = i的值，加上超维度，减去有空余，所以要挑选一部分加入
%     %策略：（v2挑选决策变量中0的个数多的）
%     %策略： NGSA II
%     selected_A = A(find(Front_rank==selected_rank),:);
%     selected_vA = value_A(find(Front_rank==selected_rank),:);
%     CrowdDis = CrowdingDistance(selected_vA,ones(1,rank_size(selected_rank)));
%     % CrowdDis越大，越容易被选中
%     [~,index_selected] = sort(CrowdDis,'descend');
%     sizeOfSelected_Ai = NA - size(Pruned_A,1);
%     Pruned_A = [Pruned_A;selected_A(index_selected(1:sizeOfSelected_Ai),:)];
%     Pruned_vA = [Pruned_vA;selected_vA(index_selected(1:sizeOfSelected_Ai),:)];
%     A = Pruned_A;
%     value_A = Pruned_vA;
%% SPEA2
%SPEA是通过逐层聚类的方式对A进行惩罚的
%首先每个点都是一类
%类内平均距离最小的两类并做一类
%达到数目后，将每一类中平均距离最小的一类作为保留下来的点
%但是这对于在最外部的集来说是不友好的，缺乏延展性
% 故 提出 SPEA2
% SPEA2剪枝的基本思想与NSGAII类似，都是先找出前沿面上的来
% 后对于rank=max的那部分粒子，也是选最稀疏的
% NSGA-II:选择最稀疏的，前后两个粒子各维度之差的加和
% SPEA2：去掉最密的，即距离最小的，按照与各粒子之间的欧式距离
[FrontNo,MaxFNo] = NDSort(value_A,100);
    Next = false(1,length(FrontNo));%全部是false，在100之前的MaxFNo-1全部置true，之后最后一个，有选择性的找
    %与NGSAII的思路相同
    Next(FrontNo<MaxFNo) = true;%找出满足最后一个等级的之前所有等级的数
    
    PopObj = value_A;
    fmax   = max(PopObj(FrontNo==1,:),[],1);
    fmin   = min(PopObj(FrontNo==1,:),[],1);
    PopObj = (PopObj-repmat(fmin,size(PopObj,1),1))./repmat(fmax-fmin,size(PopObj,1),1);%为了选择最后的一个，归一

    %% Select the solutions in the last front
    Last = find(FrontNo==MaxFNo);
    del  = Truncation(PopObj(Last,:),length(Last)-100+sum(Next));
    Next(Last(~del)) = true;%del==0的被选入
    % Population for next generation
    %Population = Population(Next);
    A = A(Next,:);
    value_A= value_A(Next,:);
end

end
function Del = Truncation(PopObj,K)
%截断
%每次找出与其他粒子之间距离最小的一个去掉
% Select part of the solutions by truncation

    N = size(PopObj,1);
    
    %% Truncation
    Distance = pdist2(PopObj,PopObj);%对称矩阵
    Distance(logical(eye(length(Distance)))) = inf;%自己跟自己的距离是inf
    Del = false(1,N);%全部都是false
    while sum(Del) < K %K是要去掉的个数
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);%按照一行排序，对称矩阵OK
        %Temp的意义在于：将每一行中最小的元素找出来，这是i跟其他元素的最小距离
        
        [~,Rank] = sortrows(Temp);%将一行中距离最小的编程true
        %Rank的意义在于，将所有元素中的最小距离中的最小距离找出，那么这一个元素就是要被移除的（置True）
        
        Del(Remain(Rank(1))) = true;
    end
end
