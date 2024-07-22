function CrowdDis = CrowdingDistance(PopObj,FrontNo)
% Calculate the crowding distance of each solution front by front

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M]    = size(PopObj);%N是粒子个数，M是目标函数的个数
    CrowdDis = zeros(1,N);
    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts) %对每个rank等级的粒子都有
        Front = find(FrontNo==Fronts(f));
        Fmax  = max(PopObj(Front,:),[],1);% 目标函数值
        Fmin  = min(PopObj(Front,:),[],1);
        for i = 1 : M  % 对f等级的粒子的每一个目标函数都有
            [~,Rank] = sortrows(PopObj(Front,i)); % 将等级f内的所有粒子按照目标函数M进行排序
            CrowdDis(Front(Rank(1)))   = inf;  %保留其中最大的，和最小的
            CrowdDis(Front(Rank(end))) = inf;
            for j = 2 : length(Front)-1
                CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
            end
        end
    end
end