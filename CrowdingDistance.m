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

    [N,M]    = size(PopObj);%N�����Ӹ�����M��Ŀ�꺯���ĸ���
    CrowdDis = zeros(1,N);
    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts) %��ÿ��rank�ȼ������Ӷ���
        Front = find(FrontNo==Fronts(f));
        Fmax  = max(PopObj(Front,:),[],1);% Ŀ�꺯��ֵ
        Fmin  = min(PopObj(Front,:),[],1);
        for i = 1 : M  % ��f�ȼ������ӵ�ÿһ��Ŀ�꺯������
            [~,Rank] = sortrows(PopObj(Front,i)); % ���ȼ�f�ڵ��������Ӱ���Ŀ�꺯��M��������
            CrowdDis(Front(Rank(1)))   = inf;  %�����������ģ�����С��
            CrowdDis(Front(Rank(end))) = inf;
            for j = 2 : length(Front)-1
                CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
            end
        end
    end
end