function [x,Obj_fitness,pBest,value_pBest] = pBest_get(x,Obj_fitness,pBest,value_pBest)
% Ѱ��pBest������õ��������x����p���Ǿ͸���
% relat_px��p��x�Ĺ�ϵ�����Ժ���²���ʱ���õ�
[N,~] = size(Obj_fitness);%201
N = N/2;
Compare_Obj = [Obj_fitness;value_pBest];%ǰN�����µ�x����x��FrontNoֵС�ں����FrontNoֵ����ô�͸���pBest
%�������Ƚ�
[FrontNo,~] = NDSort(Compare_Obj,Inf);
relat_px = [];
for i = 1:N
    if FrontNo(i) <= FrontNo(i+2*N) %x����p 1 %%%%%%%%%%%ֻ�����ڲŸ�ֵ�������ᵼ��pBestһ���ҵ�ǰ����Ͳ����ƶ�
        pBest(i,:) = x(i,:);
        value_pBest(i,:) = Obj_fitness(i,:);
        if FrontNo(i+N) <= FrontNo(i)
           pBest(i,:) = x(i+N,:);
           value_pBest(i,:) = Obj_fitness(i+N,:);
            x(i,:) =  x(i+N,:);
            Obj_fitness(i,:) = Obj_fitness(i+N,:);
        end
    else
        if FrontNo(i+N) <= FrontNo(i+2*N)
            pBest(i,:) = x(i+N,:);
            value_pBest(i,:) = Obj_fitness(i+N,:);
             x(i,:) =  x(i+N,:);
            Obj_fitness(i,:) = Obj_fitness(i+N,:);
        end
    end  
end
x = x(1:N,:);
Obj_fitness = Obj_fitness(1:N,:);
end