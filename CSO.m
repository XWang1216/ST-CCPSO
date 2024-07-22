function [x,v] = CSO(x,v,w,c1,c2,fi,pBest,gBest,func_dim,xrange,Obj_fitness)
%% 分种群
% 使用kmeans进行聚类
[idx,~] = kmeans(x,3);
fd1 = find(idx==1);
fd2 = find(idx==2);
fd3 = find(idx==3);
%在每个种群中进行竞争，在更新x v
%注意需要计算lbest 和种群均值
%先计算均值和选取lbest
mean_sub(1,:) = mean(x(idx==1,:),1);%每一个维度的平均
mean_sub(2,:) = mean(x(idx==2,:),1);
mean_sub(3,:) = mean(x(idx==3,:),1);

%% 竞争与更新
% 种群1
% 分成种群
subswarm1_x = x(fd1,:);
subswarm10 = sum(subswarm1_x==0,2);
subswarm1_v = v(fd1,:);
subswarm1_pbest = pBest(fd1,:);

[Front_rank1,~]= NDSort(Obj_fitness(fd1,:),inf);
subswa = subswarm1_x(find(Front_rank1==1),:);
lbest(1,:) = subswa(randi(size(subswa,1)),:);

% 随机排序randperm
randswarm1 = randperm(size(subswarm1_x,1));%0的数量
% NDSort(Front_rank1) 0的数量
sub1_pair_size = floor(size(subswarm1_x,1)/2);
% 更新种群 pbest gbest（按照原序列）
for i = 1:sub1_pair_size
    if Front_rank1(randswarm1(i)) > Front_rank1(randswarm1(i+sub1_pair_size))
        loser_index = randswarm1(i);
        winner_index = randswarm1(i+sub1_pair_size);
    elseif Front_rank1(randswarm1(i)) < Front_rank1(randswarm1(i+sub1_pair_size))
        winner_index = randswarm1(i);
        loser_index = randswarm1(i+sub1_pair_size);
    else
        if subswarm10(randswarm1(i)) > subswarm10(randswarm1(i+sub1_pair_size))
            winner_index = randswarm1(i);
            loser_index = randswarm1(i+sub1_pair_size);
        elseif subswarm10(randswarm1(i)) < subswarm10(randswarm1(i+sub1_pair_size))
            loser_index = randswarm1(i);
            winner_index = randswarm1(i+sub1_pair_size);
        else
            loser_index = randswarm1(i);
            winner_index = randswarm1(i+sub1_pair_size);
        end
    end
    %winner
    R2 = rand;
    subswarm1_v(winner_index,:) = w*subswarm1_v(winner_index,:) + c1*rand*(subswarm1_pbest(winner_index,:)-subswarm1_x(winner_index,:))...
                        + c2*R2*(gBest-subswarm1_x(winner_index,:));
    %loser
    subswarm1_v(loser_index,:) = rand*subswarm1_v(loser_index,:) + rand*(subswarm1_x(winner_index,:)-subswarm1_x(loser_index,:))...
                                 + fi*rand*(mean_sub(1,:)-subswarm1_x(loser_index,:));
    % x
    subswarm1_x(winner_index,:) = subswarm1_x(winner_index,:) + subswarm1_v(winner_index,:);
    subswarm1_x(loser_index,:) = subswarm1_x(loser_index,:) + subswarm1_v(loser_index,:);
                                 
end


% 种群2
% 分成种群
subswarm2_x = x(fd2,:);
subswarm20 = sum(subswarm2_x==0,2);
subswarm2_v = v(fd2,:);
subswarm2_pbest = pBest(fd2,:);

[Front_rank2,~]= NDSort(Obj_fitness(fd2,:),inf);
subswa = subswarm2_x(find(Front_rank2==1),:);
lbest(2,:) = subswa(randi(size(subswa,1)),:);

% 随机排序randperm
randswarm1 = randperm(size(subswarm2_x,1));
% NDSort(Front_rank1) 0的数量
sub1_pair_size = floor(size(subswarm2_x,1)/2);
% 更新种群 pbest gbest（按照原序列）
for i = 1:sub1_pair_size
    if Front_rank2(randswarm1(i)) > Front_rank2(randswarm1(i+sub1_pair_size))
        loser_index = randswarm1(i);
        winner_index = randswarm1(i+sub1_pair_size);
    elseif Front_rank2(randswarm1(i)) < Front_rank2(randswarm1(i+sub1_pair_size))
        winner_index = randswarm1(i);
        loser_index = randswarm1(i+sub1_pair_size);
    else
        if subswarm20(randswarm1(i)) > subswarm20(randswarm1(i+sub1_pair_size))
            winner_index = randswarm1(i);
            loser_index = randswarm1(i+sub1_pair_size);
        elseif subswarm20(randswarm1(i)) < subswarm20(randswarm1(i+sub1_pair_size))
            loser_index = randswarm1(i);
            winner_index = randswarm1(i+sub1_pair_size);
        else
            loser_index = randswarm1(i);
            winner_index = randswarm1(i+sub1_pair_size);
        end
    end
    %winner
    R2 = rand;
    subswarm2_v(winner_index,:) = w*subswarm2_v(winner_index,:) + c1*rand*(subswarm2_pbest(winner_index,:)-subswarm2_x(winner_index,:))...
                        + c2*R2*(gBest-subswarm2_x(winner_index,:));
    %loser
    subswarm2_v(loser_index,:) = rand*subswarm2_v(loser_index,:) + rand*(subswarm2_x(winner_index,:)-subswarm2_x(loser_index,:))...
                                 + fi*rand*(mean_sub(2,:)-subswarm2_x(loser_index,:));
    % x
    subswarm2_x(winner_index,:) = subswarm2_x(winner_index,:) + subswarm2_v(winner_index,:);
    subswarm2_x(loser_index,:) = subswarm2_x(loser_index,:) + subswarm2_v(loser_index,:);
                                 
end
% 种群3
% 分成种群
subswarm3_x = x(fd3,:);
subswarm30 = sum(subswarm3_x==0,2);
subswarm3_v = v(fd3,:);
subswarm3_pbest = pBest(fd3,:);

[Front_rank3,~]= NDSort(Obj_fitness(fd3,:),inf);
subswa = subswarm3_x(find(Front_rank3==1),:);
lbest(3,:) = subswa(randi(size(subswa,1)),:);

% 随机排序randperm
randswarm1 = randperm(size(subswarm3_x,1));
% NDSort(Front_rank1) 0的数量
sub1_pair_size = floor(size(subswarm3_x,1)/2);
% 更新种群 pbest gbest（按照原序列）
for i = 1:sub1_pair_size
    if Front_rank3(randswarm1(i)) > Front_rank3(randswarm1(i+sub1_pair_size))
        loser_index = randswarm1(i);
        winner_index = randswarm1(i+sub1_pair_size);
    elseif Front_rank3(randswarm1(i)) < Front_rank3(randswarm1(i+sub1_pair_size))
        winner_index = randswarm1(i);
        loser_index = randswarm1(i+sub1_pair_size);
    else
        if subswarm30(randswarm1(i)) > subswarm30(randswarm1(i+sub1_pair_size))
            winner_index = randswarm1(i);
            loser_index = randswarm1(i+sub1_pair_size);
        elseif subswarm30(randswarm1(i)) < subswarm30(randswarm1(i+sub1_pair_size))
            loser_index = randswarm1(i);
            winner_index = randswarm1(i+sub1_pair_size);
        else
            loser_index = randswarm1(i);
            winner_index = randswarm1(i+sub1_pair_size);
        end
    end
    %winner
    R2 = rand;
    subswarm3_v(winner_index,:) = w*subswarm3_v(winner_index,:) + c1*rand*(subswarm3_pbest(winner_index,:)-subswarm3_x(winner_index,:))...
                        + c2*R2*(gBest-subswarm3_x(winner_index,:));
    %loser
    subswarm3_v(loser_index,:) = rand*subswarm3_v(loser_index,:) + rand*(subswarm3_x(winner_index,:)-subswarm3_x(loser_index,:))...
                                 + fi*rand*(mean_sub(3,:)-subswarm3_x(loser_index,:));
    % x
    subswarm3_x(winner_index,:) = subswarm3_x(winner_index,:) + subswarm3_v(winner_index,:);
    subswarm3_x(loser_index,:) = subswarm3_x(loser_index,:) + subswarm3_v(loser_index,:);
                                 
end
% 将种群123汇总
x(fd1,:) = subswarm1_x;
x(fd2,:) = subswarm2_x;
x(fd3,:) = subswarm3_x;

v(fd1,:) = subswarm1_v;
v(fd2,:) = subswarm2_v;
v(fd3,:) = subswarm3_v;
% x = [subswarm1_x;subswarm2_x;subswarm3_x];
% v = [subswarm1_v;subswarm2_v;subswarm3_v];

subswarm1_x = []; subswarm2_x = []; subswarm3_x = [];
subswarm1_v = []; subswarm2_v =[]; subswarm3_v =[];
subswarm1_pbest =[]; subswarm2_pbest =[]; subswarm3_pbest = [];
%% 限制范围
N = size(x,1);
for irange = 1:N %对每一个粒子
      Upper_flag = xrange(:,2)'<x(irange,:);%都是行向量， 寻找超出最大值的粒子，如果有的话，就是1
      Upper_flag_T = sum(Upper_flag);
      if Upper_flag_T > 0  %说明存在超出上限的粒子
          Upper_index = find(Upper_flag == 1);
          x(irange,Upper_index) = xrange(Upper_index,2);
      end
      
      Low_flag = xrange(:,1)' > x(irange,:); %行向量，寻找低于最小值的粒子，有的话，是1
      Low_flag_T = sum(Low_flag);
      if Low_flag_T > 0
          Low_index = find(Low_flag == 1);
          x(irange,Low_index) = xrange(Low_index,1);
      end 
end
 x(find(abs(x)<0.005))=0;
end