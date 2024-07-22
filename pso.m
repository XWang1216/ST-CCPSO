function [x,v] =  pso(x,v,w,c1,c2,pBest,gBest,func_dim,xrange,Obj_fitness)
[idx,~] = kmeans(x,3);
fd1 = find(idx==1);
fd2 = find(idx==2);
fd3 = find(idx==3);

[Front_rank1,~]= NDSort(Obj_fitness(fd1,:),inf);
subswa = x(find(Front_rank1==1),:);
lbest(1,:) = subswa(randi(size(subswa,1)),:);

[Front_rank2,~]= NDSort(Obj_fitness(fd2,:),inf);
subswa = x(find(Front_rank2==1),:);
lbest(2,:) = subswa(randi(size(subswa,1)),:);

[Front_rank3,~]= NDSort(Obj_fitness(fd3,:),inf);
subswa = x(find(Front_rank3==1),:);
lbest(3,:) = subswa(randi(size(subswa,1)),:);



[N,~] = size(x);

%PSO经典更新公式
%w是列，v是一个矩阵  rand是否应该一样？不一样
% v = w*v + c1*rand*(pBest-x)+ c2*(rand*(repmat(gBest,N,1)-x)+(1-rand)*(lbest()-x));
% x = x+v;

for i = 1:N
    R2=rand;
    v = w*v + c1*rand*(pBest-x)+ c2*(R2*(repmat(gBest,100,1)-x)+(1-R2).*(lbest(idx(i),:)-x));
x = x+v;
end
 %判断x是否超出种群范围
 for irange = 1:N
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
end