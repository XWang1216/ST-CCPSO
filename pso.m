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

%PSO������¹�ʽ
%w���У�v��һ������  rand�Ƿ�Ӧ��һ������һ��
% v = w*v + c1*rand*(pBest-x)+ c2*(rand*(repmat(gBest,N,1)-x)+(1-rand)*(lbest()-x));
% x = x+v;

for i = 1:N
    R2=rand;
    v = w*v + c1*rand*(pBest-x)+ c2*(R2*(repmat(gBest,100,1)-x)+(1-R2).*(lbest(idx(i),:)-x));
x = x+v;
end
 %�ж�x�Ƿ񳬳���Ⱥ��Χ
 for irange = 1:N
      Upper_flag = xrange(:,2)'<x(irange,:);%������������ Ѱ�ҳ������ֵ�����ӣ�����еĻ�������1
      Upper_flag_T = sum(Upper_flag);
      if Upper_flag_T > 0  %˵�����ڳ������޵�����
          Upper_index = find(Upper_flag == 1);
          x(irange,Upper_index) = xrange(Upper_index,2);
      end
      
      Low_flag = xrange(:,1)' > x(irange,:); %��������Ѱ�ҵ�����Сֵ�����ӣ��еĻ�����1
      Low_flag_T = sum(Low_flag);
      if Low_flag_T > 0
          Low_index = find(Low_flag == 1);
          x(irange,Low_index) = xrange(Low_index,1);
      end 
 end
end