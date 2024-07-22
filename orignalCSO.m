function [x,v] = orignalCSO(x,v,fi,Obj_fitness,xrange)
[Front_rank,~]= NDSort(Obj_fitness,inf);
R1 = rand;
R2 = rand;
R3 = rand;
mean_sub= mean(x,1);
%% ���ѡ����������
randswarm = randperm(size(x,1));
pair_size = floor(size(x,1)/2);
for i = 1:pair_size
    if Front_rank(randswarm(i)) > Front_rank(randswarm(i+pair_size))
        loser_index = randswarm(i);
        winner_index = randswarm(i+pair_size);
    else
        winner_index = randswarm(i);
        loser_index = randswarm(i+pair_size);
    end
    %loser
    v(loser_index,:) = R1*v(loser_index,:) + R2*(x(winner_index,:)-x(loser_index,:))...
                                 + fi*R3*(mean_sub-x(loser_index,:));
    % x
    x(loser_index,:) = x(loser_index,:) + v(loser_index,:);
                                 
end
%% ���Ʒ�Χ
N = size(x,1);
for irange = 1:N %��ÿһ������
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