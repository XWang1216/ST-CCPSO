function [x,Obj_fitness,g] = RDA_x(x,Obj_fitness, xrange,funcname,theta,L,eta,g,iter)
%每个粒子用过RDA生成一个新的x，共生成2*Pop_size个x，返回
[~,f_num] = size(Obj_fitness);
%% RDA,使决策变量稀疏
[row_x,col_x]=size(x);
%对每一个粒子x，变异出新的粒子new_x.新的粒子通过RDA对其中一个目标函数进行优化
%通常new_x要比x更加稀疏
new_x = [];

if funcname=='SMOP3'
    Tempp = SMOPP3(x,theta);
else
    Tempp = zeros(row_x,1);
end
% 新函数 func
if iter == 1
for i = 1:row_x%对每一个粒子而言
    j = ceil(f_num*rand(1));%随机选择一个目标函数进行优化
    [df] = get_gg(j,funcname, x(i,:),theta,col_x,Tempp(i,:));
    
     new_x(i,:) = zeros(1,col_x);
     g(i,:) = df;
     
    for k = 1:col_x
        if k==1
        if L(1) < abs(df(k))
            new_x(i,k) = x(i,k) - eta*df(k);
        end
        else
            if L(2) < abs(df(k))
            new_x(i,k) = x(i,k) - eta*df(k);
            end
        end
    end
end
else
    for i = 1:row_x%对每一个粒子而言
    j = ceil(f_num*rand(1));%随机选择一个目标函数进行优化
    [df] = get_gg(j,funcname, x(i,:),theta,col_x,Tempp(i,:));
    
     new_x(i,:) = zeros(1,col_x);
    for k = 1:col_x
   
            if abs(g(i,k)) < abs(df(k))
            new_x(i,k) = x(i,k) - eta*df(k);
            end
      
    end
    g(i,:) = ((iter-1).*g(i,:) + df) ./iter;
end
end

%% 限制范围
%超出范围的，在维度内随机生成
x_maxmin = (xrange(:,2)-xrange(:,1))';
for irange = 1:row_x %对每一个粒子
      Upper_flag = xrange(:,2)'<new_x(irange,:);%都是行向量， 寻找超出最大值的粒子，如果有的话，就是1
      Upper_flag_T = sum(Upper_flag);
      if Upper_flag_T > 0  %说明存在超出上限的粒子
          Upper_index = find(Upper_flag == 1);
          new_x(irange,Upper_index) = rand(1,size(Upper_index,2)).*x_maxmin(Upper_index);
         %new_x(irange,Upper_index)=x(irange,Upper_index);
      end
      
      Low_flag = xrange(:,1)' > new_x(irange,:); %行向量，寻找低于最小值的粒子，有的话，是1
      Low_flag_T = sum(Low_flag);
      if Low_flag_T > 0
          Low_index = find(Low_flag == 1);
          new_x(irange,Low_index) = rand(1,size(Low_index,2)).*x_maxmin(Low_index);
          %new_x(irange,Low_index)=x(irange,Low_index);
      end 
      new_x(find(abs(new_x)<0.005))=0;
end
 %% 计算new_x的目标函数值
                     switch funcname
                            case {'SMOP1'}
                                 [nObj_fitness,~] = SMOP1(new_x,theta);
                            case {'SMOP2'}
                                 [nObj_fitness,~] = SMOP2(new_x,theta);
                         case{'SMOP3'}
                             [nObj_fitness,~] = SMOP3(new_x,theta);
                         case {'SMOP4'}
                             [nObj_fitness,~] = SMOP4(new_x,theta);
                         case {'SMOP5'}
                             [nObj_fitness,~] = SMOP5(new_x,theta);
                             case {'SMOP6'}
                             [nObj_fitness,~] = SMOP6(new_x,theta);
                         case {'SMOP7'}
                             [nObj_fitness,~] = SMOP7(new_x,theta);
                             case {'SMOP8'}
                             [nObj_fitness,~] = SMOP8(new_x,theta);
                     end
                       
%% 生成新的x
x = [x;new_x];
Obj_fitness = [Obj_fitness;nObj_fitness];
end
function [df] = get_gg(j,funcname,a,theta,D,Temp)
switch funcname
    case {'SMOP1'}
        [~,g] = SMOP1(a,theta);
        K = ceil(theta*(D-2+1));
        if j==1
            df(1) = 1+g/(D-2+1);
            df(2:2+K-1) = 2*a(1)*(a(2:2+K-1)-pi/3)/(D-2+1);
            %df(2:6) = 2*a(1)*(a(2:6)-pi/3)/9;
            df(2+K:D) = 4*a(1).*(a(2+K:D)+pi.*sin(2*pi*a(2+K:D)).*cos(2*pi*a(2+K:D)))/(D-2+1);
        elseif j==2
            df(1) = -1-g/(D-2+1);
            df(2:2+K-1) = 2*(1-a(1))*(a(2:2+K-1)-pi/3)/(D-2+1);
            %df(2:6) = 2*(1-a(1))*(a(2:6)-pi/3)/9;
            df(2+K:D) = 4*(1-a(1)).*(a(2+K:D)+pi.*sin(2*pi*a(2+K:D)).*cos(2*pi*a(2+K:D)))/(D-2+1);
        end
        case {'SMOP2'}
        [~,g] = SMOP2(a,theta);
        K = ceil(theta*(D-2+1));
        if j==1
            df(1) = 1+g/(D-2+1);
            df(2:2+K-1) = 4*a(1).*(a(2:2+K-1)-pi/3+pi*sin(2*pi*(a(2:2+K-1)-pi/3)).*cos(2*pi*(a(2:2+K-1)-pi/3)))/(D-2+1);
            %df(2:6) = 2*a(1)*(a(2:6)-pi/3)/9;
            df(2+K:D) = a(1).*(-1+4*(200*a(2+K:D).*exp(-100*a(2+K:D).^2)))/(D-2+1);
        elseif j==2
            df(1) = -1-g/(D-2+1);
            df(2:2+K-1) = 4.*(1-a(1)).*(a(2:2+K-1)-pi/3+pi.*sin(2*pi*(a(2:2+K-1)-pi/3)).*cos(2*pi*(a(2:2+K-1)-pi/3)))/(D-2+1);
            %df(2:6) = 2*(1-a(1))*(a(2:6)-pi/3)/9;
            df(2+K:D) = (1-a(1)).*(-1+4*(200*a(2+K:D).*exp(-100*a(2+K:D).^2)))/(D-2+1);
        end
    case {'SMOP3'}
        [~,g] = SMOP3(a,theta);
        K = ceil(theta*(D-2+1));
        if j==1
            df(1) = 1+g/(D-2+1);
            df(2:2+K-1) = 2*a(1).*(a(2:2+K-1)-pi/3)/(D-2+1);
            %df(2:6) = 2*a(1)*(a(2:6)-pi/3)/9;
            %df(2+K:D) = 0;
            M=2;
            df(2+K:D) = -2*a(1).*a(2+K:D)/(D-2+1);
            for i = 1 : ceil((D-M-K+1)/10)
                if Temp(i) == true
                   df(M+K+(i-1)*10:min(M+K+i*10-1,end)) =0;
                end
             end
        elseif j==2
            df(1) = -1-g/(D-2+1);
            df(2:2+K-1) = 2*(1-a(1)).*(a(2:2+K-1)-pi/3)/(D-2+1);
            %df(2:6) = 2*(1-a(1))*(a(2:6)-pi/3)/9;
            df(2+K:D) = -2*(1-a(1)).*a(2+K:D)/(D-2+1);
             M=2;
             for i = 1 : ceil((D-M-K+1)/10)
                if Temp(i) == true
                   df(M+K+(i-1)*10:min(M+K+i*10-1,end)) =0;
                end
             end
        end
        case {'SMOP4'}
        [~,g,sortg] = SMOP4(a,theta);
        K = ceil(theta*(D-2+1));
        if j==1
            df(1) = (1+g/(D-2+1))*sin(a(1)*pi/2)*pi/2;
            df(2:D) = (1-cos(a(1)*pi/2)).*(-1+800.*exp(-100*((a(2:D)-pi/3).^2)).*(a(2:D)-pi/3))/(D-2+1);
        elseif j==2
            df(1) = -(1+g/(D-2+1))*cos(a(1)*pi/2)*pi/2;
            df(2:D) = (1-sin(a(1)*pi/2)).*(-1+800.*exp(-100*((a(2:D)-pi/3).^2)).*(a(2:D)-pi/3))/(D-2+1);
        end
        df(sortg(D-K:end)) = 0;
    case {'SMOP5'}
        [~,g] = SMOP5(a,theta);
        K = ceil(theta*(D-2+1));
        if j==1
            df(1) = (1+g/(D-2+1))*sin(a(1)*pi/2)*pi/2;
            df(2:D) = (1-cos(a(1)*pi/2)).*(2*(a(2:D)-pi/3).*(2*a(2:D).^2+sin(2*pi*a(2:D)).^2)+...
                      4*(a(2:D)-pi/3).^2.*(a(2:D).^2+pi*sin(2*pi*a(2:D)).*cos(2*pi*a(2:D))));
            %df(2+K:D) = 0;
            df(2:D) = df(2:D)./(D-2+1);
        elseif j==2
            df(1) = -(1+g/(D-2+1))*cos(a(1)*pi/2)*pi/2;
           df(2:D) = (1-cos(a(1)*pi/2)).*(2*(a(2:D)-pi/3).*(2*a(2:D).^2+sin(2*pi*a(2:D)).^2)+...
                      4*(a(2:D)-pi/3).^2.*(a(2:D).^2+pi*sin(2*pi*a(2:D)).*cos(2*pi*a(2:D))));
            %df(2+K:D) = 0;
            df(2:D) = df(2:D)./(D-2+1);
        end
        
        case {'SMOP6'}
        [temp,g] = SMOPP6(a,theta);
        K = ceil(theta*(D-2+1));
        if j==1
            df(1) = (1+g/(D-2+1))*sin(a(1)*pi/2)*pi/2;
            df(2:D) = 2*(a(2:D)-pi/3)+...
                      72*pi*pi*linspace(0,1,D-2+1).*(a(2:D)-pi/3).*cos(power((6*pi*(a(2:D)-pi/3)),2));
            df(2:D) =  df(2:D).* (1-cos(a(1)*pi/2))./(D-2+1);
            df(temp) = 0;
        elseif j==2
            df(1) = -(1+g/(D-2+1))*cos(a(1)*pi/2)*pi/2;
            df(2:D) = 2*(a(2:D)-pi/3)+...
                      72*pi*pi*linspace(0,1,D-2+1).*(a(2:D)-pi/3).*cos(power((6*pi*(a(2:D)-pi/3)),2));
            df(2:D) = df(2:D).* (1-sin(a(1)*pi/2))./(D-2+1);
            df(temp) = 0;
        end
        
        case {'SMOP7'}
        [~,g] = SMOP7(a,theta);
        K = ceil(theta*(D-2+1));
        if j==1
            df(1) = -(1+g/(D-2+1))*sin(a(1)*pi/2)*pi/2;
            df(2:2+K-1) = cos(a(1)*pi/2)*4*(a(2:2+K-1)-pi/3+pi*sin(2*pi*(a(2:2+K-1)-pi/3)).*cos(2*pi*(a(2:2+K-1)-pi/3)));
            %少乘4
            for i = 2+K+1:D-1
                df(i) = a(i)-0.9*a(i+1)+pi*sin(2*pi*(a(i)-0.9*a(i+1)))*cos(2*pi*(a(i)-0.9*a(i+1)))...
                    -0.9*(a(i-1)-0.9*a(i)+pi*sin(2*pi*(a(i-1)-0.9*a(i)))*cos(2*pi*(a(i-1)-0.9*a(i))));
            end
            df(D) = a(D)-0.9*a(2+K)+pi*sin(2*pi*(a(D)-0.9*a(2+K)))*cos(2*pi*(a(D)-0.9*a(2+K)))...
                -0.9*(a(D-1)-0.9*a(D)+pi*sin(2*pi*(a(D-1)-0.9*a(D)))*cos(2*pi*(a(D-1)-0.9*a(D))));
            df(2+K) = a(2+K)-0.9*a(2+K+1)+pi*sin(2*pi*(a(2+K)-0.9*a(2+K+1)))*cos(2*pi*(a(2+K)-0.9*a(2+K+1)))...
                -0.9*(a(D)-0.9*a(2+K)+pi*sin(2*pi*(a(D)-0.9*a(2+K)))*cos(2*pi*(a(D)-0.9*a(2+K))));
            %乘回来
            df(2+K:D) = 4*cos(a(1)*pi/2).*df(2+K:D);
            df(2:D) = df(2:D)./(D-2+1);
        elseif j==2
            df(1) = (1+g/(D-2+1))*cos(a(1)*pi/2)*pi/2;
            df(2:2+K-1) = sin(a(1)*pi/2)*4*(a(2:2+K-1)-pi/3+pi*sin(2*pi*(a(2:2+K-1)-pi/3)).*cos(2*pi*(a(2:2+K-1)-pi/3)));
            for i = 2+K+1:D-1
                df(i) = a(i)-0.9*a(i+1)+pi*sin(2*pi*(a(i)-0.9*a(i+1)))*cos(2*pi*(a(i)-0.9*a(i+1)))...
                    -0.9*(a(i-1)-0.9*a(i)+pi*sin(2*pi*(a(i-1)-0.9*a(i)))*cos(2*pi*(a(i-1)-0.9*a(i))));
            end
            df(D) = a(D)-0.9*a(2+K)+pi*sin(2*pi*(a(D)-0.9*a(2+K)))*cos(2*pi*(a(D)-0.9*a(2+K)))...
                -0.9*(a(D-1)-0.9*a(D)+pi*sin(2*pi*(a(D-1)-0.9*a(D)))*cos(2*pi*(a(D-1)-0.9*a(D))));
            df(2+K) = a(2+K)-0.9*a(2+K+1)+pi*sin(2*pi*(a(2+K)-0.9*a(2+K+1)))*cos(2*pi*(a(2+K)-0.9*a(2+K+1)))...
                -0.9*(a(D)-0.9*a(2+K)+pi*sin(2*pi*(a(D)-0.9*a(2+K)))*cos(2*pi*(a(D)-0.9*a(2+K))));
            df(2+K:D) = 4*sin(a(1)*pi/2).*df(2+K:D);
            df(2:D) = df(2:D)./(D-2+1);
        end
        
         case {'SMOP8'}
        [~,g] = SMOP8(a,theta);
        K = ceil(theta*(D-2+1));
        if j==1
            df(1) = -(1+g/(D-2+1))*sin(a(1)*pi/2)*pi/2;
            df(2) = -1+800.*(a(2)-a(3)).*exp(-100.*(a(2)-a(3)).^2);
            for i = 3:2+K-1
                df(i) = -1+800.*(a(i)-mod(a(i+1)+pi,2)).*exp(-100.*(a(i)-mod(a(i+1)+pi,2)).^2)...
                    + 1-800.*(a(i-1)-mod(a(i)+pi,2)).*exp(-100.*(a(i-1)-mod(a(i)+pi,2)).^2);
            end
            %df(3:2+K-1) = 4*(a(2:2+K-1)-pi/3+pi*sin(2*pi*(a(2:2+K-1)-pi/3)).*cos(2*pi*(a(2:2+K-1)-pi/3)))/(D-2+1);
            for i = 2+K+1:D-1
                df(i) = -1+800.*(a(i)-a(i+1)).*exp(-100.*(a(i)-a(i+1)).^2)...
                    + 0.9*(1-800.*(a(i-1)-a(i)).*exp(-100.*(a(i)-a(i-1)).^2));
            end
            df(D) = -1+800.*(a(D)-a(2+K)).*exp(-100.*(a(D)-a(2+K)).^2)...
                    + 0.9*(1-800.*(a(D-1)-a(D)).*exp(-100.*(a(D)-a(D-1)).^2));
                
            df(2+K) = -1+800.*(a(2+K)-a(2+K+1)).*exp(-100.*(a(2+K)-a(2+K+1)).^2)...
                    + 0.9*(1-800.*(a(D)-a(2+K)).*exp(-100.*(a(2+K)-a(D)).^2))...
                    + 1-800.*(a(2+K-1)-mod(a(2+K)+pi,2)).*exp(-100.*(a(2+K-1)-mod(a(2+K)+pi,2)).^2);
            %df(2+K) = a(2+K)-0.9*a(2+K+1)+pi*sin(2*pi*(a(2+K)-0.9*a(2+K+1)))*cos(2*pi*(a(2+K)-0.9*a(2+K+1)))...
                %-0.9*(a(D)-0.9*a(2+K)+pi*sin(2*pi*(a(D)-0.9*a(2+K)))*cos(2*pi*(a(D)-0.9*a(2+K))));
            %df(2+K:D) = 4*cos(a(1)*pi/2);
            df(2:D) = cos(a(1)*pi/2).*df(2:D)./(D-2+1);
        elseif j==2
            df(1) = (1+g/(D-2+1))*cos(a(1)*pi/2)*pi/2;
            df(2) = -1-800.*(a(2)-a(3)).*exp(100.*(a(2)-a(3)).^2)./power(exp(100*(a(2)-a(3)).^2),2);
           for i = 3:2+K-1
                df(i) = -1+800.*(a(i)-mod(a(i+1)+pi,2)).*exp(-100.*(a(i)-mod(a(i+1)+pi,2)).^2)...
                    + 1-800.*(a(i-1)-mod(a(i)+pi,2)).*exp(-100.*(a(i-1)-mod(a(i)+pi,2)).^2);
            end
            %df(3:2+K-1) = 4*(a(2:2+K-1)-pi/3+pi*sin(2*pi*(a(2:2+K-1)-pi/3)).*cos(2*pi*(a(2:2+K-1)-pi/3)))/(D-2+1);
            for i = 2+K+1:D-1
                df(i) = -1+800.*(a(i)-a(i+1)).*exp(-100.*(a(i)-a(i+1)).^2)...
                    + 0.9*(1-800.*(a(i-1)-a(i)).*exp(-100.*(a(i)-a(i-1)).^2));
            end
            df(D) = -1+800.*(a(D)-a(2+K)).*exp(-100.*(a(D)-a(2+K)).^2)...
                    + 0.9*(1-800.*(a(D-1)-a(D)).*exp(-100.*(a(D)-a(D-1)).^2));
                
            df(2+K) = -1+800.*(a(2+K)-a(2+K+1)).*exp(-100.*(a(2+K)-a(2+K+1)).^2)...
                    + 0.9*(1-800.*(a(D)-a(2+K)).*exp(-100.*(a(2+K)-a(D)).^2))...
                    + 1-800.*(a(2+K-1)-mod(a(2+K)+pi,2)).*exp(-100.*(a(2+K-1)-mod(a(2+K)+pi,2)).^2);
            %df(2+K) = a(2+K)-0.9*a(2+K+1)+pi*sin(2*pi*(a(2+K)-0.9*a(2+K+1)))*cos(2*pi*(a(2+K)-0.9*a(2+K+1)))...
                %-0.9*(a(D)-0.9*a(2+K)+pi*sin(2*pi*(a(D)-0.9*a(2+K)))*cos(2*pi*(a(D)-0.9*a(2+K))));
            %df(2+K:D) = 4*cos(a(1)*pi/2);
            df(2:D) = sin(a(1)*pi/2).*df(2:D)./(D-2+1);
        end
end
end

function g = g1(x,t)
    g = (x-t).^2;
end