%% -----------------------------
%% �ļ����ƣ����㷨���� SMOP
%% �޸����ڣ�2020.04.15
%% �����޸ģ�������
%% -----------------------------
clc
clear
close all;
%% ���������
rng('default'); 
for index_problem = [1 2 3 4 5 6 7 8]
%% ��������
w_initial                         = 0.7968;       % ����Ȩ��0.8
c1_initial                        = 1.4962;       % ���¹�ʽ�в���
c2_initial                        = 1.4962;       % ���¹�ʽ�в���
fi                                = 0;            % CSO��loser���¹�ʽ0.2
loop                              = 30;           % ��������30
NA                                = 100;          % ����A����󴢴�������
Pop_size                          = 100;          % ��Ⱥ��ģ 100
Max_iter                          = 5000;         % һ�ε����������´���5000����1000ά������ԣ�
rec                               = [0 0 800 620];
theta                             = 0.05;         % ϡ��̶�
theta_w                           = '05';
D                                 = 1000;         % ����ά��
how_update                        = 'CCPSO' ;     % CCPSO CSO PSO
set_L_max                         = 0.05;         % ռ������ı���
L_min                             = 0;            % ռ������ı���
set_eta                           = 0;  % �ɺ���
%% ����
Problem = {'SMOP1','SMOP2','SMOP3','SMOP4','SMOP5','SMOP6','SMOP7','SMOP8'};
funcname = Problem{index_problem};
%% ������front set
switch funcname
    case {'SMOP1'}
        PF = UniformPoint(500,2);%500
    case {'SMOP2'}
        PF = UniformPoint(500,2);
    case {'SMOP3'}
        PF = UniformPoint(500,2);
    case {'SMOP4'}
        P = UniformPoint(500,2);
        c = ones(size(P,1),2);
        for i = 1 : size(P,1)
            for j = 2 : 2
                temp = P(i,j)/P(i,1)*prod(1-c(i,2-j+2:2-1));
                c(i,2-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
            end
        end
        x = acos(c)*2/pi;
        PF = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
        x=[];
        P=[];
    case {'SMOP5'}
        P = UniformPoint(500,2);
        c = ones(size(P,1),2);
        for i = 1 : size(P,1)
            for j = 2 : 2
                temp = P(i,j)/P(i,1)*prod(1-c(i,2-j+2:2-1));
                c(i,2-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
            end
        end
        x = acos(c)*2/pi;
        PF = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
        x=[];
        P=[];
    case {'SMOP6'}
        P = UniformPoint(500,2);
        c = ones(size(P,1),2);
        for i = 1 : size(P,1)
            for j = 2 : 2
                temp = P(i,j)/P(i,1)*prod(1-c(i,2-j+2:2-1));
                c(i,2-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
            end
        end
        x = acos(c)*2/pi;
        PF = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
        x=[];
        P=[];
    case {'SMOP7'}
        PF = UniformPoint(500,2);
        PF = PF./repmat(sqrt(sum(PF.^2,2)),1,2);
         case {'SMOP8'}
        PF = UniformPoint(500,2);
        PF = PF./repmat(sqrt(sum(PF.^2,2)),1,2);
end
%% ������·������
rep = '.\\results';
if not(exist(rep))
    mkdir(rep);
end

rep = sprintf('./results/%s',funcname);
if not(exist(rep))
    mkdir(rep);
end

 for L_max = set_L_max
    str = sprintf('lambda_%04d',floor(10000*L_max));
    resultsSaveDir1 = [rep '\\' str];
    if not(exist(resultsSaveDir1))
        mkdir(resultsSaveDir1);
    end
    for eta = set_eta
        str = sprintf('eta_%04d',floor(10000* eta));
        resultsSaveDir2 = [resultsSaveDir1 '\\' str];
        if not(exist(resultsSaveDir2))
            mkdir(resultsSaveDir2);
        end
                
                resultsSaveFigDir = [resultsSaveDir2 '\\fig'];
                if not(exist(resultsSaveFigDir))
                    mkdir(resultsSaveFigDir);
                end
                
                resultsSaveEmfDir = [resultsSaveDir2  '\\emf'];
                if not(exist(resultsSaveEmfDir))
                    mkdir(resultsSaveEmfDir);
                end
                
                resultsSaveMatDir = [resultsSaveDir2  '\\mat'];
                if not(exist(resultsSaveMatDir))
                    mkdir(resultsSaveMatDir);
                end
                
                %% ��ʼѭ������
                IGD = [];
                SP = [];
                A_record = [];
                value_A_record = [];
               
                
                for iloop = 1:loop %��loop�ζ�����ʵ��
                    iloop;
                    
                    rng(iloop);
                    %�����ڴ�
                    A =[];
                    value_A = [];
                    A_Front = [];
                    value_A_Front = [];
                    pBest = [];
                    value_pBest=[];
                    iter = 1;
                    gBest_rec =[];
                    gBest = [];
                    %% ��Ⱥ��ʼ��
                    w  = w_initial;
                    c1 = c1_initial;
                    c2 = c2_initial;
                    
                            func_dim = D;
                             xrange(1,1) = 0;
                             xrange(1,2) = 1;
                             xrange(2:D,1) = ones(1,D-1)*(-1);
                             xrange(2:D,2) = ones(1,D-1)*2;

                      step = (L_max-L_min).*1/(Max_iter-1);
                      lamb(1,:) = 1*L_max : -step :1*L_min;
                      step = (L_max-L_min).*3/(Max_iter-1);
                      lamb(2,:) = 3*L_max : -step :3*L_min;
                             
                  
                    x = rand(Pop_size,func_dim);
                    v = rand(Pop_size,func_dim);
                    for k=1:func_dim
                        x(:,k) = x(:,k)*(xrange(k,2)-xrange(k,1))+xrange(k,1);
                        v(:,k) = v(:,k)*(xrange(k,2)-xrange(k,1))+xrange(k,1);
                    end
                  
                    %���Ŀ�꺯��ֵ
                    switch funcname
                        case {'SMOP1'}
                            [Obj_fitness,~] = SMOP1(x,theta);
                        case {'SMOP2'}
                            [Obj_fitness,~] = SMOP2(x,theta);
                        case {'SMOP3'}
                            [Obj_fitness,~] = SMOP3(x,theta);
                        case {'SMOP4'}
                            [Obj_fitness,~,~] = SMOP4(x,theta);
                        case {'SMOP5'}
                            [Obj_fitness,~] = SMOP5(x,theta);
                        case {'SMOP6'}
                            [Obj_fitness,~] = SMOP6(x,theta);
                        case {'SMOP7'}
                            [Obj_fitness,~] = SMOP7(x,theta);
                            case {'SMOP8'}
                            [Obj_fitness,~] = SMOP8(x,theta);
                    end                                  
                    

                    
                    %% ��һ�ε�������ʼֵ��ΪpBest
                    pBest = x; %��һ��ֱ�Ӵ���
                    value_pBest = Obj_fitness;
                    g = [];
                    [x,Obj_fitness,g] = RDA_x(x,Obj_fitness,xrange,funcname,theta,lamb(:,iter),eta,g,iter);
                    
                   
                    [x,Obj_fitness,pBest, value_pBest] = pBest_get(x,Obj_fitness,pBest,value_pBest);
                    %% ���A
                
                    [A,value_A] = A_get(pBest,value_pBest,A,value_A,NA,xrange,funcname,theta);

%                         plot(value_A(:,1),value_A(:,2),'r.')
%                         hold on
%                         pause(0.1)
%                         clf
                   
                    igdigd(1) = get1_IGD(value_A,PF);
                    %��pBest���ҵ�gBest
                    [gBest] = gBest_get(A,value_A);
               
                    gBest_rec = [gBest_rec;gBest];
                    
                    switch how_update
                        case {'CCPSO'}
                               %����x �� v,�ں����м���lbest��ƽ��ֵ
                               [x,v] = CSO(x,v,w,c1,c2,fi,pBest,gBest,func_dim,xrange,Obj_fitness);
                        case {'CSO'}
                            [x,v] = orignalCSO(x,v,fi,Obj_fitness,xrange);
                        case {'PSO'}
                            [x,v] =  pso(x,v,w,c1,c2,pBest,gBest,func_dim,xrange,Obj_fitness);
                    end
                    
                    %�ӵڶ�����ʼ����ѭ��
                    for iter = 2:Max_iter
                        %iter
                        %% ���Ŀ�꺯��ֵ
                        switch funcname
                            case {'SMOP1'}
                                [Obj_fitness,~] = SMOP1(x,theta);
                            case {'SMOP2'}
                                [Obj_fitness,~] = SMOP2(x,theta);
                            case {'SMOP3'}
                                [Obj_fitness,~] = SMOP3(x,theta);
                            case {'SMOP4'}
                                [Obj_fitness,~,~] = SMOP4(x,theta);
                            case {'SMOP5'}
                                [Obj_fitness,~] = SMOP5(x,theta);
                                case {'SMOP6'}
                            [Obj_fitness,~] = SMOP6(x,theta);
                            case {'SMOP7'}
                                [Obj_fitness,~] = SMOP7(x,theta);
                                case {'SMOP8'}
                                [Obj_fitness,~] = SMOP8(x,theta);
                        end
                       
                        %% ��RDA���Ӹ���x,���ɺ��x��2*Pop_size����Ⱥ��С������ѡȡpBest
                        [x,Obj_fitness,g] = RDA_x(x,Obj_fitness,xrange,funcname,theta,lamb(:,iter),eta,g,iter);
                        %% �Ƚϣ��ҵ�����ֵ�⣬����pBest(�ɲο�pso�����е���������gBest,��A
                        %����pBest
                        [x,Obj_fitness,pBest, value_pBest] = pBest_get(x,Obj_fitness,pBest,value_pBest);
                        
                        %% ���£���֦��gBest
                        % ����A % ��֦
                        [A,value_A] = A_get(pBest,value_pBest,A,value_A,NA,xrange,funcname,theta);
                  
                        %% ����
%                             plot(value_A(:,1),value_A(:,2),'r.')
%                             hold on
%                             pause(0.1)
%                             clf
                        %% ����igd��ֵ
                        igdigd(iter) = get1_IGD(value_A,PF);
                        igdigd(iter);
                        %% ��gBest
                        [gBest] = gBest_get(A,value_A);
                        gBest_rec = [gBest_rec;gBest];
                        
                        %% ����x v CSO 
                        switch how_update
                        case {'CCPSO'}
                               %����x �� v,�ں����м���lbest��ƽ��ֵ
                               [x,v] = CSO(x,v,w,c1,c2,fi,pBest,gBest,func_dim,xrange,Obj_fitness);
                        case {'CSO'}
                            [x,v] = orignalCSO(x,v,fi,Obj_fitness,xrange);
                        case {'PSO'}
                            [x,v] =  pso(x,v,w,c1,c2,pBest,gBest,func_dim,xrange,Obj_fitness);
                    end
                       
                    end
                    %% A(t)������
                    
                    A_record = [A_record;A];  %һ�ε���������ĵ���A��������
                    value_A_record = [value_A_record;value_A];
                    A_index(iloop) = size(A,1);
                
                    %% ���۱�׼ Performance Metrics
                    
                    igd1 = get1_IGD(value_A,PF);
                    IGD(iloop) = igd1;
                   
                    
                    %Spacing (SP)
                    [sp] = get_SP(value_A);
                    SP(iloop) = sp;
                    
                    
                    
                end
                %% IGD���۱�׼���ۺ�
                %�����ֵ��������ֵ������ֵ����¼
                IGD_best = min(IGD);
                IGD_worst = max(IGD);
                IGD_mean = mean(IGD);
                IGD_std = std(IGD);
                disp('Done');
                
                SP_best = min(SP);
                SP_worst = max(SP);
                SP_mean = mean(SP);
                SP_std  = std(SP);
                disp('Done');
                
                
                %% ��ͼ
                %����
                str1 = sprintf(' theta: %2f, update method: %s, funcname: %s\n', theta, how_update, funcname);
                str2 = sprintf('IGD best:%.6f,  IGD worst:%.6f,  IGD mean:%.6f,  IGD std:%.7f\n',...
                    IGD_best,IGD_worst,IGD_mean,IGD_std);
                str3 = sprintf('SP best:%.6f,  SP worst:%.6f,  SP mean:%.6f,  SP std:%.7f\n',...
                    SP_best,SP_worst,SP_mean,SP_std);
                str4 = sprintf('maximum size of A: %1f, Population size: %1f\n',...
                    NA,Pop_size);
                str6 = sprintf('fi:%2f,L_max:%3f,L_min:%3f,loop:%2f\n',...
                    fi,L_max,L_min,loop);
                str5 = sprintf('w: %2f, c1: %4f, c2: %4f \n',...
                    w,c1,c2);
                title_graph = [str1 str4 str5 str6 str2 str3];
                disp('----------------------------');
                disp(title_graph);
                %% ��ͼ
                %% ��һ��ͼ������ѭ�����������ϵ�ͼ
                filename = sprintf('Loop_%d_Update_%s_theta_%s',...
                      loop,how_update,theta_w);
                
                h1 = figure('Name','all_loop','Position',rec);
                hold on
                plot(value_A_record(:,1),value_A_record(:,2),'r.');
                %plot(A_final(:,1),value_A_final(:,2),'r.');
                hold on
                plot(PF(:,1),PF(:,2));
                
                xlabel('f1')
                ylabel('f2')
                box on;
                title(title_graph);
                saveas(h1, [resultsSaveFigDir '\\' filename '_all_loop_record'], 'fig');
                close(figure(h1));
                
%                  %% �ڶ���ͼ��Front
%                filename = sprintf('Loop_%d_Update_%s_theta_%s',...
%                       loop,how_update,theta_w);
%                 
%                 h1 = figure('Name','all_loop','Position',rec);
%                 hold on
%                 plot(value_A_Front(:,1),value_A_Front(:,2),'r.');
%                 %plot(A_final(:,1),value_A_final(:,2),'r.');
%                 hold on
%                 plot(PF(:,1),PF(:,2));
%                 
%                 xlabel('f1')
%                 ylabel('f2')
%                 box on;
%                 title(title_graph);
%                 saveas(h1, [resultsSaveFigDir '\\' filename '_all_loop_Front'], 'fig');
%                 close(figure(h1));
                
                                              %% ������ͼ Median
                 filename = sprintf('Loop_%d_Update_%s_theta_%s',...
                      loop,how_update,theta_w);
                
                  % median
                  [~,value_index] = sort(IGD);
                  value_index_m = value_index(ceil(loop/2));
                  
                  value_A_index1 = sum(A_index(1:value_index_m-1))+1;
                  value_A_index2 = sum(A_index(1:value_index_m));
                  
                h2 = figure('Name','Median','Position',rec);
                hold on
                plot(value_A_record(value_A_index1:value_A_index2,1),value_A_record(value_A_index1:value_A_index2,2),'r.');
                hold on
                plot(PF(:,1),PF(:,2));
                value_A_median = value_A_record(value_A_index1:value_A_index2,:);
                xlabel('f1')
                ylabel('f2')
                box on;
                title(title_graph);
                saveas(h2, [resultsSaveFigDir '\\' filename '_median_loop'], 'fig');
                close(figure(h2));
                
                %% ���ķ�ͼ igd
                filename = sprintf('Loop_%d_Update_%s_theta_%s',...
                      loop,how_update,theta_w);
                
                h2 = figure('Name','Train MSE','Position',rec);
                hold on
                plot(igdigd);
                
                xlabel('time')
                ylabel('IGD')
                box on;
                title(title_graph);
                saveas(h2, [resultsSaveFigDir '\\' filename '_IGD'], 'fig');
                close(figure(h2));
                %% ��������
                %IGD
                save(fullfile([resultsSaveMatDir '\\' filename '_IGD_best.mat']),'IGD_best');
                save(fullfile([resultsSaveMatDir '\\' filename '_IGD_worst.mat']),'IGD_worst');
                save(fullfile([resultsSaveMatDir '\\' filename '_IGD_mean.mat']),'IGD_mean');
                save(fullfile([resultsSaveMatDir '\\' filename '_IGD_std.mat']),'IGD_std');
                save(fullfile([resultsSaveMatDir '\\' filename '_IGD.mat']),'IGD');
               save(fullfile([resultsSaveMatDir '\\' filename '_IGD.txt']),'IGD','-ascii');
                
                save(fullfile([resultsSaveMatDir '\\' filename '_results.txt']),'title_graph','-ascii');
                %SP
                save(fullfile([resultsSaveMatDir '\\' filename '_SP_best.mat']),'SP_best');
                save(fullfile([resultsSaveMatDir '\\' filename '_SP_worst.mat']),'SP_worst');
                save(fullfile([resultsSaveMatDir '\\' filename '_SP_mean.mat']),'SP_mean');
                save(fullfile([resultsSaveMatDir '\\' filename '_SP_std.mat']),'SP_std');
                save(fullfile([resultsSaveMatDir '\\' filename '_SP.mat']),'SP');
                
                %A
                save(fullfile([resultsSaveMatDir '\\' filename '_A_record.mat']),'A_record');
                save(fullfile([resultsSaveMatDir '\\' filename '_value_A_record.mat']),'value_A_record');
                save(fullfile([resultsSaveMatDir '\\' filename '_A_index.mat']),'A_index');
                save(fullfile([resultsSaveMatDir '\\' filename '_value_A_median.mat']),'value_A_median');
                disp('Done');
            
        
    end
    end

end