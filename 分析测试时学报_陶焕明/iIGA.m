%clear all;
clc;

popsize=50;% the value of population 种群取复数
%CodeL=2;%基因个数

G=100;  % the max generation
mnum=8; % mnum---记忆细胞数目
Pc=0.85;
Pm=0.05;
Tacl=0.85;     %计算抗体浓度

X=xlsread('文件位置');%导入数据
Y = X(:,1);
X2=xlsread('文件位置');%导入数据
Y1 = X2(:,1);
[m,n] = size(X);
len=ceil((n-1)/4);
K1 = 50;
T1 = 0;
T2 = 0;
T3 = 0;
T4 = 0;
for k1 = 1:K1
%---------------------初始化种群-------------------------%
for i=1:popsize
    x(i,:)=Initial(len);
    fx(i)=fitness(x(i,:),X,len);
end

%--------------------------------------------------------%
[Ag,RealValue]=affinity1(x,popsize,X,len);%抗体与抗原的亲和力
%-------------------产生初始记忆细胞----------------------%
for i=1:mnum 
    mcell(i,:)=Initial(len);
    mfx(i)=fitness(mcell(i,:),X,len);
end
%--------------------------------------------------------%
[mx,mAg]=memorycell(x,mcell,mnum,popsize,X,len);%更新初始记忆细胞
[Order_Ag,Index_Ag]=sort(Ag);%初始记忆细胞更新初始种群
x(Index_Ag(1:mnum),:)=mx;
[Ag,RealValue]=affinity1(x,popsize,X,len);%抗体与抗原的亲和力
[Ab]=affinity2(x,popsize,len);%抗体与抗体之间的亲和力

%----------------------------改进的地方----------------------------%
Ab_min = min(min(Ab));
Ab_max = max(max(Ab));
C_std_max = 0;
C_max = zeros(1,popsize);
for Tacl = Ab_min:0.001:Ab_max
    C=(sum(Ab>Tacl,2)/popsize)';%计算抗体浓度
    C_temp = std(C);
    if(C_temp > C_std_max)
        Tacl_great = Tacl;
        C_std_max = C_temp;
        C_max = C;
    end
end
C = C_max;
%C=(sum(Ab>Tacl,2)/popsize)';%计算抗体浓度
[Best_Ag,Index]=max(Ag);%求初始代精英抗体
Best_Value(1)=RealValue(Index);
%Temp_Value=Best_Value(1);
Best_gene=x(Index,:);
% [Worst_Ag,Index1]=min(Ag);%求初始代最差抗体
% Worst_Value(1)=RealValue(Index1);
% Worst_gene=x(Index1,:);
avg_Ag=0;
for k=1:popsize
    avg_Ag=avg_Ag+RealValue(k);
end
avg_Ag=avg_Ag/popsize;%求初始代抗体与抗原亲和度平均值
fprintf(1,'\n1--->Best : %f Avg : %f \n',Best_Value(1),avg_Ag)
[mx,mAg]=memorycell(x,mcell,mnum,popsize,X,len);%更新初始记忆细胞
%-------------------------------开始进化----------------------------------%
time(1)=1;
m_mag(1) = mean(mAg);
m_ag(1) = mean(Ag);
arr_best_v(1) = Best_Ag;
arr_best_ag(1) = Best_Ag;
for gen=2:G
  time(gen)=gen;  
%---------------------1:Select Cross Mutation Operation----------------%
[E]=Select_Antibody(C,Ag,popsize,gen,G);
[newench2select]=Select(x,E,popsize);%复制个体
[newench2cross]=Cross(newench2select,Pc);%交叉
[ench2]=Mutate(newench2cross,Pm);%变异
%[chrom]=decode(ench2,len,MinX,MaxX);%解码
%-----------------------------2:免疫操作-------------------------------%
[Ag,RealValue]=affinity1(ench2,popsize,X,len);%抗体与抗原的亲和力
[Ab]=affinity2(ench2,popsize,len);%抗体与抗体之间的亲和力

%----------------------------改进的地方----------------------------%
Ab_min = min(min(Ab));
Ab_max = max(max(Ab));
C_std_max = 0;
C_max = zeros(1,popsize);
arr_i = 1;
for Tacl = Ab_min:0.001:Ab_max
    C=(sum(Ab>Tacl,2)/popsize)';%计算抗体浓度
    C_temp = std(C);
    arr_c_temp(arr_i) = C_temp;
    if(C_temp > C_std_max)
        Tacl_great = Tacl;
        C_std_max = C_temp;
        C_max = C;
    end
    arr_i = arr_i + 1;
end
C = C_max;

% plot(Ab_min:0.001:Ab_max,arr_c_temp);
% xlim([Ab_min,Ab_max])
% ylim([0,0.1])
% xlabel('Similarity threshold','Fontname','Times New Roman','FontSize',13)
% ylabel('STD','Fontname','Times New Roman','FontSize',13)
% box off
% set(gca, 'Fontname', 'Times New Roman','FontSize',13);

%C=(sum(Ab>Tacl,2)/popsize)';%计算抗体浓度
[Best_cur_Ag,Index]=max(Ag);
Best_cur_gene=ench2(Index,:);
[Worst_cur_Ag,index]=min(Ag);
Worst_Value=RealValue(index);
Worst_cur_gene=ench2(index);
avg_Value=0;
for k=1:popsize
    avg_Value=avg_Value+RealValue(k);
end
avg_Value=avg_Value/popsize;
%贪心策略

%Best_cur_gene1 = Best_cur_gene;
for i = 1:20
    random = round(175*rand())+1;
    if random == 176
        random = 175;
    end
    if Best_cur_gene(random) == 0
        Best_cur_gene(random) = 1;
    elseif Best_cur_gene(random) == 1
        Best_cur_gene(random) = 0;
    end
    best_v=fitness(Best_cur_gene,X,len);
    if best_v > Best_cur_Ag
        arr_ii(gen) = i;
        arr_best_v(gen) = best_v;
        arr_best_ag(gen) = Best_cur_Ag;
        Best_cur_Ag = best_v;
        break;
    elseif best_v <= Best_cur_Ag
        if Best_cur_gene(random) == 0
            Best_cur_gene(random) = 1;
        elseif Best_cur_gene(random) == 1
            Best_cur_gene(random) = 0;
        end
    end
    if i == 20
        arr_best_v(gen) = Best_cur_Ag;
        arr_best_ag(gen) = Best_cur_Ag;
    end
        
end
ench2(index,:)=Best_cur_gene;
%精英抗体保留策略
if Best_cur_Ag < Best_Ag
    Best_Value(gen)=Best_Ag;
    %ench2(index,:)=Best_gene;
elseif Best_cur_Ag >= Best_Ag
    Best_Value(gen)=Best_cur_Ag;
    Best_Ag=Best_cur_Ag;
    %Temp_Value=RealValue(Index);
    Best_gene=Best_cur_gene;
end
% fprintf(1,'\n%d--->BEST : %f AVG : %f \n',gen,Best_Value(gen),avg_Value);
% fprintf(1,'--->Value_F = %f',Best_Value(gen));
disp(Best_gene)
% if abs(Best_Value(gen)-Best_result)<0.001
%     fprintf(1,'在第%d代收敛到最优解 ',gen);
%     return;
% end
[mx,mAg]=memorycell(x,mcell,mnum,popsize,X,len);%更新初始记忆细胞
[~,Index_Ag]=sort(Ag);%更新本代的基因
ench2(Index_Ag(1:mnum),:)=mx;
x=ench2;
Ag(Index_Ag(1:mnum))=mAg;

m_mag(gen) = mean(mAg);
m_ag(gen) = mean(Ag);
end

bar(arr_best_v);
hold on; 
bar(arr_best_ag,'w');
xlim([1,100]);
ylim([0.7,0.8]);
xlabel('Number of iterations','Fontname','Times New Roman','FontSize',13)
ylabel('The highest affinity','Fontname','Times New Roman','FontSize',13)
box off
set(gca, 'Fontname', 'Times New Roman','FontSize',13);
h = legend('After explore','Before explore','location','NorthWest');
set(h,'Fontname','Times New Roman','fontsize',9,'box','off');


% plot(1:100,m_mag,'-r',1:100,m_ag,'-g');
% ylim([0.68,0.8])
% xlabel('Number of iterations','Fontname','Times New Roman','FontSize',13)
% ylabel('Average affinity','Fontname','Times New Roman','FontSize',13)
% box off
% set(gca, 'Fontname', 'Times New Roman','FontSize',13);
% h = legend('Elite group','General population','location','NorthWest');
% set(h,'Fontname','Times New Roman','fontsize',9,'box','off');


[Rc,RMSEC,beta,yc]= fitaaa(Best_gene,X,len);
j=1;
[n,m]=size(X);
for i=1:len
    
    if Best_gene(i) == 1 && i == len
            X_end(:,j:j+m-(4*i-2)) = X(:,4*i-2:end);
    elseif Best_gene(i) == 1
            X_end(:,j:j+3) = X(:,4*i-2:4*i+1);
            j = j+4;
    end
    
end
% X_end_1 = UVE([Y,X_end]);
% [Rc1,RMSEC1,beta1,yc1] = fitaaa_uve(X_end_1);

% figure(1);
% plot(Y,yc,'ro');
% xlabel('Actual Value');
% ylabel('Predictive Value');
% text(7.7,9.8,['Rc = ',num2str(Rc)]);
% text(7.7,9.65,['RMSEC = ',num2str(RMSEC)]);
% hold on;
% plot([7.5,10],[7.5,10],'linewidth',1.5);
[Rp,RMSEP,yp]= fitbbb(Best_gene,X2,len,beta);
Rp3(k1) = RMSEP;
T1 = T1 + Rc;
T2 = T2 + RMSEC;
T3 = T3 + Rp;
T4 = T4 + RMSEP;
% figure(2);
% plot(Y1,yp,'ro');
% xlabel('Actual Value');
% ylabel('Predictive Value');
% text(7.7,9.8,['Rp = ',num2str(Rp)]);
% text(7.7,9.65,['RMSEP = ',num2str(RMSEP)]);
% hold on;
% plot([7.5,10],[7.5,10],'linewidth',1.5);
end
disp(T1/K1);
disp(T2/K1);
disp(T3/K1);
disp(T4/K1);
%------------------------------画图--------------------------------%
%save IMGA_MAX_F1 time Best_Value
figure(3);
plot(time,Best_Value,'k');
xlabel('迭代次数','Fontname','Times New Roman','FontSize',13);
ylabel('亲和度','Fontname','Times New Roman','FontSize',13);% the value of objective
box off
set(gca, 'Fontname', 'Times New Roman','FontSize',13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%免疫遗传算法子函数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
%---------------------计算抗体与抗原的亲和度----------------%

    
