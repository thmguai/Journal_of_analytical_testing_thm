%clear all;
clc;

popsize=50;% the value of population ��Ⱥȡ����
%CodeL=2;%�������

G=100;  % the max generation
mnum=8; % mnum---����ϸ����Ŀ
Pc=0.85;
Pm=0.05;
Tacl=0.85;     %���㿹��Ũ��

X=xlsread('�ļ�λ��');%��������
Y = X(:,1);
X2=xlsread('�ļ�λ��');%��������
Y1 = X2(:,1);
[m,n] = size(X);
len=ceil((n-1)/4);
K1 = 50;
T1 = 0;
T2 = 0;
T3 = 0;
T4 = 0;
for k1 = 1:K1
%---------------------��ʼ����Ⱥ-------------------------%
for i=1:popsize
    x(i,:)=Initial(len);
    fx(i)=fitness(x(i,:),X,len);
end

%--------------------------------------------------------%
[Ag,RealValue]=affinity1(x,popsize,X,len);%�����뿹ԭ���׺���
%-------------------������ʼ����ϸ��----------------------%
for i=1:mnum 
    mcell(i,:)=Initial(len);
    mfx(i)=fitness(mcell(i,:),X,len);
end
%--------------------------------------------------------%
[mx,mAg]=memorycell(x,mcell,mnum,popsize,X,len);%���³�ʼ����ϸ��
[Order_Ag,Index_Ag]=sort(Ag);%��ʼ����ϸ�����³�ʼ��Ⱥ
x(Index_Ag(1:mnum),:)=mx;
[Ag,RealValue]=affinity1(x,popsize,X,len);%�����뿹ԭ���׺���
[Ab]=affinity2(x,popsize,len);%�����뿹��֮����׺���

%----------------------------�Ľ��ĵط�----------------------------%
Ab_min = min(min(Ab));
Ab_max = max(max(Ab));
C_std_max = 0;
C_max = zeros(1,popsize);
for Tacl = Ab_min:0.001:Ab_max
    C=(sum(Ab>Tacl,2)/popsize)';%���㿹��Ũ��
    C_temp = std(C);
    if(C_temp > C_std_max)
        Tacl_great = Tacl;
        C_std_max = C_temp;
        C_max = C;
    end
end
C = C_max;
%C=(sum(Ab>Tacl,2)/popsize)';%���㿹��Ũ��
[Best_Ag,Index]=max(Ag);%���ʼ����Ӣ����
Best_Value(1)=RealValue(Index);
%Temp_Value=Best_Value(1);
Best_gene=x(Index,:);
% [Worst_Ag,Index1]=min(Ag);%���ʼ������
% Worst_Value(1)=RealValue(Index1);
% Worst_gene=x(Index1,:);
avg_Ag=0;
for k=1:popsize
    avg_Ag=avg_Ag+RealValue(k);
end
avg_Ag=avg_Ag/popsize;%���ʼ�������뿹ԭ�׺Ͷ�ƽ��ֵ
fprintf(1,'\n1--->Best : %f Avg : %f \n',Best_Value(1),avg_Ag)
[mx,mAg]=memorycell(x,mcell,mnum,popsize,X,len);%���³�ʼ����ϸ��
%-------------------------------��ʼ����----------------------------------%
time(1)=1;
m_mag(1) = mean(mAg);
m_ag(1) = mean(Ag);
arr_best_v(1) = Best_Ag;
arr_best_ag(1) = Best_Ag;
for gen=2:G
  time(gen)=gen;  
%---------------------1:Select Cross Mutation Operation----------------%
[E]=Select_Antibody(C,Ag,popsize,gen,G);
[newench2select]=Select(x,E,popsize);%���Ƹ���
[newench2cross]=Cross(newench2select,Pc);%����
[ench2]=Mutate(newench2cross,Pm);%����
%[chrom]=decode(ench2,len,MinX,MaxX);%����
%-----------------------------2:���߲���-------------------------------%
[Ag,RealValue]=affinity1(ench2,popsize,X,len);%�����뿹ԭ���׺���
[Ab]=affinity2(ench2,popsize,len);%�����뿹��֮����׺���

%----------------------------�Ľ��ĵط�----------------------------%
Ab_min = min(min(Ab));
Ab_max = max(max(Ab));
C_std_max = 0;
C_max = zeros(1,popsize);
arr_i = 1;
for Tacl = Ab_min:0.001:Ab_max
    C=(sum(Ab>Tacl,2)/popsize)';%���㿹��Ũ��
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

%C=(sum(Ab>Tacl,2)/popsize)';%���㿹��Ũ��
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
%̰�Ĳ���

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
%��Ӣ���屣������
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
%     fprintf(1,'�ڵ�%d�����������Ž� ',gen);
%     return;
% end
[mx,mAg]=memorycell(x,mcell,mnum,popsize,X,len);%���³�ʼ����ϸ��
[~,Index_Ag]=sort(Ag);%���±����Ļ���
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
%------------------------------��ͼ--------------------------------%
%save IMGA_MAX_F1 time Best_Value
figure(3);
plot(time,Best_Value,'k');
xlabel('��������','Fontname','Times New Roman','FontSize',13);
ylabel('�׺Ͷ�','Fontname','Times New Roman','FontSize',13);% the value of objective
box off
set(gca, 'Fontname', 'Times New Roman','FontSize',13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%�����Ŵ��㷨�Ӻ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
%---------------------���㿹���뿹ԭ���׺Ͷ�----------------%

    
