function [Ag,Real]=affinity1(x,popsize,X,len) 
for i=1:popsize
    x1=x(i,:);
    Real(i) = fitness(x1,X,len);
    Ag(i)=Real(i)/1;%����ȡ����ֵʱ��������ȣ���Сֵʱ��תΪ�����ֵ���⣬����һ��
end
