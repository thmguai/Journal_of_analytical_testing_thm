function [Ag,Real]=affinity1(x,popsize,X,len) 
for i=1:popsize
    x1=x(i,:);
    Real(i) = fitness(x1,X,len);
    Ag(i)=Real(i)/1;%在求取极大值时，两者相等；求极小值时，转为求最大值问题，并归一化
end
