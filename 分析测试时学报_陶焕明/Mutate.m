function [newchmutate]=Mutate(newchcross,Pm)
point=find(rand(size(newchcross))<Pm);%发现变异点
newchmutate=newchcross;
newchmutate(point)=1-newchcross(point);   