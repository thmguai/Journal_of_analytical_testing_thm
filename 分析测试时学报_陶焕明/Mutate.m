function [newchmutate]=Mutate(newchcross,Pm)
point=find(rand(size(newchcross))<Pm);%���ֱ����
newchmutate=newchcross;
newchmutate(point)=1-newchcross(point);   