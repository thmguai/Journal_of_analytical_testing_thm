function [mx,mAg]=memorycell(x,mcell,mnum,popsize,X,len)
snum=mnum+popsize;
for i=1:popsize
    x1=x(i,:);
    GxAg(i)=fitness(x1,X,len);
end
for i=1:mnum
    GxmAg(i)=fitness(mcell(i,:),X,len);
end
TAg(1:popsize)=GxAg;
Tch(1:popsize,:)=x;
TAg(popsize+1:popsize+mnum)=GxmAg;
Tch(popsize+1:popsize+mnum,:)=mcell;
[OrderTAg,IndexTAg]=sort(TAg);
for i=1:mnum
    mx(i,:)=Tch(IndexTAg(snum+1-i),:);
    mAg(i)=OrderTAg(snum+1-i);
end

