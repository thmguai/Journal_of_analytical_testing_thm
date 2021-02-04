function [newench2cross]=Cross(newench2select,Pc)
[Ordersj,Indexsj]=sort(rand(size(newench2select,1),1));%随机交叉
newench2select=newench2select(Indexsj,:);
lchrom=size(newench2select,2);
poscut=ceil(rand(size(newench2select,1)/2,1)*(lchrom-1));%选好随机交叉点
poscut=poscut.*(rand(size(poscut))<Pc); %根据交叉概率进行交叉
for i=1:length(poscut),%length（A）表示A的最大行/列数
     newench2cross([2*i-1 2*i],:)=[newench2select([2*i-1 2*i],1:poscut(i)) newench2select([2*i 2*i-1],poscut(i)+1:lchrom)];%单点交叉
end
