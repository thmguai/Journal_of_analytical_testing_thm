function [newench2cross]=Cross(newench2select,Pc)
[Ordersj,Indexsj]=sort(rand(size(newench2select,1),1));%�������
newench2select=newench2select(Indexsj,:);
lchrom=size(newench2select,2);
poscut=ceil(rand(size(newench2select,1)/2,1)*(lchrom-1));%ѡ����������
poscut=poscut.*(rand(size(poscut))<Pc); %���ݽ�����ʽ��н���
for i=1:length(poscut),%length��A����ʾA�������/����
     newench2cross([2*i-1 2*i],:)=[newench2select([2*i-1 2*i],1:poscut(i)) newench2select([2*i 2*i-1],poscut(i)+1:lchrom)];%���㽻��
end
