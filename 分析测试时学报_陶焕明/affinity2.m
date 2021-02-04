function [Ab]=affinity2(ench2,popsize,len)
for i=1:popsize
    for j=i:popsize
        gene1=ench2(i,:);
        gene2=ench2(j,:);
        s=0;
        for k=1:len
            if gene1(k)==gene2(k)
                Hj(k)=0;
            else
                Hj(k)=log10(2);
            end
            s=s+Hj(k);
        end
        H_2=s/len;
        Ab(i,j)=1/(1+H_2);
        Ab(j,i)=Ab(i,j);
    end
end