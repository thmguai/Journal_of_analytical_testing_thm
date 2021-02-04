function [newchselect]=Select(ench2,E,popsize)
sumE=sum(E);
pE=E/sumE;
psE=0;
psE(1)=pE(1);
for i=2:popsize
    psE(i)=psE(i-1)+pE(i);
end
for i=1:popsize
    sita=rand;
    for g=1:popsize
        if sita<=psE(g)
            n=g;
            break;
        end
    end
    newchselect(i,:)=ench2(n,:);
end

