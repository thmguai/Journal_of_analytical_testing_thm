function [E]=Select_Antibody(C,Ag,popsize,k,K)
lanm=0.7;
miu=1.25;
%------------·½·¨1------------
for i=1:popsize
    E(i)=lanm*Ag(i)+(1-lanm)*exp(-miu*C(i));
end
