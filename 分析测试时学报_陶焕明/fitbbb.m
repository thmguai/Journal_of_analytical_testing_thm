function [R,RMSEP,y] = fitbbb(x,X,L,beta)     %适应度函数

    [n,m]=size(X);
    %X=xlsread('C:\Users\c430\Desktop\GA\HY1.xlsx');%导入数据
    j = 1;
    %f = 7;
    for i=1:L
        if x(i) == 1 && i == L
            XX(:,j:j+m-(4*i-2)) = X(:,4*i-2:end);
        elseif x(i) == 1
            XX(:,j:j+3) = X(:,4*i-2:4*i+1);
            j = j+4;
        end
    end
    YY = X(:,1);
    RMSEP = sqrt(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2)/28);
    y = XX*beta(2:end,:)+beta(1,:);
    R = sqrt(1-(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2))/(sum((YY-mean(YY)).^2)));

    z = R/(1+RMSEP);
    %disp(z)
end