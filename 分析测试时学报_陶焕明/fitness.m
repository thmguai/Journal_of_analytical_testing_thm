function z = fitness(x,X,L)
    [n,m]=size(X);
    j = 1;
    f = 7;
    for i=1:L
        if x(i) == 1 && i == L
            XX(:,j:j+m-(4*i-2)) = X(:,4*i-2:end);
        elseif x(i) == 1
            XX(:,j:j+3) = X(:,4*i-2:4*i+1);
            j = j+4;
        end
    end
    YY = X(:,1);
    [xl,yl,xs,ys,beta,pctvar,mse]=plsregress(XX,YY,f);%对xr和Y进行pls回归
    RMSEC = sqrt(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2)/52);

    R = sqrt(1-(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2))/(sum((YY-mean(YY)).^2)));

    z = R/(1+RMSEC);