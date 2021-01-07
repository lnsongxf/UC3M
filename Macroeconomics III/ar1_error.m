function error = ar1_error(para,data)

[parcol,~]=size(para);
trans=(parcol-2)/2;
perms=(parcol-2)/2;

sigma_i = para(1:trans);
sigma_e = para(trans+1:trans+perms);
theta = para(trans+perms+1);
rho = para(trans+perms+2);

col = length(data);
model = zeros(col,1);

% When kappa=0 
kappa=0;
model(kappa+1,1)=sigma_e(1)+sigma_i(1);
index0=1+kappa+trans;
model(index0,1)=sigma_e(2)+(rho-1)^2*sigma_e(1)+sigma_i(2)+(theta-1)*sigma_i(1);
for jj=3:trans-kappa
    ee=1;
    SUM0=0;
    while ee<jj
        SUM0=rho^(2*(ee-1))*sigma_e(jj-ee)+SUM0;
        ee=ee+1;
    end    
    index0=trans-(jj-2)+index0;
    model(index0,1)=sigma_e(jj)+(rho-1)^2*SUM0+sigma_i(jj)+(theta-1)*sigma_i(jj-1)+theta^2*sigma_i(jj-2);
end
% When kappa=1
kappa=1;
model(kappa+1,1)=(rho-1)*sigma_e(1)+(theta-1)*sigma_i(1);
index1=1+kappa+trans;
model(index1,1)=(rho-1)^2*rho*sigma_e(1)+(rho-1)*sigma_e(2)+sigma_i(2)-(theta-1)*theta*sigma_i(1);
for jj=3:trans-kappa
    ee=1;
    SUM1=0;
    while ee<jj
        SUM1=rho^(2*(ee-1))*sigma_e(jj-ee)+SUM1;
        ee=ee+1;
    end
    index1=trans-(jj-2)+index1;
    model(index1,1)=(rho-1)*sigma_e(jj)+(rho-1)^2*rho*SUM1+sigma_i(jj)+(theta-1)*sigma_i(jj)-theta*(theta-1)*sigma_i(jj-1);  
end
% When kappa=2
kappa=2;
model(kappa+1,1)=rho*(rho-1)*sigma_e(1)-theta*sigma_i(1);
index2=1+kappa+trans;
model(index2,1)=(rho-1)^2*rho^2*sigma_e(1)+rho*(rho-1)*sigma_e(2)-theta*sigma_i(2);
for jj=3:trans-kappa
    ee=1;
    SUM2=0;
    while ee<jj
        SUM2=rho^(2*(ee-1))*sigma_e(jj-ee)+SUM2;
        ee=ee+1;
    end
    index2=trans-(jj-2)+index2;
    model(index2,1)=rho*(rho-1)*sigma_e(jj)+(rho-1)^2*rho^2*SUM2-theta*sigma_i(jj);  
end
% When kappa=3
kappa=3;
model(kappa+1,1)=rho^2*(rho-1)*sigma_e(1);
index3=1+kappa+trans;
model(index3,1)=(rho-1)^2*rho^3*sigma_e(1)+rho^2*(rho-1)*sigma_e(2);
for jj=3:trans-kappa
    ee=1;
    SUM3=0;
    while ee<jj
        SUM3=rho^(2*(ee-1))*sigma_e(jj-ee)+SUM3;
        ee=ee+1;
    end
    index3=trans-(jj-2)+index3;
    model(index3,1)=rho^2*(rho-1)*sigma_e(jj)+(rho-1)^2*rho^3*SUM3;  
end
% When kappa>3
for kappa=4:30
    model(kappa+1,1)=rho^(kappa-1)*(rho-1)*sigma_e(1);
    indexn=1+kappa+trans;
    model(indexn,1)=(rho-1)^2*rho^kappa*sigma_e(1)+rho^(kappa-1)*(rho-1)*sigma_e(2);
    for jj=3:trans-kappa
    ee=1;
    SUMN=0;
    while ee<jj
        SUMN=rho^(2*(ee-1))*sigma_e(jj-ee)+SUMN;
        ee=ee+1;
    end
    indexn=trans-(jj-2)+indexn;
    model(indexn,1)=rho^(kappa-1)*(rho-1)*sigma_e(jj)+(rho-1)^2*rho^kappa*SUMN;  
    end
end    
I=eye(col);
error=(-data+model)'*I*(-data+model);
end