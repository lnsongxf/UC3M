function error = ar1_error2(para,data)

[parcol,~]=size(para);
trans=(parcol-2)/2;
perms=(parcol-2)/2;

sigma_i = para(1);
sigma_e = para(2);
theta = para(3);
rho = para(4);

col = length(data);
model = zeros(col,1);

% When kappa=0 
kappa=0;
model(kappa+1,1)=sigma_e+(rho-1)^2*sigma_e/(1-rho^2)+sigma_i*(1+(theta-1)^2+theta^2);
index0=1+kappa+trans;
model(index0,1)=sigma_e+(rho-1)^2*sigma_e/(1-rho^2)+sigma_i*(1+(theta-1)^2+theta^2);
for jj=3:trans-kappa   
    index0=trans-(jj-2)+index0;
    model(index0,1)=sigma_e+(rho-1)^2*sigma_e/(1-rho^2)+sigma_i*(1+(theta-1)^2+theta^2);
end
% When kappa=1
kappa=1;
model(kappa+1,1)=(rho-1)*sigma_e+rho*(rho-1)^2*sigma_e/(1-rho^2)+sigma_i*((theta-1)*(1-theta));
index1=1+kappa+trans;
model(index1,1)=(rho-1)*sigma_e+rho*(rho-1)^2*sigma_e/(1-rho^2)+sigma_i*((theta-1)*(1-theta));
for jj=3:trans-kappa
    index1=trans-(jj-2)+index1;
    model(index1,1)=(rho-1)*sigma_e+rho*(rho-1)^2*sigma_e/(1-rho^2)+sigma_i*((theta-1)*(1-theta));  
end
% When kappa=2
kappa=2;
model(kappa+1,1)=rho*(rho-1)*sigma_e+rho^2*(rho-1)^2*sigma_e/(1-rho^2)-theta*sigma_i;
index2=1+kappa+trans;
model(index2,1)=rho*(rho-1)*sigma_e+rho^2*(rho-1)^2*sigma_e/(1-rho^2)-theta*sigma_i;
for jj=3:trans-kappa
    index2=trans-(jj-2)+index2;
    model(index2,1)=rho*(rho-1)*sigma_e+rho^2*(rho-1)^2*sigma_e/(1-rho^2)-theta*sigma_i;  
end
% When kappa=3
kappa=3;
model(kappa+1,1)=rho^2*(rho-1)*sigma_e+rho^3*(rho-1)^2*sigma_e/(1-rho^2);
index3=1+kappa+trans;
model(index3,1)=rho^2*(rho-1)*sigma_e+rho^3*(rho-1)^2*sigma_e/(1-rho^2);
for jj=3:trans-kappa
    index3=trans-(jj-2)+index3;
    model(index3,1)=rho^2*(rho-1)*sigma_e+rho^3*(rho-1)^2*sigma_e/(1-rho^2);  
end
% When kappa>3
for kappa=4:30
    model(kappa+1,1)=rho^(kappa-1)*(rho-1)*sigma_e+rho^kappa*(rho-1)^2*sigma_e/(1-rho^2);
    indexn=1+kappa+trans;
    model(indexn,1)=rho^(kappa-1)*(rho-1)*sigma_e+rho^kappa*(rho-1)^2*sigma_e/(1-rho^2);
    for jj=3:trans-kappa
    indexn=trans-(jj-2)+indexn;
    model(indexn,1)=rho^(kappa-1)*(rho-1)*sigma_e+rho^kappa*(rho-1)^2*sigma_e/(1-rho^2);  
    end
end 

I=eye(col);
error=(-data+model)'*I*(-data+model);

end