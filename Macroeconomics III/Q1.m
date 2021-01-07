clear all
close all
clc

%% 

load('mvh.mat','mvh'); %covariance matrix 

%number of years
[col,~]=size(mvh); 

%write lower diagonel matrix into vector
mvh=vech(mvh);

%number of estimated parameters
trans = col; %number of transitory variances
perms = col; %number of persistent variances
NPARMS=trans+perms+2;

initial=zeros(NPARMS,1);
initial(1:trans)=0.047*ones(trans,1); %variance transitory
initial(trans+1:trans+perms)=0.027*ones(trans,1); %variance persistent
initial(trans+perms+1)=0.0439; %MA(1)
initial(trans+perms+2)=0.99; %rho
lb=zeros(NPARMS,1);
lb(1:trans)=0;
lb(trans+1:trans+perms)=0;
lb(trans+perms+1)=-1;
lb(trans+perms+2)=0;
ub=zeros(NPARMS,1);
ub(1:trans)=Inf;
ub(trans+1:trans+perms)=Inf;
ub(trans+perms+1)=1;
ub(trans+perms+2)=1;
options =optimoptions('fmincon','Algorithm','sqp','TolCon',1e-10,'TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',10000);
[paraopt,fval]=fmincon(@(x) ar1_error(x,mvh),initial,[],[],[],[],lb,ub,[],options);
age_varying = paraopt';
%% Bootstrap Standard Errors

load('mvh2.mat','mvh2'); %covariance matrix

boots = zeros(NPARMS,100);
for j = 1:100
    
    %write lower diagonel matrix into vector
    mvh=vech(mvh2(:,:,j));
    
    initial=zeros(NPARMS,1);
    initial(1:trans)=0.047*ones(trans,1); %sigma transitory
    initial(trans+1:trans+perms)=0.027*ones(trans,1); %sigma persistent
    initial(trans+perms+1)=0.0439; %MA(1)
    initial(trans+perms+2)=0.99; %rho
    lb=zeros(NPARMS,1);
    lb(1:trans)=0;
    lb(trans+1:trans+perms)=0;
    lb(trans+perms+1)=-1;
    lb(trans+perms+2)=0;
    ub=zeros(NPARMS,1);
    ub(1:trans)=Inf;
    ub(trans+1:trans+perms)=Inf;
    ub(trans+perms+1)=1;
    ub(trans+perms+2)=1;
    options =optimoptions('fmincon','Algorithm','sqp','TolCon',1e-10,'TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',10000);
    
    [paraopt,fval]=fmincon(@(x) ar1_error(x,mvh),initial,[],[],[],[],lb,ub,[],options);
    boots(:,j) = paraopt';
end

%obtain 5th and 95th point estimates
perc5 = prctile(boots,5,2);
perc95 = prctile(boots,95,2);

%% Age invariant moments

load('mvh.mat','mvh'); %covariance matrix 

[col,~]=size(mvh); 

%write lower diagonel matrix into vector
mvh=vech(mvh);

NPARMS=4;

initial=zeros(NPARMS,1);
initial(1)=0.047; %sigma transitory
initial(2)=0.027; %sigma persistent
initial(3)=0.0439; %MA(1)
initial(4)=0.99; %rho
lb=zeros(NPARMS,1);
lb(1)=0;
lb(2)=0;
lb(3)=-1;
lb(4)=0;
ub=zeros(NPARMS,1);
ub(1)=Inf;
ub(2)=Inf;
ub(3)=1;
ub(4)=1;
options =optimoptions('fmincon','Algorithm','sqp','TolCon',1e-10,'TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',10000);

[paraopt,fval]=fmincon(@(x) ar1_error2(x,mvh),initial,[],[],[],[],lb,ub,[],options);
age_invariant = paraopt';

%% Graphs

x=[25:55];

figure(1)
plot(x,age_varying(1:trans),'r-o','LineWidth',2.5)
hold on
plot(x,age_invariant(1,1)*ones(trans,1),'b-d','LineWidth',1.5)
plot(x,perc5(1:trans),'r--','LineWidth',1.0)
plot(x,perc95(1:trans),'r--','LineWidth',1.0)
xlim([25 55])
grid on
xlabel('Age','FontSize',14,'FontWeight','bold')
ylabel('Variance','FontSize',14,'FontWeight','bold')
legend('Age-variant','Age-invariant','FontSize',14,'FontWeight','bold')
hold off
saveas(gcf,'Transitory','epsc')

figure(2)
plot(x,age_varying(trans+1:trans+perms),'r-o','LineWidth',2.5)
hold on
plot(x,age_invariant(1,2)*ones(trans,1),'b-d','LineWidth',1.5)
plot(x,perc5(trans+1:trans+perms),'r--','LineWidth',1.0)
plot(x,perc95(trans+1:trans+perms),'r--','LineWidth',1.0)
xlim([25 55])
grid on
xlabel('Age','FontSize',14,'FontWeight','bold')
ylabel('Variance','FontSize',14,'FontWeight','bold')
legend('Age-variant','Age-invariant','FontSize',14,'FontWeight','bold')
hold off
saveas(gcf,'Persistent','epsc')
