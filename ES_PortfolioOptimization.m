%==========================================================================
% Quantitative Asset & Risk Managment
% Asset Project
%Portfolio optimization under Expected Shortfall
%
% Maxime Borel, Joe Habib, Yoann Joye & Axel Bandelier
% Date: 21 April 2020
%==========================================================================
clear; 
clc; 
close all;

restoredefaultpath      % restore search path to default
addpath(genpath('/Users/maxime/Documents/Universite?/HEC/MscF/4.2/QARM I/QARM_Assignment_AM/Asset project v4/function'))
clear RESTOREDEFAULTPATH_EXECUTED
%% import data
data = importdata...
    ('/Users/maxime/Documents/Universite?/HEC/MscF/4.2/QARM I/QARM_Assignment_AM/Asset project v3/data/17_industry.xlsx');
GS102 = importdata....
    ('/Users/maxime/Documents/Universite?/HEC/MscF/4.2/QARM I/QARM_Assignment_AM/Asset project v3/data/GS10-2.xls');
marketportfolio=importdata...
    ('/Users/maxime/Documents/Universite?/HEC/MscF/4.2/QARM I/QARM_Assignment_AM/Asset project v3/data/Marketreturn.xlsx');

date=datetime(data.data(:,1),'ConvertFrom','excel');
Return=data.data(:,2:end);
[T,N]=size(Return);
names=data.textdata(1,2:end);
e = ones(N,1);
rollingwindow=360; %set up the rolling window for 360 months
optimperiod=T-rollingwindow; %this is is the number of day in the optimisation period

Marketreturn=marketportfolio.data(:,6)./100; %we import the market return from FF
Marketreturn=Marketreturn(optimperiod:end); %we consider the performance only on the optimisation period
Rf=GS102(:,2); %import the riskfree rate

%clean data, replace the -99.99 and -999 by 0 
for i=1:N
    idx=Return(:,i)==-99.9900;
    Return(idx,i)=0;
    idx=Return(:,i)==-999;
    Return(idx,i)=0;
end

Return=Return./100;
Rf=Rf./1200; % depercentage and deannualize whole set of data
RF=Rf(optimperiod:end); %consider rf only on the optimization period
%% Market Portfolio
marketportfolio.perf=Marketreturn;
marketportfolio.cumperf=cumprod(1+marketportfolio.perf); %compute the cumulative performance
figure
plot(date(optimperiod:end,1),marketportfolio.cumperf)
%% optimsetting
w0=ones(N,1)./N; %initial value for optim as equaly weighted
Const.lb = zeros(N +rollingwindow +1,1);  % Lowerbound, no short sale 
Const.ub =  [0.2*ones(N,1) ;inf ;inf*ones(rollingwindow,1)]; % Upperbound
Const.sumW=1;
Const.conflvl95=0.95; % confidence level
Const.CVaR=[0.095;zeros(rollingwindow,1)]; %inequality restriction 
Const.CVaR2=zeros(rollingwindow,1); %need this for model 3 as no restriction on CVaR
x0=[w0;zeros(rollingwindow+1,1)]; %total inital value for optimization
Aeq=[ones(1,N) zeros(1,rollingwindow+1)]; %linear constraint that we want the sum of weight equal 1
ratio=1/(1-Const.conflvl95)*1/rollingwindow;

options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'MaxIterations',10000,'ConstraintTolerance',1.0000e-6, ...
    'OptimalityTolerance',1.0000e-6,'MaxFunctionEvaluations',...
    100000,'Display','none');

% setting all the variable we want to compute
modeltest.W=zeros(N+rollingwindow+1,optimperiod); %in this vector we have all weight the VaR and all the auxilary variables
modeltest.perf=zeros(optimperiod,1);
modeltest.fval=zeros(optimperiod,1);
modeltest.CVaR=zeros(optimperiod,1);

modeltest3.W=zeros(N+rollingwindow+1,optimperiod);
modeltest3.perf=zeros(optimperiod,1);
modeltest3.CVaR=zeros(optimperiod,1);

modeltest5.W=zeros(N+rollingwindow+1,optimperiod);
modeltest5.perf=zeros(optimperiod,1);
modeltest5.CVaR=zeros(optimperiod,1);

MV.W=zeros(N,optimperiod);
MV.perf=zeros(optimperiod,1);

GMVP.W=zeros(N,optimperiod);
GMVP.perf=zeros(optimperiod,1);
%% optimization
%paper model
for t=1:optimperiod
    ReturnInRollingwindow=Return(t:t+rollingwindow-1,:);
    A=[zeros(1,N),1, ratio.*ones(1,rollingwindow);-ReturnInRollingwindow ,(-1).*ones(rollingwindow,1) ...
        ,(-1).*eye(rollingwindow)];% define linear restriction for CVaR and weight, we want u>=0 and the cvar<=treshold
    
    objfun1=@(W) -(mean(ReturnInRollingwindow)*W(1:N,1)); %minimize the loss function
    [modeltest.W(:,t),modeltest.fval(t,1)]=fmincon(objfun1,x0,...
        A,Const.CVaR,Aeq,Const.sumW,Const.lb,Const.ub,[],options);
    
    modeltest.perf(t,1)=Return(t+rollingwindow,:)*modeltest.W(1:N,t); %compute the performance for next period
    modeltest.CVaR(t,1)=modeltest.W(N+1,t)+ratio*sum(modeltest.W(N+2:end,t)); %CVaR= VaR + sum of auxilary variable
end
modeltest.fval(:,1)=modeltest.fval(:,1)*-1; %expected return for each period

%sharpe ratio
for t=1:optimperiod
    ReturnInRollingwindow=Return(t:t+rollingwindow-1,:);
    expectedreturn=mean(ReturnInRollingwindow);
    sigma=cov(ReturnInRollingwindow);
    A=[zeros(1,N),1, 1/(rollingwindow*(1-Const.conflvl95)).*ones(1,rollingwindow);...
        -ReturnInRollingwindow ,(-1).*ones(rollingwindow,1) ,(-1).*eye(rollingwindow)]; %same as before 
    
    objfun3=@(W) -((expectedreturn*W(1:N,1)-Rf(t+rollingwindow-1))/portfoliostd(W(1:N,1),sigma)); %max the sharpe ratio we use r for the actual period as we know it now
    modeltest3.W(:,t)=fmincon(objfun3,x0,...
        A,Const.CVaR,Aeq,Const.sumW,Const.lb,Const.ub,[],options);
    
    modeltest3.perf(t,1)=Return(t+rollingwindow,:)*modeltest3.W(1:N,t); 
    modeltest3.CVaR(t,1)=modeltest3.W(N+1,t)+ratio*sum(modeltest3.W(N+2:end,t)); 
end
%ES in the mean variance
for t=1:optimperiod
    ReturnInRollingwindow=Return(t:t+rollingwindow-1,:);
    expectedreturn=mean(ReturnInRollingwindow);
    A2=[-ReturnInRollingwindow ,(-1).*ones(rollingwindow,1) ,(-1).*eye(rollingwindow)]; %we change the formualtion as no more 
                                                                                        %restriction on the CVaR just previous
                                                                                        %juste u>=0 and sum of weight
    objfun5=@(W) -(expectedreturn*W(1:N,1)-(W(N+1,1)+ratio*sum(W(N+2:end,1))));
    modeltest5.W(:,t)=fmincon(objfun5,x0,...
        A2,Const.CVaR2,Aeq,Const.sumW,Const.lb,Const.ub,[],options);
    
    modeltest5.perf(t,1)=Return(t+rollingwindow,:)*modeltest5.W(1:N,t);
    modeltest5.CVaR(t,1)=modeltest5.W(N+1,t)+ratio*sum(modeltest5.W(N+2:end,t));
end

%MV
for t=1:optimperiod
    ReturnInRollingwindow=Return(t:t+rollingwindow-1,:);
    expectedreturn=mean(ReturnInRollingwindow);
    sigma=cov(ReturnInRollingwindow);
    
    objfun6=@(W) -(expectedreturn*W-0.5.*W'*sigma*W); %classic mean variance optimization lambda equal 1
    MV.W(:,t)=fmincon(objfun6,w0,...
        [],[],e',Const.sumW,Const.lb(1:N,1),Const.ub(1:N,1),[],options);
    
    MV.perf(t,1)=Return(t+rollingwindow,:)*MV.W(:,t);
end

for t=1:optimperiod
    ReturnInRollingwindow=Return(t:t+rollingwindow-1,:);
    sigma=cov(ReturnInRollingwindow);
    
    objfun7=@(W) W'*sigma*W;%minimum variance portfolio
    GMVP.W(:,t)=fmincon(objfun7,w0,...
        [],[],e',Const.sumW,Const.lb(1:N,1),Const.ub(1:N,1),[],options);
    
    GMVP.perf(t,1)=Return(t+rollingwindow,:)*GMVP.W(:,t);
end
MV.TO=turnoverp(MV.W,MV.perf,Return(rollingwindow:end,:)); %compute the turnover of all portfolio
GMVP.TO=turnoverp(GMVP.W,GMVP.perf,Return(rollingwindow:end,:));
Model1TO=turnoverp(modeltest.W(1:17,:),modeltest.perf,Return(rollingwindow:end,:));
Model2TO=turnoverp(modeltest3.W(1:17,:),modeltest3.perf,Return(rollingwindow:end,:));
Model3TO=turnoverp(modeltest5.W(1:17,:),modeltest5.perf,Return(rollingwindow:end,:));
TOmodel=table(Model1TO,Model2TO,Model3TO,'Rownames',{'Turnover'},'VariableNames',{'Model1','Model2','Model3'});
%% Performance
%cumulative performance
modeltest.cumperf=cumprod(1+modeltest.perf);
modeltest3.cumperf=cumprod(1+modeltest3.perf);
modeltest5.cumperf=cumprod(1+modeltest5.perf);
GMVP.cumperf=cumprod(1+GMVP.perf);
MV.cumperf=cumprod(1+MV.perf);

%Comparative table
Names_ESportfolio={'Model 1','Model 2','Model 3'};
ReturnAllPortfolio=[modeltest.perf,modeltest3.perf,modeltest5.perf];
SumStat_ESp=Portfoliostatistics(ReturnAllPortfolio,RF,Const.conflvl95,Names_ESportfolio);
SumStat_ESp=table([SumStat_ESp;TOmodel]);

Names_M1={'Model 1','Market','GMVP','MV'};
ReturnModel1=[modeltest.perf,marketportfolio.perf,GMVP.perf,MV.perf];
SumStat_M1=Portfoliostatistics(ReturnModel1,RF,Const.conflvl95,Names_M1);
TOM1=[Model1TO,'-',GMVP.TO,MV.TO];

Names_M2={'Model 2','Market','GMVP','MV'};
ReturnModel2=[modeltest3.perf,marketportfolio.perf,GMVP.perf,MV.perf];
SumStat_M2=Portfoliostatistics(ReturnModel2,RF,Const.conflvl95,Names_M2);
TOM2=[Model2TO,'-',GMVP.TO,MV.TO];

Names_M3={'Model 3','Market','GMVP','MV'};
ReturnModel3=[modeltest5.perf,marketportfolio.perf,GMVP.perf,MV.perf];
SumStat_M3=Portfoliostatistics(ReturnModel3,RF,Const.conflvl95,Names_M3);
TOM3=[Model3TO,'-',GMVP.TO,MV.TO];
%% Weights
modeltest.avgweight=mean(modeltest.W(1:N),1); % take the mean over time for each asset class of the weight
modeltest3.avgweight=mean(modeltest3.W(1:N),1);
modeltest5.avgweight=mean(modeltest5.W(1:N),1);

p1=figure; %show it as the a pie graph
pie(modeltest.avgweight)
legend(names,'Location','northwest')
title('Average weights Model 1 at 95%')
saveas(p1,'results/Average weights Model 1 at 95','png');

p2=figure;
pie(modeltest3.avgweight)
legend(names,'Location','northwest')
title('Average weights Model 2 at 95%')
saveas(p2,'results/Average weights Model 2 at 95','png');

p3=figure;
pie(modeltest5.avgweight)
legend(names,'Location','northwest')
title('Average weights Model 3 at 95%')
saveas(p3,'results/Average weights Model 3 at 95','png');

%% graph
%we represent the evolution of the weight over time for each portfolio
%max return st ES<9.5%
f1=figure;
area(date((1+rollingwindow):end),modeltest.W(1:N,:)','FaceColor','flat');
ylim([0,1])
legend(names)
title('max return st ES<9.5%')
saveas(f1,'results/max return st ES','png');

%max SR st ES<9.5%
f2=figure;
area(date((1+rollingwindow):end),modeltest3.W(1:N,:)','FaceColor','flat');
ylim([0,1])
legend(names)
title('max SR st ES<9.5%')
saveas(f2,'results/max SR st ES','png');

%max expected return -ES
f3=figure;
area(date((1+rollingwindow):end),modeltest5.W(1:N,:)','FaceColor','flat');
ylim([0,1])
legend(names)
title('max expected return -ES')
saveas(f3,'results/max expected return -ES','png');

%MV
f4=figure;
area(date((1+rollingwindow):end),MV.W','FaceColor','flat');
ylim([0,1])
legend(names)
title('MV')
saveas(f4,'results/MV','png');

%GMVP
f5=figure;
area(date((1+rollingwindow):end),GMVP.W','FaceColor','flat');
ylim([0,1])
legend(names)
title('GMVP')
saveas(f5,'results/GMVP','png');

%Expected Shortfall illustration
f6=figure;
Varhist=quantile(modeltest.perf,0.05);
mask=modeltest.perf < Varhist;
histogram(modeltest.perf(mask),25,'FaceColor','r')
hold on
histogram(modeltest.perf(~mask),35,'FaceColor','b')
xline(Varhist,'--r');
xline(mean(modeltest.perf(modeltest.perf<Varhist)),'--r');
title('Expected Shortfall illustration')
saveas(f6,'results/Expected Shortfall illustration','png');

%plot all performance
f7=figure;
plot(date(optimperiod:end),modeltest.cumperf,...
    date(optimperiod:end),modeltest3.cumperf,...
    date(optimperiod:end),modeltest5.cumperf,...
    date(optimperiod:end),GMVP.cumperf,...
    date(optimperiod:end),marketportfolio.cumperf,...
    date(optimperiod:end),MV.cumperf,...
    'LineWidth',1)
legend('Model 1','Model 2','Model 3','GMVP',...
    'Market','MV','Location','northwest')
title('Cumulative performance of portfolio at 95%')
saveas(f7,'results/Cumulative performance of portfolio at 95%','png');
%% Last decade 2010-2020
%same as upper but we consider only part of the time the below is almost
%copy past
date2=date(optimperiod:end);
LD=length(date2)-10*12; %we have 30y of perf and we want the last 10 so 10y *12 give starting point 

modeltest.perfLD=modeltest.perf(LD:end);
modeltest3.perfLD=modeltest3.perf(LD:end);
modeltest5.perfLD=modeltest5.perf(LD:end);
GMVP.perfLD=GMVP.perf(LD:end);
MV.perfLD=MV.perf(LD:end);
marketportfolio.perfLD=marketportfolio.perf(LD:end);

%turnonver last decade
MV.TOLD=turnoverp(MV.W(:,LD:end),MV.perfLD,Return(LD:end,:));
GMVP.TOLD=turnoverp(GMVP.W(:,LD:end),GMVP.perfLD,Return(LD:end,:));
Model1TOLD=turnoverp(modeltest.W(1:17,LD:end),modeltest.perfLD,Return(LD:end,:));
Model2TOLD=turnoverp(modeltest3.W(1:17,LD:end),modeltest3.perfLD,Return(LD:end,:));
Model3TOLD=turnoverp(modeltest5.W(1:17,LD:end),modeltest5.perfLD,Return(LD:end,:));

%cumulative performance last decade
modeltest.cumperfLD=cumprod(1+modeltest.perfLD);
modeltest3.cumperfLD=cumprod(1+modeltest3.perfLD);
modeltest5.cumperfLD=cumprod(1+modeltest5.perfLD);
GMVP.cumperfLD=cumprod(1+GMVP.perfLD);
MV.cumperfLD=cumprod(1+MV.perfLD);
marketportfolio.cumperfLD=cumprod(1+marketportfolio.perfLD);

%Comparative table last decade
ReturnAllPortfolio=[modeltest.perfLD,modeltest3.perfLD,modeltest5.perfLD];
SumStat_ESp_LD=Portfoliostatistics(ReturnAllPortfolio,...
    RF(LD:end),Const.conflvl95,Names_ESportfolio);

ReturnModel1=[modeltest.perfLD,marketportfolio.perf(LD:end),...
    GMVP.perfLD,MV.perfLD];
SumStat_M1_LD=Portfoliostatistics(ReturnModel1,RF(LD:end),...
    Const.conflvl95,Names_M1);

ReturnModel2=[modeltest3.perfLD,marketportfolio.perf(LD:end),...
    GMVP.perfLD,MV.perfLD];
SumStat_M2_LD=Portfoliostatistics(ReturnModel2,RF(LD:end),...
    Const.conflvl95,Names_M2);

ReturnModel3=[modeltest5.perfLD,marketportfolio.perf(LD:end),...
    GMVP.perfLD,MV.perfLD];
SumStat_M3_LD=Portfoliostatistics(ReturnModel3,RF(LD:end),...
    Const.conflvl95,Names_M3);

f8=figure;
plot(date2(LD:end),modeltest.cumperfLD,...
    date2(LD:end),modeltest3.cumperfLD,...
    date2(LD:end),modeltest5.cumperfLD,...
    date2(LD:end),GMVP.cumperfLD,...
    date2(LD:end),marketportfolio.cumperfLD,...
    date2(LD:end),MV.cumperfLD,...
    'LineWidth',1)
legend('Model 1','Model 2','Model 3','GMVP',...
    'Market','MV','Location','northwest')
title('Cumulative performance of portfolio at 95% 2010-2020')
saveas(f8,'results/Cumulative performance of portfolio at 95% 2010-2020','png');

%% middle decade 2000-2010
MD=length(date2)-LD+1; %we have 30y of perf and we want the last 10 so 10y *12 give starting point 

modeltest.perfMD=modeltest.perf(MD:LD-1);
modeltest3.perfMD=modeltest3.perf(MD:LD-1);
modeltest5.perfMD=modeltest5.perf(MD:LD-1);
GMVP.perfMD=GMVP.perf(MD:LD-1);
MV.perfMD=MV.perf(MD:LD-1);
marketportfolio.perfMD=marketportfolio.perf(MD:LD-1);

%turnonver mid decade
MV.TOMD=turnoverp(MV.W(MD:LD-1),MV.perfMD,Return(MD:LD-1,:));
GMVP.TOMD=turnoverp(GMVP.W(MD:LD-1),GMVP.perfMD,Return(MD:LD-1,:));
Model1TOMD=turnoverp(modeltest.W(1:17,MD:LD-1),modeltest.perfMD,Return(MD:LD-1,:));
Model2TOMD=turnoverp(modeltest3.W(1:17,MD:LD-1),modeltest3.perfMD,Return(MD:LD-1,:));
Model3TOMD=turnoverp(modeltest5.W(1:17,MD:LD-1),modeltest5.perfMD,Return(MD:LD-1,:));

%cumulative performance last decade
modeltest.cumperfMD=cumprod(1+modeltest.perfMD);
modeltest3.cumperfMD=cumprod(1+modeltest3.perfMD);
modeltest5.cumperfMD=cumprod(1+modeltest5.perfMD);
GMVP.cumperfMD=cumprod(1+GMVP.perfMD);
MV.cumperfMD=cumprod(1+MV.perfMD);
marketportfolio.cumperfMD=cumprod(1+marketportfolio.perfMD);

%Comparative table last decade
ReturnAllPortfolio=[modeltest.perfMD,modeltest3.perfMD,...
    modeltest5.perfMD];
SumStat_ESp_MD=Portfoliostatistics(ReturnAllPortfolio,RF(MD:LD-1),...
    Const.conflvl95,Names_ESportfolio);

ReturnModel1=[modeltest.perfMD,marketportfolio.perfMD,...
    GMVP.perfMD,MV.perfMD];
SumStat_M1_MD=Portfoliostatistics(ReturnModel1,RF(MD:LD-1),...
    Const.conflvl95,Names_M1);

ReturnModel2=[modeltest3.perfMD,marketportfolio.perfMD,...
    GMVP.perfMD,MV.perfMD];
SumStat_M2_MD=Portfoliostatistics(ReturnModel2,RF(MD:LD-1),...
    Const.conflvl95,Names_M2);

ReturnModel3=[modeltest5.perfMD,marketportfolio.perfMD,...
    GMVP.perfMD,MV.perfMD];
SumStat_M3_MD=Portfoliostatistics(ReturnModel3,RF(MD:LD-1),...
    Const.conflvl95,Names_M3);

f9=figure;
plot(date2(MD:LD-1),modeltest.cumperfMD,...
    date2(MD:LD-1),modeltest3.cumperfMD,...
    date2(MD:LD-1),modeltest5.cumperfMD,...
    date2(MD:LD-1),GMVP.cumperfMD,...
    date2(MD:LD-1),marketportfolio.cumperfMD,...
    date2(MD:LD-1),MV.cumperfMD,...
    'LineWidth',1)
legend('Model 1','Model 2','Model 3','GMVP',...
    'Market','MV','Location','northwest')
title('Cumulative performance of portfolio at 95% 2000-2010')
saveas(f9,'results/Cumulative performance of portfolio at 95% 2000-2010','png');

%% First decade 1990-2000
modeltest.perfFD=modeltest.perf(1:MD-1);
modeltest3.perfFD=modeltest3.perf(1:MD-1);
modeltest5.perfFD=modeltest5.perf(1:MD-1);
GMVP.perfFD=GMVP.perf(1:MD-1);
MV.perfFD=MV.perf(1:MD-1);
marketportfolio.perfFD=marketportfolio.perf(1:MD-1);

%turnonver
MV.TOFD=turnoverp(MV.W(1:MD-1),MV.perfFD,Return(1:MD-1,:));
GMVP.TOFD=turnoverp(GMVP.W(1:MD-1),GMVP.perfFD,Return(1:MD-1,:));
Model1TOFD=turnoverp(modeltest.W(1:17,1:MD-1),modeltest.perfFD,Return(1:MD-1,:));
Model2TOFD=turnoverp(modeltest3.W(1:17,1:MD-1),modeltest3.perfFD,Return(1:MD-1,:));
Model3TOFD=turnoverp(modeltest5.W(1:17,1:MD-1),modeltest5.perfFD,Return(1:MD-1,:));

%cumulative performance last decade
modeltest.cumperfFD=cumprod(1+modeltest.perfFD);
modeltest3.cumperfFD=cumprod(1+modeltest3.perfFD);
modeltest5.cumperfFD=cumprod(1+modeltest5.perfFD);
GMVP.cumperfFD=cumprod(1+GMVP.perfFD);
MV.cumperfFD=cumprod(1+MV.perfFD);
marketportfolio.cumperfFD=cumprod(1+marketportfolio.perfFD);

%Comparative table last decade
ReturnAllPortfolio=[modeltest.perfFD,modeltest3.perfFD,...
    modeltest5.perfFD];
SumStat_ESp_FD=Portfoliostatistics(ReturnAllPortfolio,...
    RF(1:MD-1),Const.conflvl95,Names_ESportfolio);

ReturnModel1=[modeltest.perfFD,marketportfolio.perfFD,...
    GMVP.perfFD,MV.perfFD];
SumStat_M1_FD=Portfoliostatistics(ReturnModel1,RF(1:MD-1),...
    Const.conflvl95,Names_M1);

ReturnModel2=[modeltest3.perfFD,marketportfolio.perfFD,...
    GMVP.perfFD,MV.perfFD];
SumStat_M2_FD=Portfoliostatistics(ReturnModel2,RF(1:MD-1),...
    Const.conflvl95,Names_M2);

ReturnModel3=[modeltest5.perfFD,marketportfolio.perfFD,...
    GMVP.perfFD,MV.perfFD];
SumStat_M3_FD=Portfoliostatistics(ReturnModel3,RF(1:MD-1),...
    Const.conflvl95,Names_M3);

f10=figure;
plot(date2(1:MD-1),modeltest.cumperfFD,...
    date2(1:MD-1),modeltest3.cumperfFD,...
    date2(1:MD-1),modeltest5.cumperfFD,...
    date2(1:MD-1),GMVP.cumperfFD,...
    date2(1:MD-1),marketportfolio.cumperfFD,...
    date2(1:MD-1),MV.cumperfFD,...
    'LineWidth',1)
legend('Model 1','Model 2','Model 3','GMVP',...
    'Market','MV','Location','northwest')
title('Cumulative performance of portfolio at 95% 1990-2000')
saveas(f10,'results/Cumulative performance of portfolio at 95% 1990-2000','png');

%% sensitivity analysis  model1 the one of the paper
% Change conf lvl 95->99 not good idea because not comparable
% Change upper, lower bounds
% change the constrainte on ES
%rolling window 120, 240 
%change value of the lambda
%% impact on upperbond
modeltestub25.W=zeros(N+rollingwindow+1,optimperiod);
modeltestub25.perf=zeros(optimperiod,1);
modeltestub25.fval=zeros(optimperiod,1);
modeltestub25.CVaR=zeros(optimperiod,1);

modeltestub15.W=zeros(N+rollingwindow+1,optimperiod);
modeltestub15.perf=zeros(optimperiod,1);
modeltestub15.fval=zeros(optimperiod,1);
modeltestub15.CVaR=zeros(optimperiod,1);

Const.ub25=[0.25*ones(N,1) ;inf ;inf*ones(rollingwindow,1)];% change de value of the upperbond
Const.ub15=[0.15*ones(N,1) ;inf ;inf*ones(rollingwindow,1)];

for t=1:optimperiod %redo the optimization for model 1 
    ReturnInRollingwindow=Return(t:t+rollingwindow-1,:);
    A=[zeros(1,N),1, ratio.*ones(1,rollingwindow);...
        -ReturnInRollingwindow ,(-1).*ones(rollingwindow,1) ,(-1).*eye(rollingwindow)];
    
    objfun1=@(W) -(mean(ReturnInRollingwindow)*W(1:N,1));
    [modeltestub25.W(:,t),modeltestub25.fval(t,1)]=fmincon(objfun1,x0,...
        A,Const.CVaR,Aeq,Const.sumW,Const.lb,Const.ub25,[],options);
    
    modeltestub25.perf(t,1)=Return(t+rollingwindow,:)*modeltestub25.W(1:N,t);
    modeltestub25.CVaR(t,1)=modeltestub25.W(18,t)+ratio*sum(modeltestub25.W(19:end,t));
    
    [modeltestub15.W(:,t),modeltestub15.fval(t,1)]=fmincon(objfun1,x0,A,...
        Const.CVaR,Aeq,Const.sumW,Const.lb,Const.ub15,[],options);
    
    modeltestub15.perf(t,1)=Return(t+rollingwindow,:)*modeltestub15.W(1:N,t);
    modeltestub15.CVaR(t,1)=modeltestub15.W(18,t)+ratio*sum(modeltestub15.W(19:end,t));
end
modeltestub25.fval(:,1)=modeltestub25.fval(t,1)*-1;
modeltestub15.fval(:,1)=modeltestub15.fval(t,1)*-1;

f11=figure;
area(date((1+rollingwindow):end),modeltestub25.W(1:N,:)','FaceColor','flat');
ylim([0,1])
legend(names)
title('max return st ES<9.5% upper bond 25%')
saveas(f11,'results/max return st ES upper bond 25','png');

f12=figure;
area(date((1+rollingwindow):end),modeltestub15.W(1:N,:)','FaceColor','flat');
ylim([0,1])
legend(names)
title('max return st ES<9.5% upper bond 15%')
saveas(f12,'results/max return st ES upper bond 15','png');

Names_bound={'Model1 UB 15%','Model1','Model1 UB 25%'};
ReturnAllPortfolio=[modeltestub15.perf,modeltest.perf,modeltestub25.perf];
SumStat_UB=Portfoliostatistics(ReturnAllPortfolio,RF,Const.conflvl95,Names_bound); % resume in a table 
modeltestub25TO=turnoverp(modeltestub25.W(1:17,:),modeltestub25.perf,Return(rollingwindow:end,:)); %turnonver
modeltestub15TO=turnoverp(modeltestub15.W(1:17,:),modeltestub15.perf,Return(rollingwindow:end,:));
Model1TOB=[modeltestub15TO,Model1TO,modeltestub25TO];
%% change constraint on the minimum ES 
modeltestES12.W=zeros(N+rollingwindow+1,optimperiod);
modeltestES12.perf=zeros(optimperiod,1);
modeltestES12.fval=zeros(optimperiod,1);
modeltestES12.CVaR=zeros(optimperiod,1);

modeltestES145.W=zeros(N+rollingwindow+1,optimperiod);
modeltestES145.perf=zeros(optimperiod,1);
modeltestES145.fval=zeros(optimperiod,1);
modeltestES145.CVaR=zeros(optimperiod,1);

Const.CVaR12=[0.12;zeros(rollingwindow,1)]; %change the constraint vector for the ES
Const.CVaR145=[0.145;zeros(rollingwindow,1)];

for t=1:optimperiod %redo the optimization
    ReturnInRollingwindow=Return(t:t+rollingwindow-1,:);
    A=[zeros(1,N),1, ratio.*ones(1,rollingwindow);...
        -ReturnInRollingwindow ,(-1).*ones(rollingwindow,1) ,(-1).*eye(rollingwindow)];
    
    objfun1=@(W) -(mean(ReturnInRollingwindow)*W(1:N,1));
    [modeltestES12.W(:,t),modeltestES12.fval(t,1)]=fmincon(objfun1,x0,A,...
        Const.CVaR12,Aeq,Const.sumW,Const.lb,Const.ub,[],options);
    
    modeltestES12.perf(t,1)=Return(t+rollingwindow,:)*modeltestES12.W(1:N,t);
    modeltestES12.CVaR(t,1)=modeltestES12.W(18,t)+ratio*sum(modeltestES12.W(19:end,t));
    
    [modeltestES145.W(:,t),modeltestES145.fval(t,1)]=fmincon(objfun1,x0,...
        A,Const.CVaR145,Aeq,Const.sumW,Const.lb,Const.ub,[],options);
    
    modeltestES145.perf(t,1)=Return(t+rollingwindow,:)*modeltestES145.W(1:N,t);
    modeltestES145.CVaR(t,1)=modeltestES145.W(18,t)+ratio*sum(modeltestES145.W(19:end,t));
end
modeltestES12.fval(:,1)=modeltestES12.fval(t,1)*-1;
modeltestES145.fval(:,1)=modeltestES145.fval(t,1)*-1;

f13=figure;
area(date((1+rollingwindow):end),modeltestES12.W(1:N,:)','FaceColor','flat');
ylim([0,1])
legend(names)
title('max return st ES<12%')
saveas(f13,'results/max return st ES<12%','png');

f14=figure;
area(date((1+rollingwindow):end),modeltestES145.W(1:N,:)','FaceColor','flat');
ylim([0,1])
legend(names)
title('max return st ES<14.5% upper bond 15%')
saveas(f14,'results/max return st ES<145%','png');

Names_ES={'Model1 ','Model1 ES 12%','Model1 ES 145%'};
ReturnAllPortfolio=[modeltest.perf,modeltestES12.perf,modeltestES145.perf];
SumStat_ES=Portfoliostatistics(ReturnAllPortfolio,RF,Const.conflvl95,Names_ES);
modeltestES12TO=turnoverp(modeltestES12.W(1:17,:),modeltestES12.perf,Return(rollingwindow:end,:));
modeltestES145TO=turnoverp(modeltestES145.W(1:17,:),modeltestES145.perf,Return(rollingwindow:end,:));
Model1TOES=[Model1TO,modeltestES12TO,modeltestES145TO];
%% Rolling Window 120 240
%240
Return240=Return(121:end,:);% consider smaller data set we want the performance start 1990 to be compute smaller RW s
%we dont consider first 121 return.
rollingwindow=240; %new rolling window
Const.CVaR=[0.095;zeros(rollingwindow,1)]; %reset the the constraint for CVaR and all optim settign as before
x0=[w0;zeros(rollingwindow+1,1)];
Aeq=[ones(1,N) zeros(1,rollingwindow+1)];
ratio=1/(1-Const.conflvl95)*1/rollingwindow;
Const.lb = zeros(N +rollingwindow +1,1);  % Lowerbound, no short sale
Const.ub =  [0.2*ones(N,1) ;inf ;inf*ones(rollingwindow,1)]; % Upperbound

modeltestRW240.W=zeros(N+rollingwindow+1,optimperiod);
modeltestRW240.perf=zeros(optimperiod,1);
modeltestRW240.fval=zeros(optimperiod,1);
modeltestRW240.CVaR=zeros(optimperiod,1);

for t=1:optimperiod %redo optimization
    ReturnInRollingwindow=Return240(t:t+rollingwindow-1,:);
    A=[zeros(1,N),1, ratio.*ones(1,rollingwindow);...
        -ReturnInRollingwindow ,(-1).*ones(rollingwindow,1) ,(-1).*eye(rollingwindow)];
    
    objfun1=@(W) -(mean(ReturnInRollingwindow)*W(1:N,1));
    [modeltestRW240.W(:,t),modeltestRW240.fval(t,1)]=fmincon(objfun1,x0,A,...
        Const.CVaR,Aeq,Const.sumW,Const.lb,Const.ub,[],options);
    
    modeltestRW240.perf(t,1)=Return240(t+rollingwindow,:)*modeltestRW240.W(1:N,t);
    modeltestRW240.CVaR(t,1)=modeltestRW240.W(18,t)+ratio*sum(modeltestRW240.W(19:end,t));
end
modeltestRW240.fval(:,1)=modeltestRW240.fval(t,1)*-1;

rollingwindow=360; % set again the rolling window for graph purpose
f15=figure;
area(date((1+rollingwindow):end),modeltestRW240.W(1:N,:)','FaceColor','flat');
ylim([0,1])
legend(names)
title('max return st ES<9.5% rolling window 240 month')
saveas(f15,'results/max return st ES<9.5% RW 240 months%','png');

%120
Return120=Return(241:end,:); % do the same for 120 rolling window
rollingwindow=120;
Const.CVaR=[0.095;zeros(rollingwindow,1)];
x0=[w0;zeros(rollingwindow+1,1)];
Aeq=[ones(1,N) zeros(1,rollingwindow+1)];
ratio=1/(1-Const.conflvl95)*1/rollingwindow;
Const.lb = zeros(N +rollingwindow +1,1);  % Lowerbound
Const.ub =  [0.2*ones(N,1) ;inf ;inf*ones(rollingwindow,1)]; % Upperbound

modeltestRW120.W=zeros(N+rollingwindow+1,optimperiod);
modeltestRW120.perf=zeros(optimperiod,1);
modeltestRW120.fval=zeros(optimperiod,1);
modeltestRW120.CVaR=zeros(optimperiod,1);

for t=1:optimperiod
    ReturnInRollingwindow=Return120(t:t+rollingwindow-1,:);
    A=[zeros(1,N),1, ratio.*ones(1,rollingwindow);...
        -ReturnInRollingwindow ,(-1).*ones(rollingwindow,1) ,(-1).*eye(rollingwindow)];
    
    objfun1=@(W) -(mean(ReturnInRollingwindow)*W(1:N,1));
    [modeltestRW120.W(:,t),modeltestRW120.fval(t,1)]=fmincon(objfun1,x0,A,...
        Const.CVaR,Aeq,Const.sumW,Const.lb,Const.ub,[],options);
    
    modeltestRW120.perf(t,1)=Return120(t+rollingwindow,:)*modeltestRW120.W(1:N,t);
    modeltestRW120.CVaR(t,1)=modeltestRW120.W(18,t)+ratio*sum(modeltestRW120.W(19:end,t));
end
modeltestRW120.fval(:,1)=modeltestRW120.fval(t,1)*-1;

rollingwindow=360;
f16=figure;
area(date((1+rollingwindow):end),modeltestRW120.W(1:N,:)','FaceColor','flat');
ylim([0,1])
legend(names)
title('max return st ES<9.5% rolling window 120 month')
saveas(f16,'results/max return st ES<9.5% RW 120 months%','png');

Names_RW={'Model1 RW 120 ','Model1 RW 240','Model1'};
ReturnAllPortfolio=[modeltestRW120.perf,modeltestRW240.perf,modeltest.perf];
SumStat_RW=Portfoliostatistics(ReturnAllPortfolio,RF,Const.conflvl95,Names_RW);
modeltestRW240TO=turnoverp(modeltestRW240.W(1:17,:),modeltestRW240.perf,Return(rollingwindow:end,:));
modeltestRW120TO=turnoverp(modeltestRW120.W(1:17,:),modeltestRW120.perf,Return(rollingwindow:end,:));
Model1TORW=[Model1TO,modeltestRW120TO,modeltestRW240TO];
%% lambda=2 under model 3
%we reset all the perivious optimization setting as they have change. 
modeltest5l2.W=zeros(N+rollingwindow+1,optimperiod);
modeltest5l2.perf=zeros(optimperiod,1);
modeltest5l2.CVaR=zeros(optimperiod,1);
Const.lb = zeros(N +rollingwindow +1,1);  % Lowerbound
Const.ub =  [0.2*ones(N,1) ;inf ;inf*ones(rollingwindow,1)]; % Upperbound
Const.sumW=1;
Const.conflvl95=0.95;
Const.CVaR=[0.095;zeros(rollingwindow,1)];
Const.CVaR2=zeros(rollingwindow,1);
Const.lambda=2;
x0=[w0;zeros(rollingwindow+1,1)];
Aeq=[ones(1,N) zeros(1,rollingwindow+1)];
ratio=1/(1-Const.conflvl95)*1/rollingwindow;

for t=1:optimperiod
    ReturnInRollingwindow=Return(t:t+rollingwindow-1,:);
    expectedreturn=mean(ReturnInRollingwindow);
    A2=[-ReturnInRollingwindow ,(-1).*ones(rollingwindow,1) ,(-1).*eye(rollingwindow)];
    
    objfun5=@(W) -(expectedreturn*W(1:N,1)-Const.lambda.*(W(N+1,1)+ratio*sum(W(N+2:end,1))));
    modeltest5l2.W(:,t)=fmincon(objfun5,x0,...
        A2,Const.CVaR2,Aeq,Const.sumW,Const.lb,Const.ub,[],options);
    
    modeltest5l2.perf(t,1)=Return(t+rollingwindow,:)*modeltest5l2.W(1:N,t);
    modeltest5l2.CVaR(t,1)=modeltest5l2.W(18,t)+ratio*sum(modeltest5l2.W(19:end,t));
end

modeltest5l2.cumperf=cumprod(1+modeltest5l2.perf);
Names_M3={'Model 3l2','Market','GMVP','MV'};
ReturnModel3l2=[modeltest5l2.perf,marketportfolio.perf,GMVP.perf,MV.perf];
SumStat_M3l2=Portfoliostatistics(ReturnModel3l2,RF,Const.conflvl95,Names_M3);%sumarize portfolio statistics
modeltest5l2TO=turnoverp(modeltest5l2.W(1:17,:),modeltest5l2.perf,Return(rollingwindow:end,:)); % turnonver
TOM3l2=[modeltest5l2TO,'-',GMVP.TO,MV.TO];

f17=figure;
area(date((1+rollingwindow):end),modeltest5l2.W(1:N,:)','FaceColor','flat');
ylim([0,1])
legend(names)
title('max expected return -2ES')
saveas(f17,'results/max expected return -2ES','png');
%% lambda=5 under model 3
Const.lambda=5;
modeltest5l5.W=zeros(N+rollingwindow+1,optimperiod);
modeltest5l5.perf=zeros(optimperiod,1);
modeltest5l5.CVaR=zeros(optimperiod,1);

for t=1:optimperiod
    ReturnInRollingwindow=Return(t:t+rollingwindow-1,:);
    expectedreturn=mean(ReturnInRollingwindow);
    A2=[-ReturnInRollingwindow ,(-1).*ones(rollingwindow,1) ,(-1).*eye(rollingwindow)];
    
    objfun5=@(W) -(expectedreturn*W(1:N,1)-Const.lambda.*(W(N+1,1)+ratio*sum(W(N+2:end,1))));
    modeltest5l5.W(:,t)=fmincon(objfun5,x0,...
        A2,Const.CVaR2,Aeq,Const.sumW,Const.lb,Const.ub,[],options);
    
    modeltest5l5.perf(t,1)=Return(t+rollingwindow,:)*modeltest5l5.W(1:N,t);
    modeltest5l5.CVaR(t,1)=modeltest5l5.W(18,t)+ratio*sum(modeltest5l5.W(19:end,t));
end

modeltest5l5.cumperf=cumprod(1+modeltest5l5.perf);
Names_M3={'Model 3l2','Market','GMVP','MV'};
ReturnModel3l5=[modeltest5l5.perf,marketportfolio.perf,GMVP.perf,MV.perf];
SumStat_M3l5=Portfoliostatistics(ReturnModel3l5,RF,Const.conflvl95,Names_M3);
modeltest5l5TO=turnoverp(modeltest5l5.W(1:17,:),modeltest5l5.perf,Return(rollingwindow:end,:));
TOM3l5=[modeltest5l5TO,'-',GMVP.TO,MV.TO];

f21=figure;
area(date((1+rollingwindow):end),modeltest5l5.W(1:N,:)','FaceColor','flat');
ylim([0,1])
legend(names)
title('max expected return -5ES')
saveas(f21,'results/max expected return -5ES','png');