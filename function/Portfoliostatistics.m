function SumStat = Portfoliostatistics(Return,Rf,conflvl,names)
[T,N]=size(Return);
RF=mean(Rf)*1200;
portfoliomean=(prod(1+Return).^(12/T)-1).*100;
portfoliovol=sqrt(var(Return).*12).*100;
sharpe=(portfoliomean-RF)./portfoliovol;
CVaR=NaN(1,N);
for i=1:N
    CVaR(1,i)=cvarp(Return(:,i),conflvl).*100;
end

name=names;
statistics=[portfoliomean;portfoliovol;sharpe;CVaR];
SumStat=mat2dataset(statistics,'VarNames',name,'ObsNames',...
    {'Average','Volatility','Sharpe Ratio','CVaR'});
SumStat=dataset2table(SumStat);
end

%