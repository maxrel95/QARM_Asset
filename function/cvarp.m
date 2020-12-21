function [CVaR] = cvarp(Return,conflvl)
%INPUT: Return: vector of portfolio return
%       conflvl: level of confidence eg 0.95, 0.99
%OUTPUT: CVaR: give the historical CVaR of the portfolio 
VaR=quantile(Return,1-conflvl);
CVaR=-mean(Return(Return<VaR));
end

