function std = portfoliostd(W,sigma)
varportfolio=W'*sigma*W;
std=sqrt(varportfolio);
end

