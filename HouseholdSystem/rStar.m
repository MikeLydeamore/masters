function r_star=rStar(k,beta,gamma,alpha)
%Solves the inter-household reproduction number R_* for households of
%size k
%If k is a vector, then the solution will be for the households of hetero
%size k(1), k(2), ..., k(end)
%Currently assumes each household size is equally likely
pi_k=1/length(k)*ones(1,length(k));
r_star=0;
for i=1:length(k)
    [Q stateList]=genQ(beta,gamma,k);
    gam=expecGamma(Q,stateList,alpha);
    r_star=r_star+(pi_k*gam);
end
    

