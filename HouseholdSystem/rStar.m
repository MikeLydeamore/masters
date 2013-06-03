function r_star=rStar(k,beta,gamma,alpha,distribution)
%Solves the inter-household reproduction number R_* for households of
%size k
%If k is a vector, then the solution will be for the households of hetero
%size k(1), k(2), ..., k(end)
%Optional 5th argument 'distribution' which allows for a predefined
%distribution of household sizes. Must be the same length as k

if nargin<5
    pi_k=1/length(k)*ones(1,length(k));
else
    pi_k=distribution;
end
r_star=0;
for i=1:length(k)
    [Q stateList]=genQ(beta,gamma,k(i));
    gam=expecGamma(Q,stateList,alpha);
    r_star=r_star+(pi_k(i)*gam);
end
    

