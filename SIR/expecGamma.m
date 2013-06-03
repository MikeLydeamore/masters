function [gam, fullGam]=expecGamma(beta,gamma,alpha,pi_k)
%Solves for R_star over a distribution of household sizes pi_k.
%pi_k must be a probability vector which represents the probability that a
%randomly chosen person will belong to a household of size k.

gam=0;
for k=1:length(pi_k)
    [Q, stateList]=genQ(beta,gamma,k);
    
    noStates=length(Q);
    remainingStates=1:noStates;
    %Determine 0 rows (absorbing states)
    absorbing=find(diag(Q)==0)';
    transient=remainingStates;
    transient(absorbing)=[];
    
    
    %Set up system to solve
    noEqn=length(transient);
    A=zeros(noEqn);
    b=zeros(noEqn,1);
    
    for i=1:noEqn
        for j=1:noEqn
            A(i,j)=Q(transient(i),transient(j)); %Build system to solve
        end
        b(i)=-alpha*stateList(2,transient(i));
    end
    
    fullGam=A\b;
    
    %If we insist on one infected, we must start in state (k-1,1)
    %This is the gamma we want
    %This method will work even if the stateList changes
    findRoot=roots([1 1 -2*length(stateList)]);
    popSize=findRoot'*(findRoot>0)-1;
    initialState=intersect(find((stateList(1,:)==(popSize-1))),find(stateList(2,:)==1));
    
    transientState=find(transient==initialState);
    
    %R* = sum_k pi_k*E[Gamma|X(0)=1]
    gam=gam+pi_k(k)*fullGam(transientState);
end
%disp(fullGam);
end