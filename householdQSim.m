function expectedInfected=householdQSim(beta,gamma,k)
%gamma=2; beta=20; k=4;

[Q stateList]=genQ(beta,gamma,k);
%Find starting state with 1 infected, rest susceptible
findRoot=roots([1 1 -2*length(stateList)]);
popSize=findRoot'*(findRoot>0)-1;
initialState=intersect(find((stateList(1,:)==(popSize-1))),find(stateList(2,:)==1));

p_0=zeros(1,length(stateList));
p_0(initialState)=1; %Set to initial infected

t=0:0.001:5;
expectedInfected=zeros(1,length(t));
for i=1:length(t)
    pt=mexpv(t(i),Q',p_0');
    expectedInfected(i)=stateList(2,:)*pt;
end
%Testing GIT Upload
plot(t,expectedInfected,'r')
end
