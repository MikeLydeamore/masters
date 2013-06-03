function [avgInfected, t]=SEIR_Q_Sim(beta,gamma,sigma,k,tspan,graph)

[Q, stateList]=genQSEIR(beta,gamma,sigma,k);
[noStates, discard]=size(Q);
%Initial state is (k-1,0,1)
initialState=intersect(intersect(find(stateList(1,:)==k-1),find(stateList(2,:)==0)),find(stateList(3,:)==1));

t=tspan(1):0.01:tspan(end);
p0=zeros(1,noStates);
p0(initialState)=1;

avgInfected=zeros(1,length(t));
p=zeros(length(t),length(p0));
for i=1:length(t)
    p(i,:)=mexpv(t(i),Q',p0);
    
    %Get number of infected
    avgInfected(i)=stateList(3,:)*p(i,:)';
end

if graph==1
    plot(t,avgInfected);
end


end