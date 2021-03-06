function [Q stateList]=genQ(beta,gamma,k)
%Very naive way of building Q matrix

if nargin < 3
    error('Not enough input arguments');
end
%# of states = (k+1)(k+2)/2
Q=sparse(((k+1)*(k+2))/2,((k+1)*(k+2))/2);

%Generate state list
Infec=1;
S=zeros(1,k+1);
I=[];
for i=k+1:-1:1
    S=[S Infec*ones(1,i-1)];
    Infec=Infec+1;
    I=[I 0:i-1];
end
   

noStates=length(S);
t=[S;I];
%Determine when infection occurs
%Require I>0 and S>0
possibleInfections=intersect(find(t(1,:)),find(t(2,:)));

for i=1:length(possibleInfections) %Determine where the transition is to for infections
    %Searching for state (S-1,I+1)
    currentS=t(1,possibleInfections(i));
    currentI=t(2,possibleInfections(i));
    infectionGoesTo=intersect(find(t(1,:)==(currentS-1)),find(t(2,:)==currentI+1)); %Finds state (S-1,I+1)
    
    %Add infection event to Q matrix, at rate beta*S*I/(k-1)
    Q(possibleInfections(i),infectionGoesTo)=beta*currentS*currentI/(k-1);
end

%Determine when recovery occurs
%Require I>0
possibleRecovery=find(t(2,:));

for i=1:length(possibleRecovery)
    currentS=t(1,possibleRecovery(i));
    currentI=t(2,possibleRecovery(i));
    recoveryGoesTo=intersect(find(t(1,:)==currentS),find(t(2,:)==(currentI-1)));
    
    %Add recovery event to Q matrix, at rate I*gamma
    Q(possibleRecovery(i),recoveryGoesTo)=currentI*gamma;
    
end

%Add in negatives
rowSum=sum(Q,2);
for i=1:noStates
    Q(i,i)=-rowSum(i);
end

stateList=t;

end