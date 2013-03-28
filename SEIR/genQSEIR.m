function [Q stateList]=genQSEIR(beta,gamma,k)

% Hideously ugly, but does work.
counter=1;
for s=1:k+1
    for e=1:k+1
        for i=1:k+1
            if (s-1)+(e-1)+(i-1)<=k
                stateList(:,counter)=[s-1;e-1;i-1];
                counter=counter+1;
            end
        end
    end
end
[m, noStates]=size(stateList);

Q=sparse(noStates,noStates);

%Determine infection events
%Require S>0, I>0
infectionEvents=intersect(find(stateList(1,:)>0),find(stateList(3,:)>0));

%Move to state (S-1, E+1, I) at rate beta*S*I/(k-1)
for i=1:length(infectionEvents)
    moveTo=intersect(intersect(find(stateList(1,:)==stateList(1,infectionEvents(i))-1),find(stateList(2,:)==stateList(2,infectionEvents(i))+1)),find(stateList(3,:)==stateList(3,infectionEvents(i))));
    Q(infectionEvents(i),moveTo)=beta*stateList(1,infectionEvents(i))*stateList(3,infectionEvents(i))/(k-1);
end

%Determine shedding events
%Require E>0

sheddingEvents=find(stateList(2,:)>0);

%Move to state (S, E-1, I+1) at rate sigma*E
for i=1:length(sheddingEvents)
end