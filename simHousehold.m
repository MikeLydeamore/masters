clear;
gamma=2; beta=40; k=10;
[Q stateList]=genQ(beta,gamma,k);

%Find starting state with 1 infected, rest susceptible
findRoot=roots([1 1 -2*length(stateList)]);
popSize=findRoot'*(findRoot>0)-1;
initialState=intersect(find((stateList(1,:)==(popSize-1))),find(stateList(2,:)==1));

state(1)=initialState;
for i=1:1000
    t(i)=exprnd(-1/Q(state(i),state(i)));
    u=rand();
    %Generate transition probabilities
    probs=cumsum(Q(state(i),:)/-Q(state(i),state(i)));
    for ic=1:length(probs)
        if u<probs(ic)
            %Transition to state ic
            state(i+1)=ic;
            break;
        end
    end
    
    if noInfected(stateList,state(i+1))==0
        break;
    end
end

%Plot infected vs time
times=[0 cumsum(t)];
for i=1:length(state)
    infNo(i)=noInfected(stateList,state(i));
end
plot(times,infNo,'o');
