function gam=expecGamma(Q,stateList,alpha)

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

gam=fullGam(transientState);
disp(fullGam);
end