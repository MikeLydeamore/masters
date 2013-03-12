function r=malthusian(Q,stateList)

alpha=2;

noStates=length(Q);
remainingStates=1:noStates;
%Determine 0 rows (absorbing states)
absorbing=find(diag(Q)==0)';
transient=remainingStates;
transient(absorbing)=[];

noEqn=length(transient);
A=zeros(noEqn);
b=zeros(noEqn,1);

for i=1:noEqn
    for j=1:noEqn
        A(i,j)=Q(transient(i),transient(j)); %Build system to solve
    end
    b(i)=alpha*stateList(2,transient(i));
end

%Numerically solve
r=fsolve(@(x) rMat(A,x,1,b),[1;1;1]);
end


