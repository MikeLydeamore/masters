function gam=expecGamma(Q,stateList)

noStates=size(Q);
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
        A(i,j)=Q(transient(i),transient(j));
    end
    b(i)=stateList(2,i);
end
disp(A);
disp(b);

gam=A\b;
end