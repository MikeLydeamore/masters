function r=malthusian_household(k, beta, gamma, alpha)

%alpha=2;
pi_k=1/length(k)*ones(1,length(k));
for i=1:k %For each household size
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
    matricesToSolve=[matricesToSolve A b]; %At end of loop, will have an A and b for each k
end
%Numerically solve
f=@(x) pi_k.*rMat(A
r=fsolve(@(x) rMat(A,x,1,b),ones(noEqn,1));
end

function toSolve(pi_k,

sum(pi_k*rMat(A,x,1,b))