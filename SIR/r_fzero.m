function r=r_fzero(beta,gamma,alpha,pi_k)

%[Q, stateList]=genQ(beta,gamma,k);
f=@(r) sysSolve(r,beta,gamma,alpha,pi_k);
x0=2*length(pi_k);
r=fzero(f,1);

end

function s=sysSolve(r,beta,gamma,alpha,pi_k)
%Solves the system of equations from a matrix Q including a discounting
%rate r which moves to the absorbing class from the transient class

s=0;
for k=1:length(pi_k)
    [Q stateList]=genQ(beta,gamma,k);
    
    transient=find(diag(Q)~=0);
    for i=1:length(transient)
        Q(transient(i),transient(i))=Q(transient(i),transient(i))-r; %Add discounting rate, -(qii + r) = -qii - r
    end
    
    for i=1:length(transient)
        for j=1:length(transient)
            A(i,j)=Q(transient(i),transient(j));
        end
        b(i)=alpha*stateList(2,transient(i));
    end
    
    x=A\-b'; %Ax+b = 0
    
    initialState=intersect(find(stateList(1,:)==(k-1)),find(stateList(2,:)==1));
    
    %Weight by pi_k
    s=s+pi_k(k)*(x'*(initialState==transient));
end
s=s-1; %Want sum(pi_kE[s])=1
end