function r=r_fzero(beta,gamma,k,alpha)
[Q, stateList]=genQ(beta,gamma,k);
f=@(r) sysSolve(r,Q,stateList,alpha,k)-1;
r=fzero(f,0.1);


end

function s=sysSolve(r,Q,stateList,alpha,k)
%Solves the system of equations from a matrix Q including a discounting
%rate r which moves to the absorbing class from the transient class

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

s=x'*(initialState==transient);
end