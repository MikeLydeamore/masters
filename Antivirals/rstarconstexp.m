function rs=rstarconstexp(beta,gamma,sigma,ro,tau,onDist,zeta,offDist,kappa,ita,N,alpha)

%Here we only have 2 integrals.

%0-c

%System of DE's:
%No antivirals present
[Q, stateList]=genQHalf(beta,gamma,sigma,0,0,0,N);

%Initial state is (N-1,1,0)
initialstate=find(stateList(1,:)==N-1&stateList(2,:)==1&stateList(3,:)==0);

transient=find(diag(Q)~=0);

for i=1:length(transient)
    for j=1:length(transient)
        A1(i,j)=Q(transient(i),transient(j));
    end
    %The external infections happen at rate alpha
    b1(i)=alpha*stateList(3,transient(i));
    
end
psi0=zeros(length(b1),1);
psi = phiv(zeta,A1,b1',psi0);

integral1=psi'*(transient==initialstate);

%c-T1+c, T1+c-Inf is just c-Inf, which is just p(c) 0-Inf using the same
%methodology as the exponentially distributed delay in the Supp Material.

%p(c)
p0=zeros(1,length(Q));
p0(initialstate)=1;
pc=mexpv(10000,Q',p0);
%This corresponds to the starting distribution, even with the larger Q we
%have below.

%Now, set up the Q
[Q1, stateList]=genQHalf(beta,gamma,sigma,ro,tau,ita,N);

switchblock=(1/kappa)*eye(length(Q));


Q2=genQHalf(beta,gamma,sigma,0,0,0,N);

Q=[Q1 switchblock ; zeros(size(switchblock)) Q2];
stateList=[stateList stateList ; zeros(1,length(stateList)) ones(1,length(stateList))];

transient=find(diag(Q)~=0);

for i=1:length(transient)
    for j=1:length(transient)
        A2(i,j)=Q(transient(i),transient(j));
    end
    %The external infections happen at rate alpha or (1-tau)alpha
    if stateList(4,transient(i))==0
        b2(i)=alpha*stateList(4,transient(i));
    else
        b2(i)=(1-tau)*alpha*stateList(4,transient(i));
    end
end
x=A2\-b2';

%Note that our distribution p(c) has to be extended with some zeros in the
%places where the antivirals have stopped (The lowest block of Q).
pc=[pc zeros(length(pc))];
integral2=x'*pc(transient);

rs=integral1+integral2;
    
end