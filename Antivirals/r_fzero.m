function r=r_fzero(beta,gamma,sigma,ro,tau,onDist,zeta,offDist,kappa,ita,pi_k,alpha)
%Solves for the malthusian r for an SEIR system without antivirals
%[Q, stateList]=genQ(beta,gamma,sigma,ro,tau,zeta,pi_k);
switch onDist
    case 'exp'
        switch offDist
            case 'exp'
                f=@(r) sysSolve(r,beta,gamma,sigma,ro,tau,zeta,kappa,ita,alpha,pi_k);
            case 'const'
					f=@(r) r_expconst(r,beta,gamma,sigma,ro,tau,zeta,kappa,ita,alpha,pi_k);
            otherwise
                error('Unknown type of offDist');
        end
    case 'const'
        switch offDist
            case 'exp'
					
            case 'const'
            otherwise
                error('Unknown type of offDist');
        end
    otherwise
        error('Unknown type of onDist');
end
r=fzero(f,[0 5]);

end

function s=sysSolve(r,beta,gamma,sigma,ro,tau,zeta,kappa,ita,alpha,pi_k)
%Solves the system of equations from a matrix Q including a discounting
%rate r which moves to the absorbing class from the transient classs
s=0;
for k=1:length(pi_k)
    [Q, stateList]=genQ(beta,gamma,sigma,ro,tau,zeta,kappa,ita,k);
    transient=find(diag(Q)~=0);
    for i=1:length(transient)
        Q(transient(i),transient(i))=Q(transient(i),transient(i))-r; %Add discounting rate, -(qii + r) = -qii - r
    end
    
    for i=1:length(transient)
        for j=1:length(transient)
            A(i,j)=Q(transient(i),transient(j));
        end
        %The external infections happen at rate alpha when there is no antivrials
        %and rate (1-tau)*alpha when there is antivirals.
        if stateList(4,transient(i))==1
            b(i)=(1-tau)*alpha*stateList(3,transient(i));
        else
            b(i)=alpha*stateList(3,transient(i));
        end
    end
    
    x=A\-b'; %Ax+b = 0
    
    %Want state (k-1, 1, 0, 0)
    initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==1&stateList(3,:)==0&stateList(4,:)==0);
    
    s=s+(pi_k(k)*x'*(initialState==transient));
end
s=s-1;
end