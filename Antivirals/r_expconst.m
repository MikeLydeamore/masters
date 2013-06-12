function s=r_expconst(r,beta,gamma,sigma,ro,tau,zeta,kappa,ita,alpha,pi_k)


%Break up into 3 Integrals
s=0;

for k=1:length(pi_k)
	%Integral 1
	%0-T where T is ~exp(zeta).
	%Use the Q matrix with absorption after antivirals introduced.
	[Q, stateList]=genQ(beta,gamma,sigma,ro,tau,zeta,kappa,ita,k,1);
	transient=find(diag(Q)~=0);
	
	%Q(transient,transient)=Q(transient,transient)-r;
	for i=1:length(transient)
		Q(transient(i),transient(i))=Q(transient(i),transient(i))-r; %Add discounting rate, -(qii + r) = -qii - r
	end
	for i=1:length(transient)
		for j=1:length(transient)
			A1(i,j)=Q(transient(i),transient(j));
		end
		%The external infections happen at rate alpha when there is no antivrials
		%and rate (1-tau)*alpha when there is antivirals.
		if stateList(4,transient(i))==1
			b1(i)=(1-tau)*alpha*stateList(3,transient(i));
		else
			b1(i)=alpha*stateList(3,transient(i));
		end
	end
	
	x=A1\-b1'; %Ax+b = 0
	
	initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==0&stateList(3,:)==1&stateList(4,:)==0);
	integral1=x'*(initialState==transient);
	
	%Integral 2
	%First need p(T). This is just the absorbing matrix run for a long time.
	[Q, stateList]=genQ(beta,gamma,sigma,ro,tau,zeta,kappa,ita,k,1);
	p0=zeros(length(stateList),1);
	p0(initialState)=1;
	pT=mexpv(10000,Q',p0);
	middlestates= stateList(4,:)==1;
	pT=pT(middlestates);
	%Also need to premultiply the remaining 2 integrals by z/(r+z).
	
	%Integral is now from 0 to A, and we solve this using phiv.
	[Q, stateList]=genQHalf(beta,gamma,sigma,ro,tau,ita,k);
	transient=find(diag(Q)~=0);
	
	%Q(transient,transient)=Q(transient,transient)-r;
	for i=1:length(transient)
		Q(transient(i),transient(i))=Q(transient(i),transient(i))-r; %Add discounting rate, -(qii + r) = -qii - r
	end
	for i=1:length(transient)
		for j=1:length(transient)
			A2(i,j)=Q(transient(i),transient(j));
		end
		%There is always antivirals for this integral.
		b2(i)=(1-tau)*alpha*stateList(3,transient(i));
	end
	
	%Our initial state for the DE's is all zeros.
	psi0=zeros(length(b2),1);
	psi = phiv(kappa,A2,b2',psi0);
	integral2=(zeta/(r+zeta))*pT(transient)'*psi;
	
	
	%Integral 3
	%Need p(A).
	[Q, stateList]=genQHalf(beta,gamma,sigma,ro,tau,ita,k);
	pA=mexpv(kappa,Q',pT);
	%Also need to multiply this integral by exp(-rA)
	
	[Q, stateList]=genQHalf(beta,gamma,sigma,0,0,0,k);
	%Q(transient,transient)=Q(transient,transient)-r;
	for i=1:length(transient)
		Q(transient(i),transient(i))=Q(transient(i),transient(i))-r; %Add discounting rate, -(qii + r) = -qii - r
	end
	for i=1:length(transient)
		for j=1:length(transient)
			A3(i,j)=Q(transient(i),transient(j));
		end
		%There is never always antivirals for this integral.
		b3(i)=alpha*stateList(3,transient(i));
	end
	
	x=A3\-b3';
	
	integral3=(zeta/(r+zeta))*(exp(-r*kappa))*pA(transient)'*x;
	
	
	s=s+(pi_k(k)*(integral1+integral2+integral3));
end
s=s-((r+sigma)/sigma);