function gam=rStar(beta,gamma,sigma,ro,tau,onDist,zeta,offDist,kappa,ita,pi_k,alpha)

switch onDist
    case 'exp'
    switch offDist
        case 'exp' %exp on, exp off
            gam=0;
            for k=1:length(pi_k)
                [Q, stateList]=genQ(beta,gamma,sigma,ro,tau,1/zeta,1/kappa,ita,k,0);
                transient=find(diag(Q)~=0);
                
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
                
                gam=gam+(pi_k(k)*x'*(initialState==transient));
            end
            %disp(fullGam);
        case 'const' %exp on, const off
            %Have 3 integrals. From 0-T1~exp, T1-T1+c, T1+c-Inf
            gam=0;
            for k=1:length(pi_k)
                %0-T1
                %Same as running from 0-Inf but with the antiviral never
                %turning off, so use a Q that has kappa=0
                integral1=0;
                [Q, stateList]=genQ(beta,gamma,sigma,ro,tau,1/zeta,0,ita,k,1);
                transient=find(diag(Q)~=0);
                
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
                
                %Want state (k-1, 1, 0, 0)
                initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==1&stateList(3,:)==0&stateList(4,:)==0);
                
                integral1=(pi_k(k)*x'*(initialState==transient));
                
                %T1-T1+c is the same as p(T1) from 0-c
                %First get p(T1). This is just a Q with absorption as soon as
                %AVs are introduced, but the 'lower half' of this system.
                integral2=0;
                [Q, stateList]=genQ(beta,gamma,sigma,ro,tau,1/zeta,0,ita,k,1);
                %Estimate at t=500 that we are in the steady state.
                p0=zeros(1,length(stateList));
                p0(initialState)=1;
                pT1=mexpv(1000,Q',p0);
                lowerStates=(stateList(4,:)==1);
                pT1=pT1(lowerStates);
                
                %Now solve the system of DEs for constant time
                %The Q matrix here is our regular running Q. We reduce by
                %(1-tau)(1-ro)
                [Q, stateList]=genQHalf(beta,gamma,sigma,ro,tau,ita,k);
                transient=find(diag(Q)~=0);
                
                for i=1:length(transient)
                    for j=1:length(transient)
                        A2(i,j)=Q(transient(i),transient(j));
                    end
                    %The external infections happen at rate 
                    %(1-tau)*alpha when there is antivirals.
                    b2(i)=(1-tau)*alpha*stateList(3,transient(i));
                    
                end
                %DE solver
                %psi0 = pT1(transient);
                psi0=zeros(length(b2),1);
                psi = phiv(kappa,A2,b2',psi0);
                integral2=(pi_k(k)*psi'*pT1(transient));
                
                %c - Inf is the same as p(T2) from 0 to Inf
                %First get p(T2). This is p(T1)exp(Qc) where Q is the reduced
                %system.
                integral3=0;
                [Q, stateList]=genQHalf(beta,gamma,sigma,tau,ro,ita,k);
                %p(T2)=p(T1)*exp(Qc)
                pT2=mexpv(kappa,Q',pT1);
					 
					 %Need an unreduced Q for the path integral
					 [Q, stateList]=genQHalf(beta,gamma,sigma,0,0,0,k);
                
                for i=1:length(transient)
                    for j=1:length(transient)
                        A3(i,j)=Q(transient(i),transient(j));
                    end
                    b3(i)=alpha*stateList(3,transient(i));
                end
                
                x=A3\-b3'; %Ax+b = 0
                
                integral3=pi_k(k)*x'*pT2(transient);
                
                gam=gam+integral1+integral2+integral3;
            end
                
            
    end
    case 'const'
        switch offDist
            case 'exp' %const on, exp off
                %We again have 3 integrals, 0-c, c-T2, T2-Inf.
                gam=0;
                for k=1:length(pi_k)
                    %0-c
                    integral1=0;
                    %Note that we have no antivirals until c, so our Q is
                    %the normal unreduced Q.
                    [Q, stateList]=genQHalf(beta,gamma,sigma,0,0,0,k);
                    transient=find(diag(Q)~=0);
                    
                    for i=1:length(transient)
                        for j=1:length(transient)
                            A1(i,j)=Q(transient(i),transient(j));
                        end
                        %The external infections happen at rate alpha
                        b1(i)=alpha*stateList(3,transient(i));
                        
                    end
                    %DE solver
                    %Want state (k-1, 0, 1)
                    initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==0&stateList(3,:)==1);
                    psi0 = zeros(length(b1),1);
                    psi = phiv(zeta,A1,b1',psi0);
                    integral1=(pi_k(k)*psi'*(initialState==transient));
                
                    %c - T2
                    %This is the same as p(c) and then 0-T2
                    %First, p(c). After c time, the household is forced on
                    %with antivirals, so the distribution will just be
                    %p0exp(Q1*c)
                    p0=zeros(1,length(Q));
                    p0(initialState)=1;
                    pc = mexpv(zeta,Q',p0);
                    %Shift this down to the 'middle section' of the large Q
                    %matrix, where dynamics have the reduced rates.
                    pc = [zeros(length(pc),1);pc;zeros(length(pc),1)];
                    
                    %Now, solve for 0-T1 which is the same as 0-Inf but
                    %starting with AVs and never turning them off.
                    integral2=0;
                    [Q, stateList]=genQ(beta,gamma,sigma,ro,tau,0,1/kappa,ita,k,-1);
                    %[Q, stateList]=genQHalf(beta,gamma,sigma,ro,tau,ita,k);
                    transient=find(diag(Q)~=0);
                
                    for i=1:length(transient)
                        for j=1:length(transient)
                            A2(i,j)=Q(transient(i),transient(j));
                        end
                        %The external infections happen at rate alpha when there is no antivrials
                        %and rate (1-tau)*alpha when there is antivirals.
                        if stateList(4,transient(i))==1
                            b2(i)=(1-tau)*alpha*stateList(3,transient(i));
                        else
                            b2(i)=alpha*stateList(3,transient(i));
                        end
                    end
                    x=A2\-b2';
                    
                    integral2=pi_k(k)*x'*pc(transient);
                    
                    %Integral 3
                    %T2-Inf
                    %Same as 0-Inf with p(T2)
                    
                    %Estimate at t=500 that we are in the steady state.
                    [Q, stateList]=genQ(beta,gamma,sigma,ro,tau,0,1/kappa,ita,k,-1);
                    pT2=mexpv(10000,Q',pc);
                    upperStates=stateList(4,:)==2;
                    pT2=pT2(upperStates);
                    if abs(sum(pT2)-1)>1e-3
                        error('p(T2) not 1');
                    end
                    
                    %Now do the integral using the unmodified Q
                    [Q, stateList]=genQHalf(beta,gamma,sigma,0,0,0,k);
                    transient=find(diag(Q)~=0);
                
                    for i=1:length(transient)
                        for j=1:length(transient)
                            A3(i,j)=Q(transient(i),transient(j));
                        end
                            b3(i)=alpha*stateList(3,transient(i));
                    end
                    x=A3\-b3';
                    integral3=pi_k(k)*x'*pT2(transient);
                    
                    gam=gam+integral1+integral2+integral3;
                end
            case 'const' %const on, const off
                %3 Integrals, 0-c, c-d, d-inf
                gam=0;
                for k=1:length(pi_k)
                    
                    %Integral 1
                    integral1=0;
                    %Note that we have no antivirals until c, so our Q is
                    %the normal unreduced Q.
                    [Q, stateList]=genQHalf(beta,gamma,sigma,0,0,0,k);
                    transient=find(diag(Q)~=0);
                    
                    for i=1:length(transient)
                        for j=1:length(transient)
                            A1(i,j)=Q(transient(i),transient(j));
                        end
                        %The external infections happen at rate alpha
                        b1(i)=alpha*stateList(3,transient(i));
                        
                    end
                    %DE solver
                    %Want state (k-1, 0, 1)
                    initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==0&stateList(3,:)==1);
                    psi0 = zeros(1,length(b1));
                    %psi0(transient==initialState)=1;
                    psi = phiv(zeta,A1,b1',psi0');
                    integral1=(pi_k(k)*psi'*(initialState==transient));
                    
                    %Integral 2
                    integral2=0;
                    
                    %First need p(c)
                    p0=zeros(1,length(Q));
                    p0(initialState)=1;
                    pC=mexpv(zeta,Q',p0);
                    
                    %Now integrate from 0-d using the reduced Q
                    
                    [Q, stateList]=genQHalf(beta,gamma,sigma,tau,ro,ita,k);
                    transient=find(diag(Q)~=0);
                    
                    for i=1:length(transient)
                        for j=1:length(transient)
                            A2(i,j)=Q(transient(i),transient(j));
                        end
                        %The external infections happen at rate alpha
                        %b2(i)=alpha*stateList(3,transient(i));
                        b2(i) = (1-tau)*alpha*stateList(3,transient(i));
                        
                    end
                    %DE solver
                    %psi0 = pC(transient);
                    psi0=zeros(length(b2),1);
                    psi = phiv(kappa,A2,b2',psi0);
                    integral2=(pi_k(k)*psi'*pC(transient));
                    
                    %Integral 3
                    integral3=0;
                    %Now we have no antivirals again from d-Inf, so use p(d)
                    %then 0-Inf
                    
                    %Note that we use the reduced Q to get the distribution
                    %here, but our starting point is p(c)
                    
                    pD=mexpv(kappa,Q',pC);
                    
                    %Now integrate from 0-Inf, which involves no DEs.
                    [Q, stateList]=genQHalf(beta,gamma,sigma,0,0,0,k);
                    transient=find(diag(Q)~=0);
                    
                    for i=1:length(transient)
                        for j=1:length(transient)
                            A3(i,j)=Q(transient(i),transient(j));
                        end
                        b3(i)=alpha*stateList(3,transient(i));
                    end
                    x=A3\-b3';
                    integral3=pi_k(k)*x'*pD(transient);
                    
                    gam=gam+integral1+integral2+integral3;
                end
        end
end
end