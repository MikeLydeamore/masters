function [currentInfected, eventTime, noAV, totalInfected, popSize]=simulateAV(beta,gamma,sigma,alpha,onDist,zeta,offDist,kappa,ro,tau,ita,pi_k,phi_k,N)
%Simulates the Markov chain where antivirals are introduced.
%beta: Internal Infection Rate
%gamma: Recovery Rate
%sigma: Progression Rate
%alpha: External Infection Rate
%onDist: Distribution of introduction rate (const or exp)
%zeta: Antiviral Introduction Rate
%offDist: Distribution of removal rate (const or exp)
%kappa: Constant time off for AVs
%Ro: Reduction in susecptibility due to antivirals
%Tau: Reduction in infectivity due to antivirals
%Ita: Reduction in infectious period due to antivirals
%pi_k: Distribution of household sizes
%phi_k: Distribution of pre-allocation to household sizes (Must be the same length as pi_k)
%N: Number of households
%
%Note that the reduction in infectious period is equivalent to a faster
%recovery rate.

%Distribute household sizes (must add up to N)
housePop=zeros(1,length(pi_k));
for i=1:length(pi_k)-1
	housePop(i)=round(pi_k(i)*N);
end
housePop(end)=N-sum(housePop);
popSize=housePop*(1:length(housePop))';

%Maximum iterations is now 4*popSize+N, as there is a transition in each
%household where they get the antivirals.
maxIter=4*popSize+N;

%Load in the population
popStatus=zeros(N,5); %S-E-I-R + Antiviral Status

cumHouseSizes=[0 cumsum(housePop)];
for i=2:length(cumHouseSizes)
	popStatus(cumHouseSizes(i-1)+1:cumHouseSizes(i),1)=i-1;
end
init_inf=20; %I_0
popStatus(1:init_inf,1)=popStatus(1:init_inf,1)-1; %Introduce I0 infected
popStatus(1:init_inf,3)=1;

%Preallocate the antivirals
%phi_k(k) is the proportion of households of size k that are preallocated AVs.
preAllocation=round(phi_k.*housePop);
%Load in antivirals
cumSizes=[0 cumsum(housePop)];
for i=1:length(cumSizes)-1
	popStatus(cumSizes(i)+1:cumSizes(i)+preAllocation(i),5)=1;
end

eventTime=zeros(1,maxIter);

currentInfected=zeros(1,maxIter);
currentExposed=zeros(1,maxIter);
noAV=0;

switch onDist
case 'exp'
switch offDist
case 'exp'
%Rate vector goes [Infection ... Progression ... Recovery ... A.V. Intro A.V. Removal]
rates=zeros(5*N,1);
for i=1:maxIter
%     fprintf('Iteration %d of %d \n',i,maxIter);
%     fprintf('Current elapsed time: %.4f\n',timer);
%     fprintf('Expected completion time: %.4f\n',(timer/i)*maxIter);
	currentInfected(i)=sum(popStatus(:,3));
	currentExposed(i)=sum(popStatus(:,2));
	if currentInfected(i)==0 && currentExposed(i)==0 %Absorbing state
		break;
	end
			

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           Rates                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Infection - Internal
    %If household has antivirals, then beta will be modified by
    %(1-tau)(1-ro)
	modifiers = (1-popStatus(:,5).*tau).*(1-popStatus(:,5).*ro);
    
	rates(1:N)=beta./(sum(popStatus,2)-1).*popStatus(:,1).* popStatus(:,3).*modifiers;
    
	%Infection - External
	%If the external household has antivirals, then we have (1-tau)
	%If the internal household has antivirals, then we have (1-ro)
	%Need the # of infected from houses with antirivals and without
	%Then the external force is alpha*(1-tau)*I_a + alpha*I_{no a}
	noInfWithAV=popStatus(:,5)'*popStatus(:,3);
	noInfWithoutAV=not(popStatus(:,5))'*popStatus(:,3);
    
	tfi=alpha*((1-tau)*noInfWithAV+noInfWithoutAV)./popSize;
	rates(1:N)=rates(1:N)+tfi*popStatus(:,1).*(1-popStatus(:,5)*ro);
       
	%Latent Progression
	%Progression always at the same rate, sigma*E
	rates(N+1:2*N)=sigma*popStatus(:,2);
    
	%Recovery
	%Recovery is at rate gamma*I if there's no antivirals, and
	%(1+ita)*gamma*I if there is antivirals
	rates(2*N+1:3*N)=(1+popStatus(:,5)*ita).*gamma.*popStatus(:,3);
    
	%Antiviral introduction
	%Constant rate zeta into households who have at least 1 infection event
	%and don't already have antivirals
	rates(3*N+1:4*N)=zeta*((popStatus(:,3)+(popStatus(:,4)==1))>0).*not(popStatus(:,5));
    
    %Antiviral removal
    %Constant rate kappa for households who already have the antivirals
    rates(4*N+1:5*N)=kappa*(popStatus(:,5));
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%                         Choosing Events                             %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	sumRates=sum(rates);
	%Normalise
	normRateVector=rates./sumRates;
    
	%Choose event
	u=rand();
	eventOccured=find(u<cumsum(normRateVector),1);
    
	%Calculate event time
	eventTime(i+1)=exprnd2(1/sumRates);
    
	%Determine event and household
	household=mod(eventOccured,N);
	if household==0
		household=N;
	end
	if eventOccured<N+1 %Infection
		event=1;
	elseif eventOccured < 2*N+1 %Progression
		event=2;
	elseif eventOccured < 3*N+1 %Recovery
		event=3;
    elseif eventOccured < 4*N+1 %Antivirals introduced
		event=4;
    else
        event=5;
	end
    
	%State change
	switch event
		case 1 %Infection
			popStatus(household,1)=popStatus(household,1)-1; %S -> S-1
			popStatus(household,2)=popStatus(household,2)+1; %E -> E+1
        
		case 2 %Progression
			popStatus(household,2)=popStatus(household,2)-1; %E -> E-1
			popStatus(household,3)=popStatus(household,3)+1; %I -> I+1
        
		case 3 %Recovery
			popStatus(household,3)=popStatus(household,3)-1; %I -> I-1
			popStatus(household,4)=popStatus(household,4)+1; %R -> R+1
        
        case 4 %AV Introduced
			popStatus(household,5)=1;
            
        otherwise
            popStatus(household,5)=0;
		end
		noAV=noAV+popStatus(:,5)'*sum(popStatus(:,1:4),2);
end

case 'const' %Exp on, const off
    %Rate vector goes [Infection ... Progression ... Recovery ... A.V. Intro]
    rates=zeros(5*N,1);
    AVTime=zeros(N,1); %Time that antiviral has been present in each household.
for i=1:maxIter
%     fprintf('Iteration %d of %d \n',i,maxIter);
%     fprintf('Current elapsed time: %.4f\n',timer);
%     fprintf('Expected completion time: %.4f\n',(timer/i)*maxIter);
	currentInfected(i)=sum(popStatus(:,3));
	currentExposed(i)=sum(popStatus(:,2));
	if currentInfected(i)==0 && currentExposed(i)==0 %Absorbing state
		break;
	end
			

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           Rates                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Infection - Internal
    %If household has antivirals, then beta will be modified by
    %(1-tau)(1-ro)
	modifiers = (1-popStatus(:,5).*tau).*(1-popStatus(:,5).*ro);
    
	rates(1:N)=beta./(sum(popStatus,2)-1).*popStatus(:,1).* popStatus(:,3).*modifiers;
    
	%Infection - External
	%If the external household has antivirals, then we have (1-tau)
	%If the internal household has antivirals, then we have (1-ro)
	%Need the # of infected from houses with antirivals and without
	%Then the external force is alpha*(1-tau)*I_a + alpha*I_{no a}
	noInfWithAV=popStatus(:,5)'*popStatus(:,3);
	noInfWithoutAV=not(popStatus(:,5))'*popStatus(:,3);

	tfi=alpha*((1-tau)*noInfWithAV+noInfWithoutAV)./popSize;
	rates(1:N)=rates(1:N)+tfi*popStatus(:,1).*(1-popStatus(:,5)*ro);
       
	%Latent Progression
	%Progression always at the same rate, sigma*E
	rates(N+1:2*N)=sigma*popStatus(:,2);
    
	%Recovery
	%Recovery is at rate gamma*I if there's no antivirals, and
	%(1+ita)*gamma*I if there is antivirals
	rates(2*N+1:3*N)=(1+popStatus(:,5)*ita).*gamma.*popStatus(:,3);
    
	%Antiviral introduction
	%Constant rate zeta into households who have at least 1 infection event
	%and don't already have antivirals
	rates(3*N+1:4*N)=zeta*((popStatus(:,3)+(popStatus(:,4)==1))>0).*not(popStatus(:,5));
    
    %Antiviral removal
    %Constant rate kappa for households who already have the antivirals
    rates(4*N+1:5*N)=kappa*(popStatus(:,5));
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%                         Choosing Events                             %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	sumRates=sum(rates);
	%Normalise
	normRateVector=rates./sumRates;
    
	%Choose event
	u=rand();
	eventOccured=find(u<cumsum(normRateVector),1);
    
	%Calculate event time
	eventTime(i+1)=exprnd2(1/sumRates);
    
	%Determine event and household
	household=mod(eventOccured,N);
	if household==0
		household=N;
	end
	if eventOccured<N+1 %Infection
		event=1;
	elseif eventOccured < 2*N+1 %Progression
		event=2;
	elseif eventOccured < 3*N+1 %Recovery
		event=3;
    else % eventOccured < 4*N+1 %Antivirals introduced
		event=4;
	end
    
	%State change
	switch event
		case 1 %Infection
			popStatus(household,1)=popStatus(household,1)-1; %S -> S-1
			popStatus(household,2)=popStatus(household,2)+1; %E -> E+1
        
		case 2 %Progression
			popStatus(household,2)=popStatus(household,2)-1; %E -> E-1
			popStatus(household,3)=popStatus(household,3)+1; %I -> I+1
        
		case 3 %Recovery
			popStatus(household,3)=popStatus(household,3)-1; %I -> I-1
			popStatus(household,4)=popStatus(household,4)+1; %R -> R+1
        
        case 4 %AV Introduced
			popStatus(household,5)=1;
            
        otherwise
            popStatus(household,5)=0;
    end
        %Check antiviral times and remove where necessary.
        AVTime=AVTime+(eventTime*popStatus(:,5));
        turnoffAV= AVTime>kappa&popStatus(:,5)==1;
		noAV=noAV+popStatus(:,5)'*sum(popStatus(:,1:4),2);
        popStatus(turnoffAV,5)=0;
end    


end
end
totalInfected=popSize-sum(popStatus(:,1));
end
