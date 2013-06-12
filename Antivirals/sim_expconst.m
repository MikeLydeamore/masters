function [infected, eventTime, noAV, totalInfected, popSize]=sim_expconst(beta,gamma,sigma,alpha,zeta,kappa,ro,tau,ita,pi_k,phi_k,N)

global testIterations;
%Attempt at discrete event simulator

t=0;
maxTime=1000;
counter=1;
noAV=0;

%Set up households
housePop=zeros(1,length(pi_k));
for i=1:length(pi_k)-1
	housePop(i)=round(pi_k(i)*N);
end
housePop(end)=N-sum(housePop);
popSize=housePop*(1:length(housePop))';

%Load population
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

AVTime=Inf*ones(N,1);

while t<maxTime && sum(popStatus(:,2)+popStatus(:,3))>0

	infected(counter)=sum(popStatus(:,3));
	eventTime(counter)=t;
	noAV=noAV+(popStatus(:,5)==1)*sum(popStatus(1:4),2);
	
	if mod(counter,100)==0
		clc
		fprintf('Current run: %d \n',testIterations);
		fprintf('Current time: %.4f \n Current Infected: %d \n',t,infected(counter));
	end
	%Update rates
	
	%Internal infection:
	modifiers=(1-popStatus(:,5).*tau).*(1-popStatus(:,5).*ro);
	
	rates(1:N)=beta./(sum(popStatus,2)-1).*popStatus(:,1).* popStatus(:,3).*modifiers;
	
	%External infection: (Add in)
	InfWithAV=popStatus(:,5)'*popStatus(:,3);
	InfWithoutAV=not(popStatus(:,5))'*popStatus(:,3);
	tfi=alpha*((1-tau)*InfWithAV+InfWithoutAV)./popSize;
	
	rates(1:N)=rates(1:N)+(tfi*popStatus(:,1).*(1-popStatus(:,5)*ro))';
	
	%Progression:
	rates(N+1:2*N)=sigma*popStatus(:,2);
	
	%Recovery
	rates(2*N+1:3*N)=(1+popStatus(:,5)*ita).*gamma.*popStatus(:,3);
	
	%Antiviral Introduction:
	rates(3*N+1:4*N)=zeta*((popStatus(:,3)+(popStatus(:,4)==1))>0).*not(popStatus(:,5));
	
	
	%Check what happens next:
	[next_const_time, household]=min(AVTime>0);
	
	sumRates=sum(rates);
	normRates=rates./sumRates;
	%Choose event:
	u=rand();
	eventOccured=find(u<cumsum(normRates),1);
	next_other_time=exprnd2(1/sumRates);
	
	if next_const_time < next_other_time
		%Next event is antivirals being removed
		t=t+next_const_time;
		popStatus(household,5)=2;
		AVTime(household)=Inf;

	else
		%Next event is... something else.
		t=t+next_other_time;
		%Determine event and household
		household=mod(eventOccured,N);
		if household==0
			household=N;
		end
		if eventOccured<N+1 %Infection
			popStatus(household,1)=popStatus(household,1)-1; %S -> S-1
			popStatus(household,2)=popStatus(household,2)+1; %E -> E+1
			
		elseif eventOccured < 2*N+1 %Progression
			popStatus(household,2)=popStatus(household,2)-1; %E -> E-1
			popStatus(household,3)=popStatus(household,3)+1; %I -> I+1
			
		elseif eventOccured < 3*N+1 %Recovery
			popStatus(household,3)=popStatus(household,3)-1; %I -> I-1
			popStatus(household,4)=popStatus(household,4)+1; %R -> R+1
			
		elseif eventOccured < 4*N+1 %Antivirals introduced
			popStatus(household,5)=1;
			AVTime(household)=t+kappa; %Set the time of activation.
			
		else
			error('Unknown event');				
		end
	end
	counter=counter+1;
end
eventTime(counter)=t;
totalInfected=sum(infected);
end