function [totalInfected, eventTime]=simulateSEIR(beta,gamma,alpha,sigma,N,k)
%profile on
global testIterations
maxIter=4*k*N;
popStatus=zeros(N,4); %S-E-I-R
popStatus(2:end,1)=k*ones(N-1,1); %Sets N-1 houses to have k susceptibles
init_inf=20;
popStatus(1:init_inf,1)=k-1; %1 infected in 10 of the households
popStatus(1:init_inf,3)=1;

eventTime=zeros(1,maxIter);

%Rate vector will go [Infection Infection ... Progression Progression ... Recovery Recovery ...] etc
rates=zeros(3*N,1); %One row per household

popSize=N*k;

totalInfected=zeros(1,maxIter);
totalExposed=zeros(1,maxIter);

timer=0;
for i=1:maxIter
    tic;
    clc
    fprintf('r_compare iteration: %d\n',testIterations);
    fprintf('Simulation iteration: %d out of %d\n',i,maxIter);
    fprintf('Time elapsed: %.4f \n',timer);
    estimatedtime=(timer/i)*maxIter;
    fprintf('Estimated completion: %.4f \n',estimatedtime);
    %% Rates
    totalInfected(i)=sum(popStatus(:,3));
    totalExposed(i)=sum(popStatus(:,2));
    if totalInfected(i)==0 && totalExposed(i)==0 %No infected left and no exposed left, in absorbing state
        break;
    end
    %Infection (both internal and external) at rate B*I*S/k-1+alpha*totalI/N*S
    rates(1:N)=(beta/(k-1))*popStatus(:,3).*popStatus(:,1)+alpha.*(totalInfected(i)/popSize).*popStatus(:,1);
    
    %Latent progression at rate sigma*E
    rates(N+1:2*N)=sigma*popStatus(:,2);
    
    %Recovery at rate gamma*I
    rates(2*N+1:3*N)=gamma*popStatus(:,3);
    
    %% Reshape, normalise and choose which event occurs
    
    %[m, n]=size(rates);
    %rateVector=reshape(rates,1,m*n);
    
    normRateVector=rates./sum(rates);
    
    %Choose event
    u=rand();
    eventOccured=find(u<cumsum(normRateVector),1);
    
    %Calculate event time
    eventTime(i+1)=exprnd(1/sum(rates));
    
    %Determine event and household
    household=mod(eventOccured,N);
    if household==0
        household=N;
    end
    if eventOccured < N+1 %Infection
        event=1;
    elseif eventOccured < 2*N+1 %Progression
        event=2;
    else
        event=3;
    end
%       [household, event]=ind2sub([N, 4],eventOccured);  
    
    %% Make corresponding change
    
    if event==1 %Infection
        popStatus(household,1)=popStatus(household,1)-1; %S -> S-1
        popStatus(household,2)=popStatus(household,2)+1; %E -> E+1
        
    elseif event==2 %Progression
        popStatus(household,2)=popStatus(household,2)-1; %E -> E-1
        popStatus(household,3)=popStatus(household,3)+1; %I -> I+1
        
    else %event==3 Recovery
        popStatus(household,3)=popStatus(household,3)-1; %I -> I-1
        popStatus(household,4)=popStatus(household,4)+1; %R -> R+1
    end
    %sdisp(popStatus);
    timer=timer+toc;
end

%plot(cumsum(eventTime),totalInfected);
%profile viewer
end