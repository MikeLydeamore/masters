function [totalInfected, eventTime]=simulateSIR(beta,gamma,alpha,N,houseSize)
%houseSize is a distribution of house sizes, ie [0.5 0.5] implies that N/2
%houses are of size 1 and N/2 sizes are of size 2, whereas [0.3 0.3 0.4]
%implies that 0.3*N houses are of size 1, 0.3*N are of size 2 and 0.4*N are
%of size 3.

%Determine how many houses are of each type (must add up to N).
housePop=zeros(1,length(houseSize));
for i=1:length(houseSize)-1
    housePop(i)=round(N*houseSize(i));
end
housePop(end)=N-sum(housePop);

%Population size
popSize=housePop*(1:length(housePop))';

maxIter=3*popSize;

%Assign each household it's population
popStatus=zeros(N,3); %S-I-R

assigningPop=housePop;
currentRow=1;
for i=1:length(housePop)
    %Want to assign housePop(i) rows of the popStatus matrix to households
    %of size i
    while assigningPop(i)>0
        popStatus(currentRow,1) = i ; %Loads house with susceptibles
        assigningPop(i)=assigningPop(i)-1;
        currentRow=currentRow+1;
    end
end
    
%Infect some people    
% popStatus(1,1)=popStatus(1,1)-1;
% popStatus(1,2)=popStatus(1,2)+1; %Infect one
popStatus(1:10,1)=popStatus(1:10,1)-1;
popStatus(1:10,2)=popStatus(1:10,2)+1; %Infect ten.

eventTime=zeros(1,maxIter);

%Rate vector will go [Infection Infection ... Recovery Recovery ...] etc
rates=zeros(N,2); %One row per household



totalInfected=zeros(1,maxIter);
for i=1:maxIter
    %% Rates
    totalInfected(i)=sum(popStatus(:,2));
    if totalInfected(i)==0 %No infected left, in absorbing state
        break;
    end
    %Infection (both internal and external) at rate B_k*I*S+alpha*totalI/N*S
    beta_k=zeros(N,1);
    householdPop=sum(popStatus,2);
    oneOnPop=zeros(length(householdPop),1);
    oneOnPop(find(householdPop>1))=1./(householdPop(find(householdPop>1))-1);
    beta_k=beta.*oneOnPop; %beta_k = beta/(k-1) only for k>1, otherwise have = 0
    rates(:,1)=beta_k.*popStatus(:,2).*popStatus(:,1)+alpha.*(totalInfected(i)/popSize).*popStatus(:,1);
    
    %External infection at rate eps*totalI*S (not sure what eps is yet)
    %rates(:,2)=eps.*popStatus(:,1);
    
    %Recovery at rate gamma*I
    rates(:,2)=gamma*popStatus(:,2);
    
    %% Reshape, normalise and choose which event occurs
    
    [m, n]=size(rates);
    rateVector=reshape(rates,1,m*n);
    
    normRateVector=rateVector./sum(rateVector);
    
    %Choose event
    u=rand();
    eventOccured=find(u<cumsum(normRateVector),1);
    
    %Calculate event time
    if sum(rateVector)==0
        x=3;
    end
    eventTime(i+1)=exprnd(1/sum(rateVector));
    if isnan(eventTime(i+1));
        x=3;
    end
    
    %Determine event and household
    [household, event]=ind2sub([N, 3],eventOccured); %3 being the number of phases
    
    %% Make corresponding change
    
    if event==1 %Internal infection
        popStatus(household,1)=popStatus(household,1)-1; %S -> S-1
        popStatus(household,2)=popStatus(household,2)+1; %I -> I+1
        
%     elseif event==2 %External infection
%         popStatus(household,1)=popStatus(household,1)-1;
%         popStatus(household,2)=popStatus(household,2)+1;
        
    else %event==2 Recovery
        popStatus(household,2)=popStatus(household,2)-1; %I -> I-1
        popStatus(household,3)=popStatus(household,3)+1; %R -> R+1
    end
    %sdisp(popStatus);
end

%plot(cumsum(eventTime),totalInfected);

end