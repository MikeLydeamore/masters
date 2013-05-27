function [totalInfected, eventTime]=simulateSIR(beta,gamma,alpha,N,houseSize)

maxIter=3*houseSize*N;

popStatus=zeros(N,3); %S-I-R
popStatus(2:end,1)=houseSize*ones(N-1,1);
popStatus(1,1)=houseSize-1; %1 infected
popStatus(1,2)=1;

eventTime=zeros(1,maxIter);

%Rate vector will go [Infection Infection ... Recovery Recovery ...] etc
rates=zeros(N,2); %One row per household

popSize=N*houseSize;

totalInfected=zeros(1,maxIter);
for i=1:maxIter
    %% Rates
    totalInfected(i)=sum(popStatus(:,2));
    if totalInfected(i)==0 %No infected left, in absorbing state
        break;
    end
    %Infection (both internal and external) at rate B*I*S+alpha*totalI/N*S
    rates(:,1)=beta*popStatus(:,2).*popStatus(:,1)+alpha.*(totalInfected(i)/popSize).*popStatus(:,1);
    
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
    eventTime(i+1)=exprnd(1/sum(rateVector));
    
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

plot(cumsum(eventTime),totalInfected);

end