beta=1;
gamma=10;
alpha=1;
eps=2;

maxIter=3*houseSize*N;

N=10; %Number of households
houseSize=4;

popStatus=zeros(N,3); %S-I-R
popStatus(2:end,1)=houseSize*ones(N-1,1);
popStatus(1,1)=houseSize-1; %1 infected
popStatus(1,2)=1;

max_events=3*houseSize*N;

eventTime=zeros(1,max_events);

%Rate vector will go [Internal Internal ... Internal External External... External ...] etc
rates=zeros(N,3); %One row per household

totalInfected=zeros(1,maxIter);
for i=1:maxIter
    %% Rates
    totalInfected(i)=sum(popStatus(:,2));
    %Internal infection at rate BIS
    rates(:,1)=beta*popStatus(:,2).*popStatus(:,1);
    
    %External infection at rate eps*totalI*S (not sure what eps is yet)
    rates(:,2)=eps*totalInfected(i).*popStatus(:,1);
    
    %Recovery at rate gamma*I
    rates(:,3)=gamma*popStatus(:,2);
    
    %% Reshape, normalise and choose which event occurs
    
    [m, n]=size(rates);
    rateVector=reshape(rates,1,m*n);
    
    normRateVector=rateVector./sum(rateVector);
    
    %Choose event
    eventTime(i)=rand();
    eventOccured=find(eventTime(i)<cumsum(normRateVector),1);
    
    %Determine event and household
    
    [household, event]=ind2sub([N, 3],eventOccured); %3 being the number of phases
    
    %% Make corresponding change
    
    if event==1 %Internal infection
        popStatus(household,1)=popStatus(household,1)-1; %S -> S-1
        popStatus(household,2)=popStatus(household,2)+1; %I -> I+1
        
    elseif event==2 %External infection
        popStatus(household,1)=popStatus(household,1)-1;
        popStatus(household,2)=popStatus(household,2)+1;
        
    else %event==3 Recovery
        popStatus(household,2)=popStatus(household,2)-1; %I -> I-1
        popStatus(household,3)=popStatus(household,3)+1; %R -> R+1
    end
    %sdisp(popStatus);
end

plot(1:maxIter,totalInfected);