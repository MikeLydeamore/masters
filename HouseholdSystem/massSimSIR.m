function [avgInfected, t, simsUsed]=massSimSIR(beta,gamma,alpha,N,houseSize)
%Simulate SIR Model and plot average

% alpha=1;
% beta=6;
% eps=1;
% gamma=1;
% N=10;
% houseSize=4;
simsUsed=0;
iterations=20;
t=0:0.01:5;
avgInfected=zeros(1,length(t));
for i=1:iterations
    clc;fprintf('Iteration: %d \nSims Used: %d \n',i,simsUsed);
    [totalInfected, eventTime]=simulateSIR(beta,gamma,alpha,N,houseSize);
    %First, check if we have an outbreak.
    %If there has been at least 4 events, then we'll assume there's been a
    %usable outbreak.
    if sum(eventTime>0)>4
        %Check if we've gone past bounds of t
        if sum(eventTime)>max(t)
            avgInfected=[avgInfected zeros(1,length(t(end)+0.01:0.01:sum(eventTime)))];
            t=[t t(end)+0.01:0.01:sum(eventTime)];
        end
        
        cumTimeEvents=cumsum(eventTime);
        for j=1:length(t)
            if t(j) > cumTimeEvents(end) %If already absorbed at this time
                increase=length(cumTimeEvents);
            else
                increase=find(t(j)<=cumTimeEvents,1);
            end
            
            avgInfected(j)=avgInfected(j)+totalInfected(increase);
        end
        simsUsed=simsUsed+1;
    end
end

plot(t,avgInfected/simsUsed);

end