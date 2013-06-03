function rs=simH_constexp(beta,gamma,sigma,ro,tau,onDist,zeta,offDist,kappa,ita,N,alpha)
reps=200000;

rs=zeros(1,reps);
%Constant time on, Exponential off
for i=1:reps
    t=0;
    S=N-1;
    E=1;
    I=0;
    R=0;
    AntiviralArrivalTime=Inf;
    OnTheWay=0;
    AV=0;
    
    while t < AntiviralArrivalTime && (E+I)>0
        rates=zeros(1,3);
        %First run until AntiviralArrivalTime
        %There is no antivirals here.
        rates(1)=beta./(N-1) * S * I;
        rates(2)=sigma*E;
        rates(3)=gamma*I;
        
        %Determine next event time
        totalRate=sum(rates);
        eventTime=exprnd(1/totalRate);
        
        u=rand();
        %Choose event
        event=find(u<(cumsum(rates)./totalRate),1);
        
        %If the antivirals would arrive before the next event
        if t+eventTime > AntiviralArrivalTime
            %The next event actually happens at time
            %(AntiviralArriveTime-t)
            eventTime=AntiviralArrivalTime-t;
            %Add in the contribution to this event
            rs(i)=rs(i)+alpha*I*eventTime;
            %Set the antivirals to be on
            AV=1;
            %Update time
            t=t+eventTime;
            %Then stop this section
            break;
        end
        
        %Add in contribution to rStar
        rs(i)=rs(i)+alpha*I*eventTime;
        switch event
            case 1
                %Infection
                S=S-1;
                E=E+1;
                
            case 2
                %Progression
                E=E-1;
                I=I+1;
                %If antivirals aren't OnTheWay, set them to be.
                if OnTheWay==0
                    AntiviralArrivalTime=zeta;
                    OnTheWay=1;
                end
                
            case 3
                %Recovery
                I=I-1;
                R=R+1;
        end
        %Update current time
        t=t+eventTime;
        
    end
    %When we get to here, the Antivirals have started, or the epidemic
    %has ended.
    %AVs 'switch off' at rate 1/kappa, and have mean duration kappa.
    
    while (E+I)>0 && AV==1
        
        %Rates (Antivirals are present here)
        rates=zeros(1,4);
        rates(1)=beta./(N-1) * S * I * (1-tau) * (1-ro);
        rates(2)=sigma*E;
        rates(3)=gamma*I;
        rates(4)=1/kappa;
        
        %Determine next event time
        totalRate=sum(rates);
        eventTime=exprnd(1/totalRate);
        
        u=rand();
        %Choose event
        event=find(u<(cumsum(rates)./totalRate),1);
        
        %Add in contribution to rStar
        rs(i)=rs(i)+(1-tau)*alpha*I*eventTime;
        
        switch event
            case 1
                S=S-1;
                E=E+1;
                
            case 2
                E=E-1;
                I=I+1;
               
            case 3
                I=I-1;
                R=R+1;
                
            case 4
                %Here, antivirals have switched off
                AV=0;
        end
        t=t+eventTime;
    end
    
    while (E+I)>0
        %Now here, the antivirals have switch off again. This time, they
        %don't switch back on.
        rates=zeros(1,3);

        rates(1)=beta./(N-1) * S * I;
        rates(2)=sigma*E;
        rates(3)=gamma*I;
        
        %Determine next event time
        totalRate=sum(rates);
        eventTime=exprnd(1/totalRate);
        
        u=rand();
        %Choose event
        event=find(u<(cumsum(rates)./totalRate),1);
        
        %Add in contribution to rStar
        rs(i)=rs(i)+alpha*I*eventTime;
        switch event
            case 1
                S=S-1;
                E=E+1;
                
            case 2
                E=E-1;
                I=I+1;
                
            case 3
                I=I-1;
                R=R+1;
        end
        t=t+eventTime;
    end
end

rs=sum(rs)/reps;
        
        


end