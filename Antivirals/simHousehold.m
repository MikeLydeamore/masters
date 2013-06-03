function rstar=simHousehold(beta,gamma,alpha,sigma,tauin,roin,onDist,zeta,offDist,kappa,N)

reps=200000;
rstar=zeros(1,reps);
switch onDist
	case 'exp'
		switch offDist
			case 'exp'
				for i=1:reps
					infectionOccured=0;
					S=N-1;
					E=1;
					I=0;
					R=0;
					AV=0;
					
					while(E+I > 0)
						%Calculate next rates
						%Infection
						rates(1)=beta/(N-1) * S * I * ((1-(AV==1)*tauin)*(1-(AV==1)*roin));
						%rates(1)=beta/(N-1) * S * I * (1-tau2)*(1-ro2);
						%Progression
						rates(2)=sigma*E;
						%Recovery
						rates(3)=gamma*I;
						%Antiviral Introduction/removal
						if AV==0
							rates(4)=infectionOccured*zeta;
						elseif AV==1
							rates(4)=1/kappa;
						else
							rates(4)=0;
						end
						
						%Determine next event time
						totalRate=sum(rates);
						%eventTime=exprnd(1/totalRate);
						eventTime=1/totalRate * log(1/rand());
						u=rand();
						
						%Choose event
						event=find(u<(cumsum(rates)./totalRate),1);
						
						modifier=(1-(AV==1)*tauin);
						
						rstar(i)=rstar(i)+alpha*modifier*I*eventTime;
						%Change state
						switch event
							case 1
								S=S-1;
								E=E+1;
								
							case 2
								E=E-1;
								I=I+1;
								infectionOccured=1;
								
							case 3
								I=I-1;
								R=R+1;
								
							case 4
								if AV==0
									AV=1;
								elseif AV==1
									AV=2;
								end
						end
						

					end
				end
				
			case 'const' %exp on, const off
				for i=1:reps
					infectionOccured=0;
					S=N-1;
					E=1;
					I=0;
					R=0;
					AV=0;
					timeLeftForAVs=Inf;
					tau=0;
					ro=0;
					
					while(E+I > 0)
						%Calculate next rates
						%Infection
						rates(1)=beta/(N-1) * S * I * (1-tau)*(1-ro);
						%rates(1)=beta/(N-1) * S * I * (1-tau2)*(1-ro2);
						%Progression
						rates(2)=sigma*E;
						%Recovery
						rates(3)=gamma*I;
						%Antiviral Introduction
						if AV==0
							rates(4)=infectionOccured*zeta;
						else
							rates(4)=0;
						end
						
						%Determine next event time
						totalRate=sum(rates);
						%eventTime=exprnd(1/totalRate);
						eventTime=1/totalRate * log(1/rand());
						u=rand();
						
						%Choose event
						event=find(u<(cumsum(rates)./totalRate),1);
						
						modifier=(1-tau);
						
						if AV==1
							if timeLeftForAVs-eventTime < 0
								eventTime=timeLeftForAVs;
								tau=0;
								ro=0;
								event=6;
								AV=2;
							else
								timeLeftForAVs=timeLeftForAVs-eventTime;
							end
						end
						
						rstar(i)=rstar(i)+alpha*modifier*I*eventTime;
						
						
						%Change state
						switch event
							case 1
								S=S-1;
								E=E+1;
								
							case 2
								E=E-1;
								I=I+1;
								infectionOccured=1;
								
							case 3
								I=I-1;
								R=R+1;
								
							case 4
								if AV==0
									AV=1;
									timeLeftForAVs=kappa;
									tau=tauin;
									ro=roin;
								elseif AV==1
									AV=2;
								end
							otherwise
								%Do nothing, AV's finished
						end
						

					end
				end
		end
	case 'const'
		switch offDist
			case 'exp' %const on, exp off
				for i=1:reps
					S=N-1;
					E=1;
					I=0;
					R=0;
					AV=0;
					timeLeftToAVs=Inf;
					tau=0;
					ro=0;
					infectionOccured=0;
					
					while(E+I > 0)
						%Calculate next rates
						%Infection
						rates(1)=beta/(N-1) * S * I * (1-tau)*(1-ro);
						%rates(1)=beta/(N-1) * S * I * (1-tau2)*(1-ro2);
						%Progression
						rates(2)=sigma*E;
						%Recovery
						rates(3)=gamma*I;
						%Antiviral Introduction
						if AV==1
							rates(4)=1/kappa;
						else
							rates(4)=0;
						end
						
						%Determine next event time
						totalRate=sum(rates);
						%eventTime=exprnd(1/totalRate);
						eventTime=1/totalRate * log(1/rand());
						u=rand();
						
						%Choose event
						event=find(u<(cumsum(rates)./totalRate),1);
						
						modifier=(1-tau);
						
						if AV==0
							if timeLeftToAVs-eventTime < 0
								eventTime=timeLeftToAVs;
								tau=tauin;
								ro=roin;
								event=6;
								AV=1;
							else
								timeLeftToAVs=timeLeftToAVs-eventTime;
							end
						end
						
						rstar(i)=rstar(i)+alpha*modifier*I*eventTime;
						
						%Change state
						switch event
							case 1
								S=S-1;
								E=E+1;
								
							case 2
								E=E-1;
								I=I+1;
								infectionOccured=infectionOccured+1;
								
								if infectionOccured==1
									timeLeftToAVs=zeta;
								end
								
							case 3
								I=I-1;
								R=R+1;
								
							case 4
								if AV==1
									AV=2;
									tau=0;
									ro=0;
								end
							otherwise
								%Do nothing, AV's introduced
						end
						

					end
				end
				
			case 'const'
				for i=1:reps
					S=N-1;
					E=1;
					I=0;
					R=0;
					AV=0;
					timeLeftToAVs=Inf;
					timeLeftForAVs=Inf;
					tau=0;
					ro=0;
					infectionOccured=0;
					
					while(E+I > 0)
						%Calculate next rates
						%Infection
						rates(1)=beta/(N-1) * S * I * (1-tau)*(1-ro);
						%rates(1)=beta/(N-1) * S * I * (1-tau2)*(1-ro2);
						%Progression
						rates(2)=sigma*E;
						%Recovery
						rates(3)=gamma*I;
						
						%Determine next event time
						totalRate=sum(rates);
						%eventTime=exprnd(1/totalRate);
						eventTime=1/totalRate * log(1/rand());
						u=rand();
						
						%Choose event
						event=find(u<(cumsum(rates)./totalRate),1);
						
						modifier=(1-tau);
						
						if AV==0
							if timeLeftToAVs-eventTime < 0
								eventTime=timeLeftToAVs;
								tau=tauin;
								ro=roin;
								event=6;
								AV=1;
								timeLeftForAVs=kappa;
							else
								timeLeftToAVs=timeLeftToAVs-eventTime;
							end
							
						elseif AV==1
							if timeLeftForAVs-eventTime < 0
								eventTime=timeLeftForAVs;
								tau=0;
								ro=0;
								event=6;
								AV=2;
							else
								timeLeftForAVs=timeLeftForAVs-eventTime;
							end
						end
						
						rstar(i)=rstar(i)+alpha*modifier*I*eventTime;
						
						%Change state
						switch event
							case 1
								S=S-1;
								E=E+1;
								
							case 2
								E=E-1;
								I=I+1;
								infectionOccured=infectionOccured+1;
								
								if infectionOccured==1
									timeLeftToAVs=zeta;
								end
								
							case 3
								I=I-1;
								R=R+1;
								
							otherwise
								%Do nothing, AV's introduced
						end
						

					end
				end
		end
end

rstar=sum(rstar)/reps;
				
end
