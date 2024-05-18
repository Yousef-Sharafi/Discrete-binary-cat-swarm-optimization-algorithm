clc;
close all;
num=2;     
MaxIt=200;  % Maximum Number of Iterations
nPop=50; 
%% Algorithm Parameters BINARY CAT 2013    
tb=10;
bitt=20;
nVar=bitt*tb; 
BestCost1_cat=zeros(num,MaxIt);  
CostFunction=@(x,tb,bitt) cost_function(x,tb,bitt);   % Cost Function
c2_cat=1;
for ittt=1:num
    for ta=1:1 
            %  Number of Decision Variables
            alpha=0.3;
            VarSize=[1 nVar];   % Decision Variables Matrix Size
            %% PSO Parameters
            SMP=3;%0.25*nPop;
            SRD=0.2;
            CDC=0.2;
            nb=round(nVar*CDC);        
            MR=0.3;
            num_seek=round(MR*nPop);
            num_track=nPop-num_seek;
            cat=randperm(nPop);
            w_cat=0.5;
            vmax_cat=4;
            %********************************
            %% Initialization
            % Define Empty Structure to Hold Particle Data
            empty_cat.Position=[];
            empty_cat.flag=[];
            empty_cat.Velocity=[];
            empty_cat.Cost=[];
            pop=repmat(empty_cat,nPop,1);
            vel_cat=rand(nPop,nVar)-0.5;
            one_vel_cat=rand(nPop,nVar)-0.5;
            zero_vel_cat=rand(nPop,nVar)-0.5;
            % Initialize Global Best
            GlobalBest.Cost=inf;
            for i=1:nPop
                % Initialize Velocity
                pop(i).Position = round(rand(1,nVar));
                pop(i).Velocity = rand(1,nVar);
                % Evaluate Solution
                pop(i).Cost=CostFunction(pop(i).Position,tb,bitt); 
                y=find(cat==i);
                if(y<=num_seek)
                    pop(i).flag=1;
                else
                    pop(i).flag=0;
                end
                % Update Global Best
                if pop(i).Cost<=GlobalBest.Cost
                    GlobalBest=pop(i);
                end
            end
            % Define Array to Hold Best Cost Values
            BestCost=zeros(MaxIt,1);
            c1=1;
            %% PSO Main Loop
            for it=1:MaxIt
                    oneadd_cat=zeros(nPop,nVar);
                    zeroadd_cat=zeros(nPop,nVar);
                    dd3_cat=c2_cat*rand;
                    %******************************************************
                    for t_i=1:nPop
                        for g_i=1:nVar
                            if(GlobalBest.Position(g_i)==0)
                               oneadd_cat(t_i,g_i)=oneadd_cat(t_i,g_i)-dd3_cat;
                               zeroadd_cat(t_i,g_i)=zeroadd_cat(t_i,g_i)+dd3_cat;
                            else
                               oneadd_cat(t_i,g_i)=oneadd_cat(t_i,g_i)+dd3_cat; 
                               zeroadd_cat(t_i,g_i)=zeroadd_cat(t_i,g_i)-dd3_cat;
                            end
                        end
                    end
                    one_vel_cat=w_cat*one_vel_cat+oneadd_cat;
                    zero_vel_cat=w_cat*zero_vel_cat+zeroadd_cat;
                    
                    for t_i=1:nPop
                        for g_i=1:nVar
                            if(abs(vel_cat(t_i,g_i))>vmax_cat)
                                one_vel_cat(t_i,g_i)=vmax_cat*sign(one_vel_cat(t_i,g_i));
                                zero_vel_cat(t_i,g_i)=vmax_cat*sign(zero_vel_cat(t_i,g_i));
                            end
                        end
                    end
                    
                    for t_i=1:nPop
                        for g_i=1:nVar
                            if(pop(t_i).Position(g_i)==1)
                                vel_cat(t_i,g_i)=zero_vel_cat(t_i,g_i);
                            else
                                vel_cat(t_i,g_i)=one_vel_cat(t_i,g_i);
                            end
                        end
                    end
                    veln_cat=logsig(vel_cat);
                    %******************************************************
                for i=1:nPop
                    if(pop(i).flag==0)                    
                            for r2=1:nVar
                              if(rand<veln_cat(i,r2))
                                  pop(i).Position(r2)=GlobalBest.Position(r2);
                              else
                                  pop(i).Position(r2)=pop(i).Position(r2);
                              end
                            end
                            pop(i).Cost = CostFunction(pop(i).Position,tb,bitt); 
                    else         
                            copy_cat=repmat(pop(i).Position,SMP,1); 
                            pop(i).Position=mutate(copy_cat,nVar,nb,SRD,tb,bitt);                    
                    end
                    try
                    pop(i).Cost = CostFunction(pop(i).Position,tb,bitt); 
                    catch tt
                        disp('ll');
                    end
                        % Update Global Best
                        if pop(i).Cost<=GlobalBest.Cost
                            GlobalBest=pop(i);
                        end 
                end
                % Store Best Cost Ever Found
                BestCost(it)=GlobalBest.Cost;  
                % Show Iteration Information
                disp(['Iteration ' num2str(it) ': Best Cost cat= ' num2str(BestCost(it))]);
                num_seek=round(MR*nPop);
                num_track=nPop-num_seek;
                cat=randperm(nPop);
                for ii=1:nPop
                    y=find(cat==ii);
                    if(y<=num_seek)
                        pop(ii).flag=1;
                    else
                        pop(ii).flag=0;
                    end
                end
            end
    end % BINARY CAT 2013
    BestCost1_cat(ittt,:)=BestCost';
end
    BestCost_all_cat1=mean(BestCost1_cat);
    BestCost1_cat=BestCost1_cat;
    %******************************************
    std1_cat1=std(BestCost1_cat);
    std_cat1=std1_cat1(MaxIt);
    %****************************************** 
    mean1_cat1=mean(BestCost1_cat);
    mean_cat1=mean1_cat1(MaxIt);
    %******************************************
    best1_cat1=max(BestCost1_cat);
    best_cat1=best1_cat1(MaxIt);
    %******************************************
    bad1_cat1=min(BestCost1_cat);
    bad_cat1=bad1_cat1(MaxIt);
    daa=[mean_cat1;std_cat1;bad_cat1;best_cat1]
  figure;
  semilogy(BestCost_all_cat1,'r','LineWidth',2);
  xlabel('Iteration');
  ylabel('Cost');
  legend('binary cat','Location','NW');



