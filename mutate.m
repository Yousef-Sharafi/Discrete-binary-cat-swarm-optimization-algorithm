function z=mutate(x,nvar,nb,srd,tb,bitt)
  CostFunction=@(x1,tb,bitt) cost_function(x1,tb,bitt);   
  costv=zeros(1,size(x,1));
  P=costv;
  costv(1)=CostFunction(x(1,:),tb,bitt);
   for i=2:size(x,1)
        j=randsample(nvar,nb)';            
        A=zeros(nvar);
        A(j)=1;
        r=find(A==1);
%         pos_neg=rand(1,numel(r));
%         for t=1:numel(r)
%             if(pos_neg(t)<=0.5)
%                 pos_neg(t)= 1;
%             else
%                 pos_neg(t)=-1;
%             end
%         end
      x(i,:)=x(1,:);
      for t=1:numel(r)
          if(rand<srd)
             x(i,r(t))=1-x(i,r(t));
          end
      end
        costv(i)=CostFunction(x(i,:),tb,bitt);
   end
   if(max(costv) == min(costv))
     P=ones(1,size(x,1));  
   else 
     P=abs(costv-max(costv))/(max(costv)-min(costv));
   end
    P=P./sum(P);
%    [~,uu]=max(P);
%     x(1,:)
%     costv(1)
    r1=RouletteWheelSelection(P);
    z=x(r1,:);
%     z=[x(r1,:) costv(r1)] ;
end
