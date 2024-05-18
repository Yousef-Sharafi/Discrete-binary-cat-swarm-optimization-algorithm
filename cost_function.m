function z=cost_function(x,tb,bitt)

up=50;
low=-50;
upb=mybin2dec(ones(1,bitt));
lowb=mybin2dec(zeros(1,bitt));
  xn=[];
  for j=1:tb
      x1=x(1+(j-1)*bitt:j*bitt);
      x1=mybin2dec(x1)/(upb-lowb)*(up-low)+low;
      xn=[xn x1];
  end 
% function for test fw
val=0;
n=size(xn,2);
for i=1:n
    val=val+(xn(1)-xn(i)^2)^2+(xn(i)-1)^2;
end
z=val;