function [threshold]=supchi2ornstein(nbtraj,step,L,df,p)


%this function calculates the quantile of the sup of an ornstein uhlenbeck chi square process (OUCS) with
%df degrees of freedom on [0,L]

%the quantile is named threshold in this function: it is the result of this
%function


%in parameter
%nbtraj is the number of trajectories simulated
%step is the step of discretization
%L is the upper bound of the interval studied
% p is defined such as
%P(OUCS<threshold)=p


sup=zeros(1,nbtraj);
suptrie=zeros(1,nbtraj);
ind=nbtraj*p;

for i=1:nbtraj

[Y]=chi2ornstein(step,L,df);

sup(i)=max(Y);

end

suptrie=sort(sup);

threshold=suptrie(ind);

