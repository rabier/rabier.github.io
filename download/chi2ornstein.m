function [Y]=chi2ornstein(step,L,df)

%this function simulates one trajectory of an Ornstein Uhlenbeck chi square process
%with df degrees of freedom on [0,L] 
%Y is this trajectory 


%parameters
%step is the step of discretization
%L is the upper bound of the interval studied

vect=[0:step:L];

N=size(vect,2);

T=zeros(N,df);


for i=1:df

T(1:N,i)=simornstein(step,L);

end


Tcarre=T.^2;


Y=sum(Tcarre,2);

%plot(vect,Y)
