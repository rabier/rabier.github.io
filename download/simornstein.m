function [X]=simornstein(step,L)

%this function simulates one trajectory of an
% ornstein Uhlenbeck process on [0,L]
%the result is the trajectory X



%parameter:
%step is the step of discretization
%L is the upper bound of the interval studied


vect=[0:step:L];

N=size(vect,2);

X=zeros(N,1);


X(1)=randn;


for i=1:N-1

epsilon= sqrt(1-exp(-4*step))*randn;
X(i+1)=exp(-2*step)*X(i)+ epsilon;

end




%plot(vect,X)
