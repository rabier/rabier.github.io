function [threshold]=supchi2ornsteindelong(L,df,p)


%this function calculate the quantile of the sup of an ornstein uhlenbeck chi square process (OUCS) with
%df degrees of freedom on [0,L]
%it uses delong results : it is approximate results

%the quantile is named threshold in this function: it is the result of this
%function
%it is inspired of formula on page 2205 of the article of Delong whose title 
% is "Crossing probabilities for a square root boundary by a bessel
% process"
%we have just adjusted this formula for the case of the ornstein Uhlenbeck
%chi square process


%in parameter
%L is the upper bound of the interval studied
%df is the number od degrees of freedom
% p is defined such as
%P(OUCS<threshold)=p


T=exp(4*L);
q=df/2;


%remark : the threshold obtained will be an approximation
%because the formula of delong p2205 is an approximation
%T has to be large enough
%c of the formula has to be large enough (it is the same c as 
%used further)

%in order that the formula be suitable, c has to be large enough
%here we take c=2 to begin
c=2;
E=10;
alpha=1-p;

erreur=zeros(10000,1);


for i=1:10000

c=c+0.001;
E= ( ((c^2/2)^q)* exp(-c^2/2)/gamma(q)  )*  ( log(T)*(1-2*q/c^2)+2/c^2);

erreur(i)=abs(E-alpha);

end;


%plot([1:10000],erreur)
u=find(erreur<0.001);

%be careful the function E is not monotone
%It crosses 0 not just one time


thresholdroot=u(size(u,1))/1000 +2;


%plot([1:10000],E);


threshold=thresholdroot^2;