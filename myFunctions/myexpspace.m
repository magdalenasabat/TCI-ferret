function x = myexpspace(a,b,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates logarithmically spaced vector between a and b
% more intuitive than logspace, can start at 0
% 
% example 
% figure(); scatter(1:1:10,fliplr(exp(linspace(0,1,10))./cumsum(exp(linspace(0,1,10)))))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = fliplr(exp(linspace(a,b,n))./cumsum(exp(linspace(a,b,n))));


end