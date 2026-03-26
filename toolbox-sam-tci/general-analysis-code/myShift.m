function out = myShift(X,shifts)


% by magdalena
%wrapperaround circshift

%pad
x = [zeros(size(X)); X(:); zeros(size(X))];

%shift
o = circshift(x,shifts);

%upad
out = o(length(X):length(X)+length(X));





end