function [y] = myResample(x,originalFs,desiredFs)
%MYRESAMPLE simplified resample function

[p,q]  = rat(desiredFs/originalFs);

y = resample(x,p,q);

end

