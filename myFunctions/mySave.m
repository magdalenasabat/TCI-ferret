function mySave(str, var)
% wrapper on save function for parfor loops
save( str, 'var' , '-v7.3');