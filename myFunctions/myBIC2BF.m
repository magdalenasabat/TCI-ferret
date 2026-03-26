function BF = myBIC2BF(BIC_0,BIC_1)
% BIC_0 is the tested model
% BIC 1 is the null model
% Transforms BIC to BF using : BF_{10} = e^{(BIC_0 - BIC_1)/2}
% magdalena sabat madalena.sabat@ens.psl.eu
% Wagenmakers, E. J. (2007). A practical solution to the pervasive problems of p values. Psychonomic bulletin & review, 14(5), 779-804 

BF = exp((BIC_0 - BIC_1)/2);

end

