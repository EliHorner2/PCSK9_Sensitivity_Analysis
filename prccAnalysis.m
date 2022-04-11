function S = prccAnalysis(X, Y, options)
arguments
    X; %  n x k ... n input samples with k parameters
    Y; % n x 1 ... nsamples with 1 output variable of interest ...  Currently
    options = "partialcorr";
end
% PRCC  to estimate the sensitivity of
% the parameters X to the response variable Y
% Inputs:
%   X - nxk matrix .. n input samples with k parameters
%   Y - n x1 matrix .. n samples with 1 output variable of interest
%   options - string for the kind of  correlation analysis to perform. default =
%   "partialcorr", other option include "corrcoef";
% Ouputs:
 %  S - output strucutre of sensitivity analysis
 %      S.S1 - normalized values of the sensivitiy indices
% Code by Jaimit Parikh based on Gallaher et al. 2018 
% https://www.sciencedirect.com/science/article/abs/pii/S0022519318304259

X = normalize(X);
Y = normalize(Y);

switch options
    case "partialcorr"
        [prcc, pval] = partialcorr([X, Y], 'Type', 'Spearman');
        s = prcc(1:end-1, end);
        S.S1 =  s' / sum(s);
        S.prcc = prcc;
        S.pval = pval;
    case "corrcoef"
        prcc = corrcoef([X, Y]);
        s = prcc(1:end-1, end);
        S.S1 =  s' / sum(s);
end

end