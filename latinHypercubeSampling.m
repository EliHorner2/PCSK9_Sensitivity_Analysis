function samples = latinHypercubeSampling(problem, N)
arguments
    problem;
    N;
end
%  latinHypercubeSampling generates a latin hypercube sample.
%  samples=latinHypercubeSampling(problem,N) generates 
% a latin hypercube samples containing N values on each of parameters
% described in the problem and scales them by upper and lower bounds.
% Inputs
%   problem - struct
        % problem has to have the following fields:
        %    lowerBound - desired lower bound of each of the parameters
        %    upperBound - desired upper bound of each of the parameters
        %    num_vars - number of parameters 
        %    names - optionally name of the parameters 
%   N - int
% Code by Jaimit Parikh

if ~all(isfield(problem, ["lowerBound", "upperBound", "num_vars"]))
    error("The necessary fields in the problem structure are not included."...
        +" The required field names are lowerBound, upperBound, and num_vars")
end

iniSamples = lhsdesign(N, problem.num_vars);
samples = rescale(iniSamples,...
        problem.lowerBound, ...
        problem.upperBound,...
        'InputMin', min(iniSamples),...
        'InputMax', max(iniSamples));
end