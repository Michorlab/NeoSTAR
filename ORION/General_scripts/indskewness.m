function diffs = indskewness(gmm,data)
%function to generate data from each
%gaussian in a GMM and compare skewnesses

%first calculate postprob
postprob = posterior(gmm,data);
%then, generate random vector to assign cells
randdraw = rand(size(data,1),1);
%then, for each gaussian, pull out data, calculate
%univariate skewness, (which should be 0 for a Gaussian)
for i = 1:gmm.NumComponents
    sampleid = randdraw<postprob(:,i);
    sample = data(sampleid,:);
    sampleskew = skewness(sample);
    diffs(i) = max(abs(sampleskew));
end
end