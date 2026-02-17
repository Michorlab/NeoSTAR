function G = RipleyG_adj(speciesA, speciesB, r, area)
    % Compute density of species B
    %area = width * height;
    
    nB = size(speciesB, 1);
    lambdaB = nB / area;  % intensity of species B (points per unit area)
    
    % Expected distance to nearest neighbor under CSR for species B
    % CSR expected G(r): G_CSR(r) = 1 - exp(-pi * lambdaB * r^2)
    
    nA = size(speciesA, 1);
    minDist = zeros(nA, 1);

    for i = 1:nA
        dists = sqrt((speciesB(:,1) - speciesA(i,1)).^2 + (speciesB(:,2) - speciesA(i,2)).^2);
        minDist(i) = min(dists);
    end

    % Empirical G(r): fraction of speciesA points with nearest neighbor within r
    G_empirical = arrayfun(@(radius) mean(minDist <= radius), r);

    % CSR baseline G-function
    G_CSR = 1 - exp(-pi * lambdaB * r.^2);

    % Normalize by CSR expectation
    G = G_empirical ./ (G_CSR + eps);  % eps to avoid division by 0
end