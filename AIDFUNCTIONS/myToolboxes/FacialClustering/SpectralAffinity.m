function label = SpectralAffinity(A)

%A = buildRVmatrix(Module);

if length(A)==1, 
    label=1; 
    return % condition put in case of only 1 LM  m
else
    [label,~,~] = spectral_clustering_weiss_ng(A, 2);
end

end

