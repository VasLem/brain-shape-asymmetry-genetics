function data= ProcrustesAnova2WayAsymmetryOutputProcess(brainSurface, showstruct, nSamplesPerPick, showPerm, saveMatFileName, smallestPvalue)
if nargin < 4
    showPerm = true;
end
if nargin < 5
    saveMatFileName = '';
end
if nargin < 6
    smallestPvalue = 0.0001;
end
i=1;
VertexValues{i,1} = showstruct.LM.I;titlenames{i,1} = 'I';
try
VertexValues{i,2} = showstruct.LM.IF;titlenames{i,2} = 'IF';
if showPerm
    VertexValues{i,3}=showstruct.LM.permIF;
else
    VertexValues{i,3}=showstruct.LM.IP;
end
titlenames{i,3} = 'p-Value';
catch
end
i=i+1;
VertexValues{i,1} = showstruct.LM.D;titlenames{i,1} = 'D';
try
VertexValues{i,2} = showstruct.LM.DF;titlenames{i,2} = 'DF';
if showPerm
    VertexValues{i,3}=showstruct.LM.permDF;
else
    VertexValues{i,3}=showstruct.LM.DP;
end
titlenames{i,3} = 'p-Value';
catch
end
i=i+1;
VertexValues{i,1} = showstruct.LM.F;titlenames{i,1} = 'F';
try
VertexValues{i,2} = showstruct.LM.FF;titlenames{i,2} = 'FF';
if showPerm
    VertexValues{i,3}=showstruct.LM.permFF;
else
    VertexValues{i,3}=showstruct.LM.FP;
end
titlenames{i,3} = 'p-Value';

catch
end
labels = ["Individual", "Directional", "Fluctuating"];
thresholds = 10.^-(log10(1/smallestPvalue)-(2:-1:0));
thresholds = [thresholds * 5; thresholds];
thresholds = thresholds(:)';

            
try
for i=1:3
    val = VertexValues{i,3};
    res = zeros(size(val));
    c = 1;
    for t=thresholds(1:end-1)
        res(val<=t) = c;
        c = c+1;
    end
    VertexValues{i,4} = res;
    titlenames{i,4} = "Significance";
end
catch
end
data.values = VertexValues;
data.titleNames = titlenames;
data.brainSurface = brainSurface;
try
data.total = showstruct.Total; % no total for AMMI
data.nSamplesPerPick = nSamplesPerPick;
data.labels = labels;
data.thresholds = thresholds;
catch
end
if ~isempty(saveMatFileName)
save(saveMatFileName, "data",'-v7.3');
end
end