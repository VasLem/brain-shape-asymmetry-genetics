function [ data, titlenames, labels, thresholds, totalStats] = ProcrustesAnova2WayAsymmetryOutputProcess(brainSurface, showstruct, nSamplesPerPick, showPerm, saveMatFileName)
if nargin < 2
    showPerm = true;
end
if nargin < 3
    saveMatFileName = '';
end
i=1;
VertexValues{i,1} = showstruct.LM.I;titlenames{i,1} = 'I';
VertexValues{i,2} = showstruct.LM.IF;titlenames{i,2} = 'IF';
if showPerm
    VertexValues{i,3}=showstruct.LM.permIF;
else
    VertexValues{i,3}=showstruct.LM.IP;
end
titlenames{i,3} = 'p-Value';
i=i+1;
VertexValues{i,1} = showstruct.LM.D;titlenames{i,1} = 'D';
VertexValues{i,2} = showstruct.LM.DF;titlenames{i,2} = 'DF';
if showPerm
    VertexValues{i,3}=showstruct.LM.permDF;
else
    VertexValues{i,3}=showstruct.LM.DP;
end
titlenames{i,3} = 'p-Value';

i=i+1;
VertexValues{i,1} = showstruct.LM.F;titlenames{i,1} = 'F';
VertexValues{i,2} = showstruct.LM.FF;titlenames{i,2} = 'FF';
if showPerm
    VertexValues{i,3}=showstruct.LM.permFF;
else
    VertexValues{i,3}=showstruct.LM.FP;
end
titlenames{i,3} = 'p-Value';

labels = ["Individual", "Directional", "Fluctuating"];
thresholds = [0.05, 0.01, 0.005, 0.0001];
totalStats = showstruct.Total;
nValues = 3;
for i=1:nValues
    val = VertexValues{i,3};
    res = zeros(size(val));
    c = 1;
    for t=thresholds
        res(val<=t) = c;
        c = c+1;
    end
    VertexValues{i,4} = res;
    titlenames{i,4} = "Significance";
end
data = VertexValues;
if ~isempty(saveMatFileName)
save(saveMatFileName, "data","nSamplesPerPick", "titlenames","labels","thresholds","totalStats","brainSurface",'-v7.3');
end
end