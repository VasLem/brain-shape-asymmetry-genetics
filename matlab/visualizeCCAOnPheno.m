nP = length(pheno.clusterPCAPhenoFeatures);
clusterArray = load('../results/hierarchicalClusteringDemo/asymmetry_reduction10/ccPriorSegmentation/levels4_mine/segmentation.mat').clusterArray;
template = load('/home/vlemon0/home/code/results/hierarchicalClusteringDemo/asymmetry_reduction10/input.mat').preprocTemplate;
%%
[fig, handles] = paintClusters(clusterArray, template, 5, false);
%%
RESULTS_DIR = '../results/visualizeCCAOnPheno/';
%%
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
EPS_DIR = [RESULTS_DIR 'eps/'];
if ~isfolder(EPS_DIR), mkdir(EPS_DIR); end
PNG_DIR = [RESULTS_DIR 'png/'];
if ~isfolder(PNG_DIR), mkdir(PNG_DIR); end
GRAPHVIZ_DIR = [RESULTS_DIR 'graphviz/'];
if ~isfolder(GRAPHVIZ_DIR), mkdir(GRAPHVIZ_DIR); end

for i=1:length(handles)

f1 = figure('visible','off'); 
s = copyobj(handles{i}.handle,f1);
set(s(1),'position',[0.08 0.25 0.5 0.5]);
set(s(2),'position',[0.42 0.25 0.5 0.5]);
exportgraphics(f1,[EPS_DIR num2str(i) '.eps'], 'ContentType', 'vector',  'BackgroundColor', 'none');
system(['gs -dSAFER -dEPSCrop -r600 -sDEVICE=pngalpha -o ' PNG_DIR num2str(i) '.png ' EPS_DIR num2str(i) '.eps']);
close(f1);
end
%%
bCId ='Wout';

text = 'digraph G{';
leg1 = '';
leg2 = '';
leg3 = '';
cmap = jet(22);
for chr=1:22
    path = ['../results/genomeDemo/chr' num2str(chr) '_perSNP_PartitionedGTL' bCId 'BC_feats0significant_snps.csv'];
    if ~isfile(path)
        disp('Chromosome ' num2str(chr) ' snps file ' path ' not found. Skipping.');
        continue
    end
    snpsTable = readtable(['../results/genomeDemo/chr' num2str(chr) '_perSNP_PartitionedGTL' bCId 'BC_feats0significant_snps.csv']);
    chrColor = ['#' convertStringsToChars(join(convertCharsToStrings(dec2hex(round(255 * cmap(chr,:))))))];
    for i=1:length(handles)
    text = [text, '\n\t', num2str(i) '[shape=plaintext, image="' PNG_DIR num2str(i) '.png",label="",fixedsize=true,width=2,fontcolor=' chrColor ']' ];
    end
    for i=1:length(handles)
        parP = handles{i}.parent;
        counts = sum((snpsTable.PARTITION == i) & (snpsTable.CHR==chr));
        pCounts =  sum((snpsTable.PARTITION == parP) & (snpsTable.CHR==chr));

        if parP ~= 0
            text =  [text, '\n\t', num2str(parP) '->' num2str(i) '[style=invis]' ];
        end
        if counts
            if parP ~= 0
                text = [text, '\n\t', num2str(parP) '->' num2str(i) '[label="' num2str(counts-pCounts) '", color=' chrColor ', fontcolor=' chrColor ']' ];
            else
                leg1 = [leg1 '\n<tr><td align="right" port="i' num2str(chr) '" >Chr.' num2str(chr) '</td></tr>'];
                leg2 = [leg2 '\n<tr><td port="i' num2str(chr) '">&nbsp;</td></tr>'];
                leg3 = [leg3 '\nkey1:i' num2str(chr) ':e -> key2:i' num2str(chr) ':w [color=' chrColor ']'];
            end
            text = [text, '\n\t', num2str(i) '->' num2str(i) '[label="' num2str(counts) '", color=' chrColor ', fontcolor=' chrColor ',arrowhead=none]' ];
        end
    end
    break
end
legendText = ['digraph {\n'...
  'rankdir=LR\n'...
  'node [shape=plaintext]\n'...
  'subgraph cluster_01 {rank=same; key1, key2 \nlabel = "Legend";'...
    '\n key1 [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">['...
    leg1...
    '</table>>]'...
    '\n key2 [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">['...
    leg2...
    '</table>>]'...
    leg3 '\n}\n}'];
    
text = [text, '\n}\n'];
fileID = fopen([GRAPHVIZ_DIR 'graph.gv'],'w');
fprintf(fileID,text);
fclose(fileID);
fileID = fopen([GRAPHVIZ_DIR 'legGraph.gv'],'w');
fprintf(fileID,legendText);
fclose(fileID);

system(['twopi -Tpng -o ' RESULTS_DIR 'significanceSNPPropagation.png ' GRAPHVIZ_DIR 'graph.gv -Granksep=2 -Gratio=auto -Gfontsize=15'])
system(['dot -Tpng -o ' RESULTS_DIR 'legSignificanceSNPPropagation.png ' GRAPHVIZ_DIR 'legGraph.gv'])
g = imread([RESULTS_DIR 'significanceSNPPropagation.png']);
figure;
imshow(g);
title([bCId 'Bonferonni Correction']);
g = imread([RESULTS_DIR 'legSignificanceSNPPropagation.png ']);
figure;
imshow(g);
