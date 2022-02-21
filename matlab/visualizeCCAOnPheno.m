addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));
RECOMPUTE_PARTS = false;
DATASET_INDEX = 1;
switch DATASET_INDEX
    case 1
        UKBIOBANK = 'UKBIOBANK';
        DATASET = 'STAGE00DATA';
        GENO_ID = 'sel19908';
    case 2
        UKBIOBANK = 'MY_UKBIOBANK';
        DATASET = 'BATCH2_2021_DATA';
        GENO_ID = 'sel16875_rmrel';
end
GENO_DIR = ['../results/genomeDemo/' DATASET '/'];
CLUSTER_DIR = ['../results/hierarchicalClusteringDemo/' DATASET '/'];
RESULTS_DIR = fullfile(pwd, ['../results/visualizeCCAOnPheno/' DATASET '/']);
clusterArray = load([CLUSTER_DIR  'asymmetry_reduction10/ccPriorSegmentation/levels4_mine/segmentation.mat']).clusterArray;
template = load([CLUSTER_DIR 'asymmetry_reduction10/ccPriorSegmentation/levels4_mine/input_info.mat']).preprocTemplate;
%%
if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
[fig, fig2, handles] = paintClusters(clusterArray, template, 4, false);

%%

PNG_DIR_1 = [RESULTS_DIR 'png_raw/'];
if ~isfolder(PNG_DIR_1), mkdir(PNG_DIR_1); end
PNG_DIR_2 = [RESULTS_DIR 'png/'];
if ~isfolder(PNG_DIR_2), mkdir(PNG_DIR_2); end
GRAPHVIZ_DIR = [RESULTS_DIR 'graphviz/'];
if ~isfolder(GRAPHVIZ_DIR), mkdir(GRAPHVIZ_DIR); end
%% Uncomment to redraw partitions
for i=1:length(handles)
    ret_path = [PNG_DIR_1 num2str(i) '.png'];
    if ~RECOMPUTE_PARTS && isfile(ret_path)
        continue
    end
    f1 = figure('visible','off');
    s = copyobj(handles{i}.handle,f1);
    set(s(1),'position',[0.0 0.25 0.5 0.5]);
    set(s(2),'position',[0.5 0.25 0.5 0.5]);
    set(f1, 'color', 'none');
    set(s, 'color', 'none');
    print(f1, '-dpng', '-r300', ret_path);
    close(f1);
end
%%

bCIds = ["With", "Wout"];

cmap = brewermap(22, 'Spectral');
countsMat = zeros(2, length(handles), 22);
pCountsMat = zeros(2, length(handles), 22);
for j=1:length(handles)
    parP = handles{j}.parent;

    for i=1:2
        bCId = char(bCIds(i));
        for chr=1:22
            path = [ GENO_DIR 'chr' num2str(chr) '/PartitionedGTL' bCId 'BC_feats0significant_snps.csv'];
            if ~isfile(path)
%                 disp(['Chromosome ' num2str(chr) ' snps file ' path ' not found. Skipping.']);
                continue
            end
            snpsTable = readtable(path);
            if isempty(snpsTable)
                continue
            end
            
            countsMat(i, j, chr) = sum((snpsTable.PARTITION == j) & (snpsTable.CHR==chr));
            pCountsMat(i, j, chr) =  sum((snpsTable.PARTITION == parP) & (snpsTable.CHR==chr));
        end
    end
end
%%
for j=1:length(handles)
    p = [fullfile(pwd, [PNG_DIR_1 num2str(j)]) '.png'];
    im = [];
    for i =1:2
        bCId = char(bCIds(i));
        ret_path = [PNG_DIR_2 bCId '/' num2str(j) '.png'];
        if ~RECOMPUTE_PARTS && isfile(ret_path)
            continue
        end
        if isempty(im)
            im = imread(p);
            im = im(ceil(0.15*size(im,1)):size(im,1),:,:);
        end
        f = figure('Visible','off');
        imshow(im);
        annotation(f, "textbox", 'Position', [0.1, 0.75, 0.8, 0.15], 'String', num2str(j), ...
                FontUnits='Normalized', FontSize=0.1, ...
                    FontWeight='bold', ...
                    HorizontalAlignment='center', VerticalAlignment='middle', LineStyle='none');
        msk = countsMat(i, j, :) >0;
        if any(msk)
            inds = find(msk);
            cols = cmap(msk,:);
            annot_pos = linspace(0.1,0.9,sum(msk)+1);
            for s=1:length(annot_pos) - 1
                annotation(f,'textbox', 'Position', [annot_pos(s) 0.1 annot_pos(s+1) - annot_pos(s) .15], ...
                    'String', num2str(countsMat(i,j, inds(s))), 'BackgroundColor', cols(s,:), ...
                    FontUnits='Normalized',FontSize=0.06, ...
                    FontWeight='bold', ...
                    HorizontalAlignment='center', VerticalAlignment='middle');
            end
        end
        
        if ~isfolder([PNG_DIR_2 bCId]), mkdir([PNG_DIR_2 bCId]); end
        saveas(f, ret_path)

    close(f)
    end
end
%%
old_path = pwd;
cd(RESULTS_DIR)
for i=1:2
    bCId = char(bCIds(i));
    text = 'digraph G{';
    leg1 = '';
    leg2 = '';
    leg3 = '';
    for j=1:length(handles)
            parP = handles{j}.parent;
            if parP ~= 0
                text =  [text, '\n\t', num2str(parP) '->' num2str(j) '[ splines="false",style="invis"]' ];
            end
    end

    for chr=1:22
        chrColor = ['#' sprintf('%02X',round(255 * cmap(chr,:)))];
        for j=1:length(handles)
            p = ['png/' bCId '/' num2str(j)]; % path relative to the svg file location
            text = [text, '\n\t', num2str(j) '[shape=plaintext, image="' ...
                p '.png",label="",fixedsize=true,width=4,height=2,fontcolor="' chrColor '"]' ];
        end
        for j=1:length(handles)
            parP = handles{j}.parent;
            counts = countsMat(i, j, chr);
            pCounts =  pCountsMat(i, j, chr);
            if counts
                if parP ~= 0
%                     text = [text, '\n\t', num2str(parP) '->' num2str(j) '[label="' num2str(counts-pCounts) '", color="' chrColor '",fontsize=15, fontcolor="' chrColor '", splines="true"]' ];
                else
                    leg1 = [leg1 '\n<tr><td align="right" port="i' num2str(chr) '" >Chr.' num2str(chr) '</td></tr>'];
                    leg2 = [leg2 '\n<tr><td port="i' num2str(chr) '">&nbsp;</td></tr>'];
                    leg3 = [leg3 '\nkey1:i' num2str(chr) ':e -> key2:i' num2str(chr) ':w [color="' chrColor '"]'];
                end
            end
        end
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
    fileID = fopen([GRAPHVIZ_DIR  'graph' bCId 'BC.gv'],'w');
    fprintf(fileID,text);
    fclose(fileID);
    fileID = fopen([GRAPHVIZ_DIR  'legGraph' bCId 'BC.gv'],'w');
    fprintf(fileID,legendText);
    fclose(fileID);
    graph_png_path = [RESULTS_DIR 'significanceSNPPropagation' bCId 'BC.svg'];
    leg_png_path = [RESULTS_DIR 'legSignificanceSNPPropagation' bCId 'BC.svg'];
    system(['twopi -Tsvg -o ' graph_png_path ' ' GRAPHVIZ_DIR 'graph' bCId 'BC.gv -Granksep=2 -Gratio=0.5 -Gfontsize=20'])
    system(['dot -Tsvg -o ' leg_png_path ' ' GRAPHVIZ_DIR 'legGraph' bCId 'BC.gv'])
end
cd(old_path)