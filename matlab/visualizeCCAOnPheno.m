addpath(genpath('.'));
addpath(genpath('AIDFUNCTIONS'));
RECOMPUTE_PARTS = true;
DATASET_INDEX = 1;
switch DATASET_INDEX
    case 1
        DATASET = 'STAGE00DATA';
    case 2
        DATASET = 'BATCH2_2021_DATA';
end
GENO_DIR = ['../results/genomeDemo/' DATASET '/'];
CLUSTER_DIR = ['../results/hierarchicalClusteringDemo/' DATASET '/'];
RESULTS_DIR = fullfile(pwd, ['../results/visualizeCCAOnPheno/' DATASET '/']);
clusterArray = load(['../results/hierarchicalClusteringDemo/STAGE00DATA/asymmetry_reduction10/levels4/segmentation.mat']).clusterArray;
template = load([CLUSTER_DIR 'asymmetry_reduction10/levels4/input_info.mat']).preprocTemplate;
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
mask = squeeze(any(countsMat, 2));
chromosomes = find(any(mask, 1));
cmap = brewermap(length(chromosomes), 'Spectral');
%%
for i=1:2
    bCId = char(bCIds(i));
    tmask_chrs = find(mask(i,:));
    [~, j] = find(chromosomes==tmask_chrs');
    cols = cmap(j, :);
    tags = tmask_chrs;
    f = figure('Visible',false);
    strtags = cellstr(strcat('Chr', string(tags)));
    legend(strtags(:));
    ax=gca;
    h = pie(ax, 1:length(j), tags');
    lh = legend(strtags(:));
    % Delete non-patch objects
    hPatch = h(strcmp(get(h,'type'),'patch')); 
    hNotPatch = h(~strcmp(get(h,'type'),'patch')); 
    delete(hNotPatch)
    ax.Colormap = cols;
    % Remove Faces of patches
    set(hPatch,'Faces', NaN)
    lhPos = get(lh,'Position');
    lh.Units='normalized';
    f.Units ='normalized';
    fPos = get(f,'Position');
    fPos(3:4)=lhPos(3:4);
    set(gcf, 'Position',fPos);
    set(lh,'Position',[0 0 1 1])
    saveas(gcf,[RESULTS_DIR 'legend' bCId 'BC.svg'])
    close(f);
end

%%
for j=1:length(handles)
    p = [PNG_DIR_1 num2str(j) '.png'];
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
            [~, cinds] = find(chromosomes==inds);
            cols = cmap(cinds,:);
            annot_pos = linspace(0,1,sum(msk)+1);
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
    text = 'digraph G{\nbgcolor="#ffffffff" # RGBA (with alpha)';
    for j=1:length(handles)
            parP = handles{j}.parent;
            if parP ~= 0
                text =  [text, '\n\t', num2str(parP) '->' num2str(j) '[ splines="false",style="invis"]' ];
            end
    end

    for chrCnt=1:length(chromosomes)
        chr = chromosomes(chrCnt);
        chrColor = ['#' sprintf('%02X',round(255 * cmap(chrCnt,:)))];
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
                end
            end
        end
    end
    text = [text, '\n}\n'];
    fileID = fopen([GRAPHVIZ_DIR  'graph' bCId 'BC.gv'],'w');
    fprintf(fileID,text);
    fclose(fileID);
    graph_png_path = [RESULTS_DIR 'CcaCircular' bCId 'BC.svg'];
    system(['twopi -Tsvg -o ' graph_png_path ' ' GRAPHVIZ_DIR 'graph' bCId 'BC.gv -Granksep=2 -Gratio=0.5 -Gfontsize=20'])
end
cd(old_path)