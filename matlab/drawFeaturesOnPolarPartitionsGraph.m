function drawFeaturesOnPolarPartitionsGraph(featMats, featMatsIds, featsClassesNames, modality, resultsDir, reduction, recomputeParts, recomputeBaseParts)
if nargin<6
    reduction=1;
end
if nargin<7
    recomputeParts=0;
end
if nargin<8
    recomputeBaseParts=0;
end
clusterDir = ['../results/' modality '/hierarchicalClusteringDemo/STAGE00DATA/'];
clusterArray = load(['../results/' modality '/hierarchicalClusteringDemo/STAGE00DATA/' modality '_reduction' num2str(reduction) '/levels4/segmentation.mat']).clusterArray;
template = load([clusterDir modality '_reduction' num2str(reduction) '/levels4/input_info.mat']).preprocTemplate;
[f1, f2, handles] = paintClusters(clusterArray, template, 4, false,[],'white',colorcube(2));
close(f2);
f1.Visible='off';
rawDir = [resultsDir 'props/'];
if ~isfolder(rawDir), mkdir(rawDir); end
for i=1:length(handles)
    ret_path = [rawDir num2str(i) '.png'];
    if ~recomputeBaseParts && isfile(ret_path)
        continue
    end
    f1 = figure('visible','off');
    s = copyobj(handles{i}.handle,f1);
    set(s,'position',[0 0 1 1]);
    set(f1, 'color', 'none');
    set(s, 'color', 'none');
    print(f1, '-dpng', '-r300', ret_path);
    saveas(f1, [rawDir num2str(i) '.pdf']);
    close(f1);
end

for i=1:length(featMats)
    featMat = featMats{i};
    featId = char(featMatsIds(i));
    featResultsDir = [resultsDir,featId,'/'];
    if ~isfolder(featResultsDir), mkdir(featResultsDir); end
    propDirOut = [featResultsDir,'props/'];
    if ~isfolder(propDirOut), mkdir(propDirOut); end
    mask = squeeze(any(isfinite(featMat), 1));
    anyOccsPerPart = find(mask);
    cmap = brewermap(length(mask), 'Spectral');
    cols = cmap(anyOccsPerPart, :);
    % Step 1: Create legend
    f = figure('Visible',false);
    ax=gca;
    h = pie(ax, 1:length(anyOccsPerPart), anyOccsPerPart');
    lh = legend(featsClassesNames(anyOccsPerPart));
    % Delete non-patch objects
    hPatch = h(strcmp(get(h,'type'),'patch'));
    hNotPatch = h(~strcmp(get(h,'type'),'patch'));
    delete(hNotPatch)
    ax.Colormap = cols;
    % Remove Faces of patches
    set(hPatch,'Faces', NaN)
    set(lh,'Box','off')
    lh.Units='pixels';
    lhPos = get(lh,'Position');
    f.Units ='pixels';
    set(f, 'Position',lhPos);
    lh.Units='normalized';
    set(lh,'Position',[0 0 1 1])
    saveas(gcf,[featResultsDir 'legend.svg'])
    close(f);
    % Step 2: Create modified partitions figures with features inscribed as
    % status bars
    for j=1:length(handles)
        p = [rawDir num2str(j) '.png'];
        ret_path = [propDirOut num2str(j) '.png'];
        if ~recomputeParts && isfile(ret_path)
            continue
        end
        im = imread(p);
        im = im(ceil(0.15*size(im,1)):size(im,1),:,:);
        f = figure('Visible','off');
        imshow(im);
        annotation(f, "textbox",'interpreter','latex', 'Position', [0.1, 0.75, 0.8, 0.15], 'String', num2str(j), ...
            FontUnits='Normalized', FontSize=0.1, ...
            FontWeight='bold', ...
            HorizontalAlignment='center', VerticalAlignment='middle', LineStyle='none');
        msk = isfinite(featMat(j, :));
        if any(msk)
            cols = cmap(msk,:);
            inds = find(msk);
            annot_pos = linspace(0,1,sum(msk)+1);
            for s=1:length(annot_pos) - 1
                annotation(f,'textbox','interpreter','latex', 'Position', [annot_pos(s) 0.1 annot_pos(s+1) - annot_pos(s) .15], ...
                    'String', num2str(featMat(j, inds(s))), 'BackgroundColor', cols(s,:), ...
                    FontUnits='Normalized',FontSize=0.06, ...
                    FontWeight='bold', ...
                    HorizontalAlignment='center', VerticalAlignment='middle');
            end
        end
        saveas(f, ret_path)
        close(f)
    end
    % Step 3: create graph
    old_path = pwd;
    cd(featResultsDir)
    text = 'digraph G{\nbgcolor="#ffffffff" # RGBA (with alpha)';
    for j=1:length(handles)
        parP = handles{j}.parent;
        if parP ~= 0
            text =  [text, '\n\t', num2str(parP) '->' num2str(j) '[ splines="false",style="invis"]' ];
        end
    end

    for chrCnt=1:length(anyOccsPerPart)
        chrColor = ['#' sprintf('%02X',round(255 * cmap(chrCnt,:)))];
        for j=1:length(handles)
            p = ['props/' num2str(j)]; % path relative to the svg file location
            text = [text, '\n\t', num2str(j) '[shape=plaintext, image="' ...
                p '.png",label="",fixedsize=true,width=4,height=2,fontcolor="' chrColor '"]' ];
        end
    end
    text = [text, '\n}\n'];
    graphvizDir = ['graphviz/'];
    if ~isfolder(graphvizDir), mkdir(graphvizDir); end
    fileID = fopen([graphvizDir  'graph.gv'],'w');
    fprintf(fileID,text);
    fclose(fileID);
    graph_png_path = ['CcaCircular.svg'];
    system(['twopi -Tsvg -o ' graph_png_path ' ' graphvizDir 'graph.gv -Granksep=2 -Gratio=0.5 -Gfontsize=20'])
    cd(old_path)
end
end