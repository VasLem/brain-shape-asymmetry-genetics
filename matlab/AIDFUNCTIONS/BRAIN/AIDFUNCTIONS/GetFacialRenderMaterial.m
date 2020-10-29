%% ASSEMBE RENDER MATERIAL FOR FACIAL SCAN

in = load('/uz/data/avalok/mic/tmp/pclaes4/karlijneShare/GWASDATA/PHENOTYPES/RV_0727SYM.mat');
facerender.ULABELS = label;
facerender.UHI = HierarchicalInterface;
facerender.UMASK = ones(1,facerender.UHI.nLC);
in = load('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/AM/Shape3DFormat/FacialTemplate_20171201.mat');

RefScan = clone(in.Template);
RefScan.UV = [];
RefScan.TextureMap = [];
RefScan.ViewMode = 'Solid';
facerender.RefScan = RefScan;
facerender.viewval = [0 90];
save('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/AM/Shape3DFormat/facerender.mat','facerender');

