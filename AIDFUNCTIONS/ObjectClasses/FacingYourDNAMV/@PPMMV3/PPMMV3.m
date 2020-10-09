classdef PPMMV3 <superClassLight
% PROPERTIES
    properties % GENERAL INTERFACING
        Mode = [];
        AssociationTest = 'CCA';%Test can be CCA or MANCOVAN
        RIPMethod = 'plsr';
        RIPK = 10;
    end
    properties (Dependent = true)
        nLev;
        PosClassBal;
    end
    properties (Hidden = true)
       HI;
    end
    properties (Hidden = true, Dependent = true)
       nModules;
       nS;
       AvgX;
       StdX;
       PercNan;
       PercNanXREG;
    end
    properties % ASSOCIATION TESTING
       PCAPerc = 99.9;
       Center = false;
       CorrectY = true;
       F;
       pF;
       R;
       DF1;
       DF2;
       ReplicationTest;
       ReplicationRIP;
       ClassifierTest;
    end
    properties (Dependent = true)
    end
    properties (Hidden = true)
        CCADIM;
        SELECT;
        MODSELECT;
        ACTION;
        PATH;
    end
    properties (Hidden = true, Dependent = true)
       XMOD; % X used in model (excluding nan values in XREG)
       XMODInd; % Index of non nan XREG values
       nPATH;
    end
    properties % REPLICATION TESTING
    end
    properties (Dependent = true)
    end
    properties (Hidden = true)
    end
    properties(Hidden = true, Dependent = true)  
    end
    properties % PHENOMIC PREDICTOR
       Predict = 'SOFT';
       HitType = 'pT';
       ModIndex = [];
    end
    properties (Dependent = true)
    end
    properties (Hidden = true)
        fCount = 0;
    end
    properties(Hidden = true, Dependent = true)  
    end
    properties (Abstract = true) % ABSTRACT PROPERTIES
       X;
       XREG;
       %Predictor;
    end
% OBJECT METHODS    
    methods % CONSTRUCTOR
       function obj = PPMMV3(varargin)
            obj = obj@superClassLight(varargin{:});
            obj.HI = HierarchicalInterface;
        end
    end
    methods % GENERAL GETTING
       function out = get.HI(obj)
           out = obj.HI;
           if~superClassLight.isH(out), out = []; end
       end
       function out = get.nLev(obj)
           if isempty(obj.HI), out = []; return; end
           out = obj.HI.nL;
       end
       function out = get.nModules(obj)
           if isempty(obj.HI), out = []; return; end
           out = obj.HI.nLC;
       end
       function out = get.nS(obj)
           out = length(obj.X); 
       end
       function out = get.XMODInd(obj)
           out = find(~isnan(obj.XREG)); 
       end
       function out = get.XMOD(obj)
           if isempty(obj.X), out = []; return; end
           out = double(obj.XREG(obj.XMODInd));
       end
       function out = get.PercNan(obj)
           if isempty(obj.X), out = []; return; end
           out = round((sum(isnan(obj.X))/obj.nS)*100);
       end
       function out = get.PercNanXREG(obj)
           if isempty(obj.XREG), out = []; return; end
           out = round((sum(isnan(obj.XREG))/obj.nS)*100);
       end
       function out = get.AvgX(obj)
            out = nanmean(obj.X);
       end
       function out = get.StdX(obj)
            out = nanstd(obj.X);
       end
       function out = get.PosClassBal(obj)
           if isempty(obj.X), out = 0; return; end
           tmp = obj.XREG;
           tmp(isnan(tmp)) = [];
           Tind = tmp==obj.PosClass;
           out = sum(Tind)/length(tmp);
       end
       function out = get.nPATH(obj)
                out = size(obj.PATH,1); 
       end
    end
    methods % GENERAL SETTING
       function obj = set.nLev(obj,in)
           if isempty(obj.HI), obj.HI=HierarchicalInterface;end
           obj.HI.L = single(in);
       end
    end
    methods % DISCOVERY ANALYSIS
       function out = runDiscoveryTest(obj,DepVar,COV,varargin)
           switch lower(obj.Mode)
               case {'modular' 'hierarchical'}
                   out = discoveryModHier(obj,DepVar,COV,obj.Mode,varargin{:});
               case 'propagate'
                   out = discoveryProp(obj,DepVar,COV,varargin{:});
               case {'univariate' 'multivariate'}
                   out = discoveryUniMulti(obj,DepVar,COV,obj.Mode);
               case 'select'
                   out = discoverySel(obj,DepVar,COV,varargin{:});
               case 'ga'
                   out = discoveryGA(obj,DepVar,COV,varargin{:});    
           end
       end
       function out = discoveryUniMulti(obj,DepVar,COV,mode)
                if iscell(DepVar), DepVar = getCellData(obj.HI,DepVar,1);end% FULL FACE PCs
                switch lower(mode)
                    case 'multivariate'
                       [out.F,out.pF] = parAssociationTest(obj,obj.XREG,DepVar,COV,[]);% MULTIVARIATE
                    case 'univariate'
                       [out.F,out.pF] = PPMMV3.parRegUnivariate(obj.XREG,DepVar,COV);% UNIVARIATE
                end
                obj.pF = out.pF;obj.F = out.F;
       end
       function out = discoveryModHier(obj,DepVar,COV,mode,varargin)
           Input = find(strcmpi(varargin, 'display'));
           if isempty(Input),display = false;else display = varargin{Input+1};end
           F = nan*zeros(1,obj.nModules);
           pF = nan*zeros(1,obj.nModules);
           R = nan*zeros(1,obj.nModules);
           DF1 = nan*zeros(1,obj.nModules);
           DF2 = nan*zeros(1,obj.nModules);
           B = cell(1,obj.nModules);
           DIM = zeros(1,obj.nModules);
           STATS = cell(1,obj.nModules);
           %warning off;
           for i=1:obj.nModules
              % i=1;
              forDepVar = [];
              switch lower(mode)
                  case 'modular'
                      forDepVar = getCellData(obj.HI,DepVar,i);
                      pcaperc = [];
                  case 'hierarchical'
                      modindex = [i getAllChildren(obj.HI,i)];
                      forDepVar = getCellData(obj.HI,DepVar,modindex); 
                      if ~isempty(obj.PCAPerc)&&~(obj.HI.Levels(i)==obj.nLev)% No reduction on lowest level 
                         pcaperc = obj.PCAPerc;
                      else
                         pcaperc = [];
                      end
                  otherwise
              end
              [F(i),pF(i),R(i),B{i},DF1(i),DF2(i),STATS{i}] =  parAssociationTest(obj,obj.XREG,forDepVar,COV,pcaperc); %#ok<*PFBNS>
           end
           %pF(pF==0) = 1*10^-50;
           out.F = F;out.pF = pF;out.R = R;out.DF1 = DF1;out.DF2 = DF2;out.B = B;out.STATS = STATS;
           if display
              plotHierValues3Dv2(-log10(out.pF),'crit',-log10(5*10^-8),'range',[0 max([-log10(5*10^-8) min(-log10(out.pF) -log10(1*10^-10))])],'title',num2str(min(out.pF)));
           end
           obj.F = F;obj.pF = pF;obj.R = R;obj.DF1 = DF1;obj.DF2 = DF2;
       end
       function out = discoveryProp(obj,DepVar,COV,varargin)
           Input = find(strcmpi(varargin, 'display'));
           if isempty(Input),display = false;else display = varargin{Input+1};end
           Input = find(strcmpi(varargin, 'displayiter'));
           if isempty(Input),displayiter = false;else displayiter = varargin{Input+1};end
           Input = find(strcmpi(varargin, 'pstart'));
           if isempty(Input),pStart = 0.05;else pStart = varargin{Input+1};end
           Input = find(strcmpi(varargin, 'pLock'));
           if isempty(Input),pLock = 0;else pLock = varargin{Input+1};end
           Input = find(strcmpi(varargin, 'use'));
           if isempty(Input),USE = ones(1,obj.nModules);else USE = varargin{Input+1};end
           if sum(USE)==0, return;end
           SELECT = zeros(1,obj.nModules);
           LOCK = zeros(1,obj.nModules);
           ACTION = zeros(1,obj.nModules);
           % POSSIBLE ACTION CODES
           % 3 = Fully Integrated
           % -3 = Children only Integrated
           % 2 = best child Integrated
           % -2 = worst child Integrated
           % 1 = best child propagated
           % -1 = Modular
           % 0 = lowest level modules
           actions = -3:1:3;
           MODSELECT = zeros(obj.nModules,obj.nModules);
           F = zeros(1,obj.nModules);
           pF = ones(1,obj.nModules);
           % START WITH LOWEST LEVEL
           L = obj.nLev;nC = 2^(L-1);
           modindex = LC2Ind(obj.HI,L*ones(1,nC),1:nC);
           for i=1:nC
               if USE(modindex(i))==0,continue;end
               forDepVar = getCellData(obj.HI,DepVar,modindex(i));
               [F(modindex(i)),pF(modindex(i))] =  PPMMV3.parCCA(obj.XREG,forDepVar,COV,[]);
               SELECT(modindex(i)) = pF(modindex(i)) <= pStart;
               MODSELECT(modindex(i),modindex(i)) = SELECT(modindex(i));
               ACTION(modindex(i)) = actions(3);
           end
           if displayiter
               plotrow = 1;
               f = figure;f.Color = [1 1 1];
               %f.Position = [1          41        1920        1084];
               plotaxes = cell(2,6);
               for pp=1:1:3
                   plotaxes{1,pp} = subplot(2,3,pp); 
               end
               for pp=1:1:3
                   plotaxes{2,pp} = subplot(2,3,pp+3); 
               end
               plotHierValues3Dv2(-log10(pF),'crit',-log10(5*10^-8),'peer',plotaxes{plotrow,1});
               %title(axpvalue,['LEVEL: ' num2str(L) ' PVALUE']);
               plotHierValues3Dv2(SELECT,'range',[0 1],'peer',plotaxes{plotrow,2});
               %title(axsel,['LEVEL: ' num2str(L) ' SELECT']);
               plotHierValues3Dv2(ACTION,'range',[-3 3],'peer',plotaxes{plotrow,3});
               %title(axaction,['LEVEL: ' num2str(L) ' SELECT']);
               drawnow;
               plotrow= plotrow+1;
               if plotrow==3, plotrow = 1;end
           end
           % WORK YOUR WAY UP
           for L=obj.nLev-1:-1:1
               % L = 5;
               nC = 2^(L-1);
               modindex = LC2Ind(obj.HI,L*ones(1,nC),1:nC);
               for i=1:nC
                   % i=1
                  % INITIALIZE FROM CHILDREN
                       ind = modindex(i);
                       if USE(ind)==0,continue;end
                       children = getChildren(obj.HI,ind);
                       pChildren = pF(children);
                       [bestp,bestsubind] = min(pChildren);
                       worstsubind = setdiff(1:2,bestsubind);
                       bestChild = children(bestsubind);
                       worstChild = children(worstsubind);
                       childrencombined = false;
                       % single branch detection
                       ChildrenAlive = ones(1,2);
                       for r=1:1:2
                           allchildren = [children(r) getAllChildren(obj.HI,children(r))];
                           if sum(SELECT(allchildren))==0
                              ChildrenAlive(r) = 0;
                           end
                       end
                  % COMBINED CHILDREN ANALYSIS
                       if sum(ChildrenAlive) == 2
                           allchildren = getAllChildren(obj.HI,ind);
                           index = find(SELECT(allchildren)==1);
                           setindex = allchildren(index);
                           forDepVar = getCellData(obj.HI,DepVar,setindex);
                           [combF,combpF] =  PPMMV3.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc);
                           if combpF<bestp
                              bestp = combpF;
                              childrencombined = true;
                              %disp([num2str(i) ': CHILDREN COMBINED']);
                           end
                       end
                  % HIERARCHICAL ANALYSIS 
                       if sum(ChildrenAlive) == 2%
                          setindex = [ind allchildren(index)]; 
                          forDepVar = getCellData(obj.HI,DepVar,setindex);  
                          [tmpF,tmppF] =  PPMMV3.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc); 
                          if tmppF<bestp
                              SELECT(ind) = 1;
                              F(ind) = tmpF;
                              pF(ind) = tmppF;
                              MODSELECT(ind,setindex) = 1;
                              ACTION(ind) = actions(7);
                              if displayiter, disp([num2str(i) ': FULLY INTEGRATED']);end
                              continue; % DONE YOU HAVE IMPROVED THE P-VALUE
                          end
                       end
                  % MODULAR ANALYSIS 
                       forDepVar = getCellData(obj.HI,DepVar,ind);
                       [modF,modpF] =  PPMMV3.parCCA(obj.XREG,forDepVar,COV,[]); 
                       if (modpF<bestp)||(sum(ChildrenAlive)==0)
                          if modpF<=pStart; 
                             SELECT(ind) = 1;
                          else
                             SELECT(ind) = 0;                          
                          end
                          F(ind) = modF;
                          pF(ind) = modpF;
                          MODSELECT(ind,ind) = SELECT(ind);
                          ACTION(ind) = actions(3);
                          index = getAllChildren(obj.HI,ind);
                          SELECT(index) = 0;
                          if displayiter, disp([num2str(i) ': MODULAR']);end
                          continue;
                       end
                  % BRANCH WITH THE BEST CHILD;
                       if ChildrenAlive(bestsubind)
                           allchildren = [bestChild getAllChildren(obj.HI,bestChild)];
                           index = find(SELECT(allchildren)==1);
                           setindex = [ind allchildren(index)];
                           forDepVar = getCellData(obj.HI,DepVar,setindex);  
                           [tmpF,tmppF] =  PPMMV3.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc);
                           if tmppF<bestp
                              SELECT(ind) = 1; %#ok<*PROPLC>
                              F(ind) = tmpF;
                              pF(ind) = tmppF;
                              MODSELECT(ind,setindex) = 1;
                              ACTION(ind) = actions(6);
                              if displayiter, disp([num2str(i) ': BEST CHILD INTEGRATED']);end
                              % remove other branch from select
                              allchildren = [worstChild getAllChildren(obj.HI,worstChild)];
                              SELECT(allchildren) = 0;
                              %if L==1, disp(num2str(tmppF));end
                              continue; % DONE YOU HAVE IMPROVED THE P-VALUE
                           end
                       end
                   % BRANCH WITH THE BEST CHILD;
                      if ChildrenAlive(worstsubind)
                           allchildren = [worstChild getAllChildren(obj.HI,worstChild)];
                           index = find(SELECT(allchildren)==1);
                           setindex = [ind allchildren(index)];
                           forDepVar = getCellData(obj.HI,DepVar,setindex);  
                           [tmpF,tmppF] =  PPMMV3.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc);
                           if tmppF<bestp
                              SELECT(ind) = 1;
                              F(ind) = tmpF;
                              pF(ind) = tmppF;
                              MODSELECT(ind,setindex) = 1;
                              ACTION(ind) = actions(2);
                              if displayiter, disp([num2str(i) ': WORST CHILD INTEGRATED']);end
                              % remove other branch from select
                              allchildren = [bestChild getAllChildren(obj.HI,bestChild)];
                              SELECT(allchildren) = 0;
                              %if L==1, disp(num2str(tmppF));end
                              continue; % DONE YOU HAVE IMPROVED THE P-VALUE
                           end
                      end
                   % REMOVE MODULE
                       SELECT(ind) = 0;
                       if childrencombined
                          F(ind) = combF;
                          pF(ind) = combpF;
                          ACTION(ind) = actions(1);
                          allchildren = getAllChildren(obj.HI,ind);
                          MODSELECT(ind,allchildren(SELECT(allchildren)==1)) = 1;
                          if displayiter, disp([num2str(i) ': CHILDREN INTEGRATED']); end
                       else
                          F(ind) = F(bestChild); 
                          pF(ind) = pF(bestChild);
                          ACTION(ind) = actions(5);
                          SELECT(worstChild) = 0;
                          allchildren = [bestChild getAllChildren(obj.HI,bestChild)];
                          MODSELECT(ind,allchildren(SELECT(allchildren)==1)) = 1;
                          if displayiter, disp([num2str(i) ': BEST CHILD PROPAGATED']);end
                          allchildren = [worstChild getAllChildren(obj.HI,worstChild)];
                          SELECT(allchildren) = 0;
                       end
               end
               if displayiter
                   plotHierValues3Dv2(-log10(pF),'crit',-log10(5*10^-8),'peer',plotaxes{plotrow,1});
                   %title(axpvalue,['LEVEL: ' num2str(L) ' PVALUE']);
                   plotHierValues3Dv2(SELECT,'range',[0 1],'peer',plotaxes{plotrow,2});
                   %title(axsel,['LEVEL: ' num2str(L) ' SELECT']);
                   plotHierValues3Dv2(ACTION,'range',[-3 3],'peer',plotaxes{plotrow,3});
                   %title(axaction,['LEVEL: ' num2str(L) ' SELECT']);
                   drawnow;
                   plotrow= plotrow+1;
                   if plotrow==3, plotrow = 1;end
                   pause;
               end
           end
           out.F = F;out.pF = pF;out.SELECT = SELECT;out.MODSELECT = MODSELECT;out.ACTION = ACTION;
           obj.F = F;obj.pF = pF;
           obj.SELECT = SELECT;obj.MODSELECT = MODSELECT;obj.ACTION = ACTION;
           if display
               f = figure;f.Color = [1 1 1];
               %f.Position = [1          41        1920        1084];
               ax1 = subplot(1,3,1);
               plotHierValues3Dv2(-log10(out.pF),'crit',-log10(5*10^-8),'range',[0 max([-log10(5*10^-8) -log10(out.pF)])],'title',num2str(min(out.pF)),'peer',ax1);   
               ax2 = subplot(1,3,2);
               plotHierValues3Dv2(SELECT,'range',[0 1],'peer',ax2,'title','Selection');
               ax3 = subplot(1,3,3);
               plotHierValues3Dv2(ACTION,'range',[-3 3],'peer',ax3,'title','Action');
               drawnow;
           end
       end
       function out = discoverySel(obj,DepVar,COV,varargin)
           Input = find(strcmpi(varargin, 'display'));
           if isempty(Input),display = false;else display = varargin{Input+1};end
           Input = find(strcmpi(varargin, 'displayiter'));
           if isempty(Input),displayiter = false;else displayiter = varargin{Input+1};end
           Input = find(strcmpi(varargin, 'pSelect'));
           if isempty(Input),pSelect = 0.05;else pSelect = varargin{Input+1};end
           SELECT = zeros(1,obj.nModules);
           F = zeros(1,obj.nModules);
           pF = ones(1,obj.nModules);
           % START WITH LOWEST LEVEL
           L = obj.nLev;nC = 2^(L-1);
           modindex = LC2Ind(obj.HI,L*ones(1,nC),1:nC);
           for i=1:nC
               forDepVar = getCellData(obj.HI,DepVar,modindex(i));
               [F(modindex(i)),pF(modindex(i))] =  PPMMV3.parCCA(obj.XREG,forDepVar,COV,[]);
               SELECT(modindex(i)) = pF(modindex(i)) <= pSelect;
           end
           % WORK YOUR WAY UP
           for L=obj.nLev-1:-1:1
               % L = 5;
               nC = 2^(L-1);
               modindex = LC2Ind(obj.HI,L*ones(1,nC),1:nC);
               for i=1:nC
                  % i=1
                   ind = modindex(i);
                   % MODULAR TEST FIRST
                   setindex = ind;
                   forDepVar = getCellData(obj.HI,DepVar,setindex);
                   [F(ind),pF(ind)] =  PPMMV3.parCCA(obj.XREG,forDepVar,COV,[]);
                   SELECT(ind) = pF(ind) <= pSelect;
                   if SELECT(ind)==0, continue;end
                   allchildren = getAllChildren(obj.HI,ind);
                   index = find(SELECT(allchildren)==1);
                   if isempty(index), continue;end
                   setindex = [ind allchildren(index)];
                   forDepVar = getCellData(obj.HI,DepVar,setindex);
                   if length(setindex)>1
                      [F(ind),pF(ind)] =  PPMMV3.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc);
                   else
                      [F(ind),pF(ind)] =  PPMMV3.parCCA(obj.XREG,forDepVar,COV,[]);
                   end
                   %SELECT(ind) = pF(ind) <= pSelect;
               end
           end
           out.F = F;out.pF = pF;out.SELECT = SELECT;
           obj.F = F;obj.pF = pF;obj.SELECT = SELECT;
           if display
               f = figure;f.Color = [1 1 1];
               %f.Position = [1          41        1920        1084];
               ax1 = subplot(1,2,1);
               plotHierValues3Dv2(-log10(out.pF),'crit',-log10(5*10^-8),'range',[0 max([-log10(5*10^-8) -log10(out.pF)])],'title',num2str(min(out.pF)),'peer',ax1);   
               ax2 = subplot(1,2,2);
               plotHierValues3Dv2(SELECT,'range',[0 1],'peer',ax2,'title','Selection');
%                ax3 = subplot(1,3,3);
%                plotHierValues3Dv2(ACTION,'range',[-3 3],'peer',ax3,'title','Action');
               drawnow;
           end
       end
       function [F,pF,R,B,DF1,DF2,STATS] = parAssociationTest(obj,X,Y,C,pcaperc)
           switch lower(obj.AssociationTest)
               case 'cca'
                   [F,pF,B,~,R,DF1,DF2,STATS] =  PPMMV3.parCCA(X,Y,C,pcaperc,obj.CorrectY);
               case 'mancovan'
                   [F,pF] =  PPMMV3.parMANCOVAN(X,Y,C,pcaperc);
                   R = [];
                   DF1 = [];
                   DF2 = [];
               case 'rip'
                   [F,pF,R] =  PPMMV3.parRIP(X,Y,C,obj.CorrectY,obj.RIPMethod,obj.RIPK);
                   DF1 = [];
                   DF2 = [];
           end
       end
       function getPropPaths(obj,pCrit)
           if isempty(obj.pF), obj.PATH = []; return;end
           if obj.pF(1)>pCrit, obj.PATH = []; return;end
           path = obj.SELECT;% THE STRONGEST PATH
           modindex = getAllModulesFromSelection(obj.HI,path);
           index = find(obj.pF<=pCrit);
           index = setdiff(index,modindex);
           if isempty(index), obj.PATH = path;return;end
           % multiple paths exist
           while ~isempty(index)
               path = [path;obj.MODSELECT(index(1),:)]; 
               modindex = union(modindex,getAllModulesFromSelection(obj.HI,path(end,:)));
               index = setdiff(index,modindex);
           end
           obj.PATH = path;
       end
       function out = discoveryGA(obj,DepVar,COV,varargin)
           Input = find(strcmpi(varargin, 'display'));
           if isempty(Input),display = false;else display = varargin{Input+1};end
           nvars = obj.nModules;
           popSize = 200;
           initPop = [eye(nvars);round(rand(popSize-nvars,nvars))];
           options = gaoptimset('Generations',50,'InitialPopulation',initPop,'PopulationSize',popSize,'PopulationType','bitstring','Vectorized','on','Display','iter');
           minfun = @(z)popDiscovery(obj,DepVar,COV,z);
           [x,~,~,~,population] = ga(minfun,nvars,[],[],[],[],[],[],[],options);
           population = unique(population,'rows');
           %plotHierValues3Dv2(x,'range',[0 1],'title','Selection');
           %plotHierValues3Dv2(sum(population,1),'title','Population');
           obj.SELECT = x;
           modindex = find(x);
           if length(modindex)==1
              pcaperc = [];
           else
              pcaperc = obj.PCAPerc;
           end
           forDepVar = getCellData(obj.HI,DepVar,modindex);
           [obj.F,obj.pF] =  parAssociationTest(obj,obj.XREG,forDepVar,COV,pcaperc); %#ok<*PFBNS>
           out.pF = obj.pF;
           out.F = obj.F;
           out.SELECT = obj.SELECT;
           %popeval = popDiscovery(obj,DepVar,COV,population);
           if display
              plotHierValues3Dv2(obj.SELECT,'range',[0 1],'title',num2str(min(out.pF)));
           end
       end
       function out = popDiscovery(obj,DepVar,COV,z)
                nZ = size(z,1);
                out = ones(nZ,1);
                parfor i=1:nZ
                    % i=1;
                    modindex = find(z(i,:));
                    if isempty(modindex)
                       pcaperc = []; 
                       continue;
                    elseif length(modindex)==1
                       pcaperc = [];
                    else
                       pcaperc = obj.PCAPerc;
                    end
                    forDepVar = getCellData(obj.HI,DepVar,modindex);
                    [~,out(i)] =  parAssociationTest(obj,obj.XREG,forDepVar,COV,pcaperc);
                    %[~,out(i)] =  PPMMV3.parCCA(obj.XREG,forDepVar,COV,pcaperc); %#ok<*PFBNS>
                end
       end
    end
    methods % REPLICATION ANALYSIS
        function out = runReplicationTest(obj,DepVar,COV,RX,RDepVar,RCOV,varargin)
            switch lower(obj.Mode)
               case 'modular'
                   [out,obj.ReplicationRIP] = replicationMod(obj,DepVar,COV,RX,RDepVar,RCOV,varargin{:});
               case 'hierarchical'
                   out = replicationHier(obj,DepVar,COV,RX,RDepVar,RCOV,varargin{:});    
               case 'propagate'
                   out = replicationProp(obj,DepVar,COV,RX,RDepVar,RCOV,varargin{:});
               case 'univariate'
                   out = replicationUni(obj,RX,RDepVar,RCOV,varargin{:});    
               case 'multivariate'
                   out = replicationMulti(obj,DepVar,COV,RX,RDepVar,RCOV,varargin{:});
               case 'select'
                   out = replicationSel(obj,DepVar,COV,RX,RDepVar,RCOV,varargin{:});
               case 'ga'
                   out = replicationGA(obj,DepVar,COV,RX,RDepVar,RCOV,varargin{:});
            end
            obj.ReplicationTest = out;
        end
        function out = replicationUni(obj,RX,RDepVar,RCOV,varargin)
                 Input = find(strcmpi(varargin, 'pcrit'));
                 if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
                 if sum(obj.pF<=pCrit)==0,out = [];return;end
                 if iscell(RDepVar), RDepVar = getCellData(obj.HI,RDepVar,1);end% FULL FACE PCs
                 axindex = find(obj.pF<=pCrit);
                 nV = size(RDepVar,2);
                 F = zeros(1,nV);pF = ones(1,nV);
                 Robj = clone(obj);Robj.X = RX;
                 for i=1:1:length(axindex)
                     [F(axindex(i)),pF(axindex(i))] = PPMMV3.parRegUnivariate(Robj.XREG,RDepVar(:,axindex(i)),RCOV);
                 end
                 out.F = F;out.pF = pF;
        end
        function out = replicationMulti(obj,DepVar,COV,RX,RDepVar,RCOV,varargin)
                Input = find(strcmpi(varargin, 'pcrit'));
                if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
                Input = find(strcmpi(varargin, 'nrperm'));
                if isempty(Input), t = 0; else t = varargin{Input+1};end
                if sum(obj.pF<=pCrit)==0,out = [];return;end
                if iscell(DepVar), DepVar = getCellData(obj.HI,DepVar,1);end% FULL FACE PCs
                if iscell(RDepVar), RDepVar = getCellData(obj.HI,RDepVar,1);end% FULL FACE PCs
                Robj = clone(obj);Robj.X = RX;
                out = PPMMV3.combinedMVReplicationTest(obj.XREG,DepVar,COV,Robj.XREG,RDepVar,RCOV,t,[]);
        end
        function [out,ripout] = replicationMod(obj,DepVar,COV,RX,RDepVar,RCOV,varargin)
                % READING VARIABLE INPUT
                Input = find(strcmpi(varargin, 'pcrit'));
                if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
                if sum(obj.pF<=pCrit)==0,out = [];return;end
                Input = find(strcmpi(varargin, 'nrperm'));
                if isempty(Input), t = 100; else t = varargin{Input+1};end
                Input = find(strcmpi(varargin, 'display'));
                if isempty(Input), display = false; else display = varargin{Input+1};end
                Input = find(strcmpi(varargin, 'angletable'));
                if isempty(Input), angletable = []; else angletable = varargin{Input+1};end
                % INITIALIZE
                %Robj = clone(obj);
                %Robj.X = RX;
                %RX = getXREG(obj,RX); NOTE THAT THIS OBJECT CANNOT TAKE
                %CARE OF THINGS!
                res = cell(1,obj.nModules);
                ripout = cell(1,obj.nModules);
                todo = zeros(1,obj.nModules);
                todo(obj.pF<=pCrit) = 1;
                parfor i=1:obj.nModules
                    if ~todo(i), continue; end
                    forDepVar = getCellData(obj.HI,DepVar,i);
                    forRDepVar = getCellData(obj.HI,RDepVar,i);
                    %res{i} = PPMMV3.combinedMVReplicationTest(obj.XREG,forDepVar,COV,RX,forRDepVar,RCOV,t,[],angletable);
                    [res{i},ripout{i}] = PPMMV3.FullReplicationTest(obj.XREG,forDepVar,COV,RX,forRDepVar,RCOV,[],obj.CorrectY,angletable,t);
                end
                out = PPMMV3.CellList2Structure(res);
        end 
        function out = replicationHier(obj,DepVar,COV,RX,RDepVar,RCOV,varargin)
                % READING VARIABLE INPUT
                Input = find(strcmpi(varargin, 'pcrit'));
                if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
                if sum(obj.pF<=pCrit)==0,out = [];return;end
                Input = find(strcmpi(varargin, 'nrperm'));
                if isempty(Input), t = 100; else t = varargin{Input+1};end
                Input = find(strcmpi(varargin, 'display'));
                if isempty(Input), display = false; else display = varargin{Input+1};end
                % INITIALIZE
                Robj = clone(obj);
                Robj.X = RX;
                res = cell(1,obj.nModules);
                todo = zeros(1,obj.nModules);
                todo(obj.pF<=pCrit) = 1;
                parfor i=1:obj.nModules
                   if ~todo(i), continue; end
                   if ~isempty(obj.PCAPerc)&&~(obj.HI.Levels(i)==obj.nLev)% No reduction on lowest level 
                        pcaperc = obj.PCAPerc;
                   else
                        pcaperc = [];
                   end
                   useindex = [i getAllChildren(obj.HI,i)];
                   forDepVar = getCellData(obj.HI,DepVar,useindex);
                   forRDepVar = getCellData(obj.HI,RDepVar,useindex);
                   res{i} = PPMMV3.combinedMVReplicationTest(obj.XREG,forDepVar,COV,Robj.XREG,forRDepVar,RCOV,t,pcaperc);
                end
                out = PPMMV3.CellList2Structure(res);
                if display
                   SignT = double(out.RIPpT<=0.05) + 2*double(obj.pF<=5*10^-8);
                   SignF = double(out.RIPpF<=0.05) + 2*double(obj.pF<=5*10^-8); 
                   Angle = double(out.pA<=0.05) + 2*double(obj.pF<=5*10^-8); 
                   f = figure;f.Color = [1 1 1];
                   ax1 = subplot(3,2,1);
                   plotHierValues3Dv2(-log10(out.RIPpT),'crit',-log10(0.05),'range',[0 max([-log10(0.05) -log10(out.RIPpT)])],'peer',ax1,'title',num2str(min(out.RIPpT)));
                   ax2 = subplot(3,2,2);
                   plotHierValues3Dv2(SignT,'range',[0 3],'peer',ax2,'title','Sing. Test');
                   ax3 = subplot(3,2,3);
                   plotHierValues3Dv2(-log10(out.pA),'crit',-log10(0.05),'range',[0 max([-log10(0.05) -log10(out.pA)])],'peer',ax3,'title',num2str(min(out.pA)));
                   ax4 = subplot(3,2,4);
                   plotHierValues3Dv2(Angle,'range',[0 3],'peer',ax4,'title','Angle Test');
                   ax5 = subplot(3,2,5);
                   plotHierValues3Dv2(-log10(out.RIPpF),'crit',-log10(0.05),'range',[0 max([-log10(0.05) -log10(out.RIPpF)])],'peer',ax5,'title',num2str(min(out.RIPpF)));
                   ax6 = subplot(3,2,6);
                   plotHierValues3Dv2(SignF,'range',[0 3],'peer',ax6,'title','Sing. Test');
                   drawnow;
                end
        end
        function out = replicationProp(obj,DepVar,COV,RX,RDepVar,RCOV,varargin)
            % READING VARIABLE INPUT
            Input = find(strcmpi(varargin, 'pcrit'));
            if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
            if sum(obj.pF<=pCrit)==0,out = [];return;end
            Input = find(strcmpi(varargin, 'nrperm'));
            if isempty(Input), t = 100; else t = varargin{Input+1};end
            modindex = find(obj.SELECT);
            DepVar = getCellData(obj.HI,DepVar,modindex);
            RDepVar = getCellData(obj.HI,RDepVar,modindex);
            Robj = clone(obj);Robj.X = RX;
            [out.F,out.pF] =  PPMMV3.parCCA(Robj.XREG,RDepVar,RCOV,obj.PCAPerc);
            [out.A,out.pA] = PPMMV3.angleTest(obj.XREG,DepVar,COV,Robj.XREG,RDepVar,RCOV,t);
        end
        function out = replicationGA(obj,DepVar,COV,RX,RDepVar,RCOV,varargin)
            Input = find(strcmpi(varargin, 'pcrit'));
            if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
            if sum(obj.pF<=pCrit)==0,out = [];return;end
            Input = find(strcmpi(varargin, 'nrperm'));
            if isempty(Input), t = 100; else t = varargin{Input+1};end
            modindex = find(obj.SELECT);
            DepVar = getCellData(obj.HI,DepVar,modindex);
            RDepVar = getCellData(obj.HI,RDepVar,modindex);
            Robj = clone(obj);Robj.X = RX;
            [out.F,out.pF] =  PPMMV3.parCCA(Robj.XREG,RDepVar,RCOV,obj.PCAPerc);
            [out.A,out.pA] = PPMMV3.angleTest(obj.XREG,DepVar,COV,Robj.XREG,RDepVar,RCOV,t);
        end
        function out = replicationSel(obj,DepVar,COV,RX,RDepVar,RCOV,varargin)
                % READING VARIABLE INPUT
                Input = find(strcmpi(varargin, 'pcrit'));
                if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
                if sum(obj.pF<=pCrit)==0,out = [];return;end
                Input = find(strcmpi(varargin, 'nrperm'));
                if isempty(Input), t = 100; else t = varargin{Input+1};end
                Input = find(strcmpi(varargin, 'display'));
                if isempty(Input), display = false; else display = varargin{Input+1};end
                % INITIALIZE
                Robj = clone(obj);
                Robj.X = RX;
                modindex = find(obj.pF<=pCrit);
                out.F = zeros(1,obj.nModules);
                out.pF = ones(1,obj.nModules);
                out.A = zeros(1,obj.nModules);
                out.pA = ones(1,obj.nModules);
                
                tmpF = zeros(1,length(modindex));
                tmppF = ones(1,length(modindex));
                tmpA = zeros(1,length(modindex));
                tmppA = ones(1,length(modindex));
                parfor i=1:length(modindex);
                   if ~isempty(obj.PCAPerc)&&~(obj.HI.Levels(modindex(i))==obj.nLev)% No reduction on lowest level 
                        pcaperc = obj.PCAPerc;
                   else
                        pcaperc = [];
                   end
                   useindex = [modindex(i) getAllChildren(obj.HI,modindex(i))];
                   selindex = find(obj.SELECT(useindex));
                   if length(selindex)==1, pcaperc = []; end
                   forDepVar = getCellData(obj.HI,DepVar,useindex(selindex));
                   forRDepVar = getCellData(obj.HI,RDepVar,useindex(selindex));
                   [tmpF(i),tmppF(i)] =  PPMMV3.parCCA(Robj.XREG,forRDepVar,RCOV,pcaperc); %#ok<*PFBNS>
                   [tmpA(i),tmppA(i)] = PPMMV3.angleTest(obj.XREG,forDepVar,COV,Robj.XREG,forRDepVar,RCOV,t);
                end
                out.F(modindex) = tmpF;
                out.pF(modindex) = tmppF;
                out.A(modindex) = tmpA;
                out.pA(modindex) = tmppA;
                if display
                   Sign = double(out.pF<=0.05) + 2*double(obj.pF<=5*10^-8); 
                   Angle = double(out.pA<=0.05) + 2*double(obj.pF<=5*10^-8); 
                   f = figure;f.Color = [1 1 1];
                   ax1 = subplot(2,2,1);
                   plotHierValues3Dv2(-log10(out.pF),'crit',-log10(0.05),'range',[0 max([-log10(0.05) -log10(out.pF)])],'peer',ax1,'title',num2str(min(out.pF)));
                   ax2 = subplot(2,2,2);
                   plotHierValues3Dv2(Sign,'range',[0 3],'peer',ax2,'title','Sing. Test');
                   ax3 = subplot(2,2,3);
                   plotHierValues3Dv2(-log10(out.pA),'crit',-log10(0.05),'range',[0 max([-log10(0.05) -log10(out.pA)])],'peer',ax3,'title',num2str(min(out.pA)));
                   ax4 = subplot(2,2,4);
                   plotHierValues3Dv2(Angle,'range',[0 3],'peer',ax4,'title','Angle Test');
                   drawnow;
                end
        end
    end
    methods % CLASSIFIER ANALYSIS
        function out = runClassifierTest(obj,X,DepVar,RX,RDepVar,varargin)
            switch lower(obj.Mode)
               case 'modular'
                   out = classifierMod(obj,X,DepVar,RX,RDepVar,varargin{:});
               case 'hierarchical'
                   out = classifierHier(obj,X,DepVar,RX,RDepVar,varargin{:});
               case 'propagate'
                   out = classifierProp(obj,X,DepVar,RX,RDepVar,varargin{:});
               case 'multivariate'
                   out = classifierMulti(obj,X,DepVar,RX,RDepVar,varargin{:});    
               case 'univariate'
                   out = classifierUni(obj,X,DepVar,RX,RDepVar,varargin{:});
               case 'select'
                   out = classifierSel(obj,X,DepVar,RX,RDepVar,varargin{:});
               case 'ga'
                   out = classifierProp(obj,X,DepVar,RX,RDepVar,varargin{:});    
            end
            obj.ClassifierTest = out;
        end
        function out = classifierMod(obj,X,DepVar,RX,RDepVar,varargin)
                 Input = find(strcmpi(varargin, 'pcrit'));
                 if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'display'));
                 if isempty(Input),display = false;else display = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'combined'));
                 if isempty(Input),combined = false;else combined = varargin{Input+1};end
                 modindex = find(obj.pF<=pCrit);
                 if isempty(modindex),out = [];return;end
                 out = [];
                 switch combined
                     case true
                         DepVar = getCellData(obj.HI,DepVar,modindex);
                         trainClassifier(obj,X,DepVar,varargin{:});
                         RDepVar = getCellData(obj.HI,RDepVar,modindex);
                         out = testClassifier(obj,RX,RDepVar,varargin{:});
                         if ~isempty(out), out.CV = obj.CV; else return, end % added by Dzemila 18112016
                     case false
                         nV = length(modindex);
                         tmp = cell(1,nV);
                         out = cell(1,obj.nModules);
                         parfor i=1:nV
                             forobj = clone(obj);
                             forDepVar = getCellData(forobj.HI,DepVar,modindex(i));
                             trainClassifier(forobj,X,forDepVar,varargin{:});
                             forRDepVar = getCellData(forobj.HI,RDepVar,modindex(i));
                             tmp{i} = testClassifier(forobj,RX,forRDepVar,varargin{:});
                             tmp{i}.CV = forobj.CV;
                         end
                         out(modindex) = tmp;
                         if display
                            figure;hold on;grid on;grid minor;plot(0:0.1:1,0:0.1:1,'k-');
                            plot(0:0.1:1,1:-0.1:0,'k--');
                            cm = colormap('lines');
                            xlabel('FPR');ylabel('TPR');
                            str = cell(1,nV+2);
                            str{1} = '';
                            str{2} = '';
                            for i=1:nV
                                str{i+2} = [num2str(modindex(i)) ' ' num2str(out{modindex(i)}.pAUC)];
                                Color = cm(mod(modindex(i),size(cm,2)-1)+1,:);
                                plot(out{modindex(i)}.XX,out{modindex(i)}.YY,'Color',Color,'LineStyle','-','LineWidth',1.5);
                                %plot(out{modindex(i)}.XX,out{modindex(i)}.YY,'b-','LineWidth',1.5);
                                %plot(out{modindex(i)}.HX,out{modindex(i)}.HY,'r.','MarkerSize',20);   
                                %title(['G: ' num2str(out.HG) ' AUC: ' num2str(out.AUC) ' pAUC: ' num2str(out.pAUC)]);
                                
                            end
                            legend(str,'Location','SouthEast');
                            for i=1:nV
                                plot(out{modindex(i)}.HX,out{modindex(i)}.HY,'r.','MarkerSize',20);
                            end
                            drawnow;
                         end
                 end
                 if display
                    return;
                 end
        end
        function out = classifierHier(obj,X,DepVar,RX,RDepVar,varargin)
                 Input = find(strcmpi(varargin, 'pcrit'));
                 if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'display'));
                 if isempty(Input),display = false;else display = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'combined'));
                 if isempty(Input),combined = false;else combined = varargin{Input+1};end
                 modindex = find(obj.pF<=pCrit);
                 if isempty(modindex),out = [];return;end
                 switch combined
                     case true
                         modindex = expandWithChildren(obj.HI,modindex);
                         DepVar = getCellData(obj.HI,DepVar,modindex);
                         trainClassifier(obj,X,DepVar,varargin{:});
                         RDepVar = getCellData(obj.HI,RDepVar,modindex);
                         out = testClassifier(obj,RX,RDepVar,varargin{:});
                         out.CV = obj.CV;
                     case false
                         nV = length(modindex);
                         tmp = cell(1,nV);
                         out = cell(1,obj.nModules);
                         parfor i=1:nV
                             forobj = clone(obj);
                             tmpindex = [modindex(i) getAllChildren(forobj.HI,modindex(i))];
                             forDepVar = getCellData(forobj.HI,DepVar,tmpindex);
                             trainClassifier(forobj,X,forDepVar,varargin{:});
                             forRDepVar = getCellData(forobj.HI,RDepVar,tmpindex);
                             tmp{i} = testClassifier(forobj,RX,forRDepVar,varargin{:});
                             tmp{i}.CV = forobj.CV;
                         end
                         out(modindex) = tmp;
                         if display
                            figure;hold on;grid on;grid minor;plot(0:0.1:1,0:0.1:1,'k-');
                            plot(0:0.1:1,1:-0.1:0,'k--');
                            cm = colormap('lines');
                            xlabel('FPR');ylabel('TPR');
                            str = cell(1,nV+2);
                            str{1} = '';
                            str{2} = '';
                            counter = 1;
                            for i=1:nV
                                str{i+2} = [num2str(modindex(i)) ' ' num2str(out{modindex(i)}.pAUC)];
                                %Color = cm(mod(modindex(i),size(cm,2)-1)+1,:);
                                Color = cm(counter,:);
                                counter = counter +1;
                                plot(out{modindex(i)}.XX,out{modindex(i)}.YY,'Color',Color,'LineStyle','-','LineWidth',1.5);
                                %plot(out{modindex(i)}.XX,out{modindex(i)}.YY,'b-','LineWidth',1.5);
                                %plot(out{modindex(i)}.HX,out{modindex(i)}.HY,'r.','MarkerSize',20);   
                                %title(['G: ' num2str(out.HG) ' AUC: ' num2str(out.AUC) ' pAUC: ' num2str(out.pAUC)]);
                                
                            end
                            legend(str,'Location','SouthEast');
                            for i=1:nV
                                plot(out{modindex(i)}.HX,out{modindex(i)}.HY,'r.','MarkerSize',20);
                            end
                            drawnow;
                         end
                 end
                 if display
                    return;
                 end
        end
        function out = classifierSel(obj,X,DepVar,RX,RDepVar,varargin)
                 Input = find(strcmpi(varargin, 'pcrit'));
                 if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'display'));
                 if isempty(Input),display = false;else display = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'combined'));
                 if isempty(Input),combined = false;else combined = varargin{Input+1};end
                 modindex = find(obj.pF<=pCrit);
                 if isempty(modindex),out = [];return;end
                 switch combined
                     case true
                         modindex = expandWithChildren(obj.HI,modindex);
                         selindex = find(obj.SELECT(modindex));
                         DepVar = getCellData(obj.HI,DepVar,modindex(selindex));
                         trainClassifier(obj,X,DepVar,varargin{:});
                         RDepVar = getCellData(obj.HI,RDepVar,modindex(selindex));
                         out = testClassifier(obj,RX,RDepVar,varargin{:});
                         out.CV = obj.CV;
                     case false
                         nV = length(modindex);
                         tmp = cell(1,nV);
                         out = cell(1,obj.nModules);
                         parfor i=1:nV
                             forobj = clone(obj);
                             tmpindex = [modindex(i) getAllChildren(forobj.HI,modindex(i))];
                             selindex = find(obj.SELECT(tmpindex));
                             forDepVar = getCellData(forobj.HI,DepVar,tmpindex(selindex));
                             trainClassifier(forobj,X,forDepVar,varargin{:});
                             forRDepVar = getCellData(forobj.HI,RDepVar,tmpindex(selindex));
                             tmp{i} = testClassifier(forobj,RX,forRDepVar,varargin{:});
                             tmp{i}.CV = forobj.CV;
                         end
                         out(modindex) = tmp;
                         if display
                            figure;hold on;grid on;grid minor;plot(0:0.1:1,0:0.1:1,'k-');
                            plot(0:0.1:1,1:-0.1:0,'k--');
                            cm = colormap('lines');
                            xlabel('FPR');ylabel('TPR');
                            str = cell(1,nV+2);
                            str{1} = '';
                            str{2} = '';
                            counter = 1;
                            for i=1:nV
                                str{i+2} = [num2str(modindex(i)) ' ' num2str(out{modindex(i)}.pAUC)];
                                %Color = cm(mod(modindex(i),size(cm,2)-1)+1,:);
                                Color = cm(counter,:);
                                counter = counter +1;
                                plot(out{modindex(i)}.XX,out{modindex(i)}.YY,'Color',Color,'LineStyle','-','LineWidth',1.5);
                                %plot(out{modindex(i)}.XX,out{modindex(i)}.YY,'b-','LineWidth',1.5);
                                %plot(out{modindex(i)}.HX,out{modindex(i)}.HY,'r.','MarkerSize',20);   
                                %title(['G: ' num2str(out.HG) ' AUC: ' num2str(out.AUC) ' pAUC: ' num2str(out.pAUC)]);
                                
                            end
                            legend(str,'Location','SouthEast');
                            for i=1:nV
                                plot(out{modindex(i)}.HX,out{modindex(i)}.HY,'r.','MarkerSize',20);
                            end
                            drawnow;
                         end
                 end
                 if display
                    return;
                 end
        end
        function out = classifierProp(obj,X,DepVar,RX,RDepVar,varargin)
            Input = find(strcmpi(varargin, 'pcrit'));
            if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
            if sum(obj.pF<=pCrit)==0, out = []; return; end
            modindex = find(obj.SELECT);
            DepVar = getCellData(obj.HI,DepVar,modindex);
            RDepVar = getCellData(obj.HI,RDepVar,modindex);
            trainClassifier(obj,X,DepVar,varargin{:});
            out = testClassifier(obj,RX,RDepVar,varargin{:});
            out.CV = obj.CV;
        end
        function out = classifierGA(obj,X,DepVar,RX,RDepVar,varargin)
            Input = find(strcmpi(varargin, 'pcrit'));
            if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
            if sum(obj.pF<=pCrit)==0, out = []; return; end
            modindex = find(obj.SELECT);
            DepVar = getCellData(obj.HI,DepVar,modindex);
            RDepVar = getCellData(obj.HI,RDepVar,modindex);
            trainClassifier(obj,X,DepVar,varargin{:});
            out = testClassifier(obj,RX,RDepVar,varargin{:});
            out.CV = obj.CV;
        end
        function out = classifierMulti(obj,X,DepVar,RX,RDepVar,varargin)
            Input = find(strcmpi(varargin, 'pcrit'));
            if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
            if sum(obj.pF<=pCrit)==0, out = []; return; end
            if iscell(DepVar), DepVar = getCellData(obj.HI,DepVar,1);end% FULL FACE PCs
            if iscell(RDepVar), RDepVar = getCellData(obj.HI,RDepVar,1);end% FULL FACE PCs
            trainClassifier(obj,X,DepVar,varargin{:});
            out = testClassifier(obj,RX,RDepVar,varargin{:});
            out.CV = obj.CV;
        end
        function out = classifierUni(obj,X,DepVar,RX,RDepVar,varargin)
            Input = find(strcmpi(varargin, 'pcrit'));
            if isempty(Input),pCrit = 1;else pCrit = varargin{Input+1};end
            if sum(obj.pF<=pCrit)==0, out = []; return; end
            if iscell(DepVar), DepVar = getCellData(obj.HI,DepVar,1);end% FULL FACE PCs
            if iscell(RDepVar), RDepVar = getCellData(obj.HI,RDepVar,1);end% FULL FACE PCs
            out = cell(1,size(DepVar,2));
            axindex = find(obj.pF<=pCrit);
            nV = length(axindex);
            tmp = cell(1,nV);
            parfor i=1:nV
                forobj = clone(obj);
                trainClassifier(forobj,X,DepVar(:,axindex(i)),varargin{:});
                tmp{i} = testClassifier(forobj,RX,RDepVar(:,axindex(i)),varargin{:});
                tmp{i}.CV = forobj.CV;
            end
            out(axindex) = tmp;
        end
    end
    methods % EFFECT ANALYSIS
        function f = imageOverlaidPValue(obj,label,RefScan)
                 scan = clone(RefScan);
                 nL = size(label,1);
                 val = ones(nL,scan.nrV);
                 for l=1:1:nL
                    % l = 1
                    clusters = unique(label(l,:));
                    nCL = length(clusters);
                    for c=1:1:nCL
                        % c=1;
                        index = find(label(l,:)==clusters(c));
                        ind = LC2Ind(obj.HI,l,c);
                        val(l,index) = obj.pF(ind);
                    end   
                 end
                 val = min(val);
                 scan.Value = -log10(val);
                 scan.ColorMode = 'Indexed';
                 scan.Material = 'Dull';
                 f = viewer(scan);
                 f.BackgroundColor = [1 1 1];
                 f.SceneLightVisible = true;
                 colorbar('peer',f.RenderAxes);
                 set(f.RenderAxes,'clim',[0 max(scan.Value)]);
                 title(f.RenderAxes,[num2str(min(val))]);
        end
        function out = signModules(obj,pCrit)
            out = find(obj.pF<=pCrit);
        end
        function [pindex,lindex] = getRelModules(obj,pCrit)
            [~,pindex] = min(obj.pF);
            modindex = find(obj.pF<=pCrit);
            [~,maxl] = min(obj.HI.Levels(modindex));
            lindex = modindex(maxl);
        end
        function out = getShapeEffect(obj,DepVar,COV,ind,SPACE,label,t)
                 if nargin<8, t= 0; end
                 [l,c] = Ind2LC(obj.HI,ind);
                 lmindex = find(label(l,:)==c);
                 modspace = getCellData(obj.HI,SPACE,ind,'SPACE');
                 Y = getCellData(obj.HI,DepVar,ind);
                 % eliminating missing data
                 XREG = obj.XREG;
                 index = find(~isnan(XREG));
                 Y = Y(index,:);
                 XREG = XREG(index,:);
                 COV = COV(index,:);
                 % done
                 n = size(Y,1);
                 Y2 = Y.*repmat(modspace.EigStd',n,1);
                 LM = repmat(modspace.AvgVec,1,n) + modspace.EigVec*Y2';
                 %LM = reconstructTraining(modspace);
                 NX = PPMMV3.getResiduals(COV,XREG);
                 out = multiplePlsRegress(NX,LM',t,'shape');
                 B =  PPMMV3.PLSR(NX,Y,[]);
                 B = B'.*modspace.EigStd;
                 
                 scan = clone(modspace.Average);
                 
                 stepsize = 5;
                 vec = modspace.EigVec*B*stepsize;
                 diff = Vec2Struc(modspace,vec);
                 pscan = clone(scan);
                 pscan.Vertices = scan.Vertices+diff;
                 
                 mscan = clone(scan);
                 mscan.Vertices = scan.Vertices-diff;
                 mscan.ColorMode = 'Single';
                 mscan.SingleColor = [0.8 0.1 0.1];
                 
                 out.Distances = vDistances(pscan,mscan);
                 out.Differences = vDifferences(pscan,mscan);
                 out.NormalDistances = vNormalDistances(mscan,pscan);
                 out.lmindex = lmindex;
                 %out.Differences(1,:) = out.Differences(1,:).*sign(scan.Vertices(1,:));
            
        end
        function [f,outpscan,outmscan] = showEffect(obj,DepVar,COV,modindex,SPACE,type)
            if nargin < 6, type = 'reduced';end
            nMod = length(modindex);
          
            for m=1:1:nMod
                %m=1
                Y = getCellData(obj.HI,DepVar,modindex(m));
                forspace = getCellData(obj.HI,SPACE,modindex(m));
                B =  PPMMV3.PLSR(obj.XREG,Y,COV);
                B = B'.*forspace.EigStd;
                %B = B';
                switch type
                    case 'reduced'
                
                        %B = B';
                        %showShapeEffect(forspace,10*B',type);

                        scan = clone(forspace.Average);
                        scan.Material = 'Dull';
                        vec = forspace.EigVec*B;
                        diff = Vec2Struc(forspace,vec);
                        scan.Value = sqrt(sum(diff.^2));
                        scan.ColorMode = 'Indexed';
                        v1 = viewer(scan);
                        v1.Tag = ['Module ' num2str(modindex(m))];

                        %viewer(scan);
                        stepsize = 5;
                        vec = forspace.EigVec*B*stepsize;
                        diff = Vec2Struc(forspace,vec);
                        pscan = clone(scan);
                        pscan.Vertices = scan.Vertices+diff;
                        pscan.ColorMode = 'Single';
                        pscan.SingleColor = [0.1 0.8 0.1];
                        pscan.ViewMode = 'Wireframe';
                        v2 = viewer(pscan);

                        mscan = clone(scan);
                        mscan.Vertices = scan.Vertices-diff;
                        mscan.ColorMode = 'Single';
                        mscan.SingleColor = [0.8 0.1 0.1];
                        mscan.ViewMode = 'Wireframe';
                        mscan.Axes = v2.RenderAxes;
                        mscan.Visible = true;
                        v2.Tag = ['Module ' num2str(modindex(m))];
                
                    case 'full'
                        
                        scan = clone(forspace.Average);
                        
%                         for i=1:1:3
%                             scan.Vertices(i,:) = scan.Vertices(i,:)+100;
%                         end
                        
                        scan.Material = 'Dull';
                        vec = forspace.EigVec*B;
                        diff = Vec2Struc(forspace,vec);
                        scan.Value = sqrt(sum(diff.^2));
                        scan.ColorMode = 'Indexed';
                        
                        %viewer(scan);
                        stepsize = 10;
                        disp(num2str(stepsize));
                        vec = forspace.EigVec*B*stepsize;
                        diff = Vec2Struc(forspace,vec);
                        pscan = clone(scan);
                        pscan.Vertices = scan.Vertices+diff;
                        pscan.ColorMode = 'Single';
                        pscan.SingleColor = [0.1 0.8 0.1];
                        pscan.ViewMode = 'Wireframe';
                        

                        mscan = clone(scan);
                        mscan.Vertices = scan.Vertices-diff;
                        mscan.ColorMode = 'Single';
                        mscan.SingleColor = [0.8 0.1 0.1];
                        mscan.ViewMode = 'Wireframe';
                        
                        f = figure;
                        f.Color = [1 1 1];
                        
%                         ax = subplot(2,3,6);axis equal;
%                         pscan.Axes = ax;pscan.Visible = true;
%                         mscan.Axes = ax;mscan.Visible = true;
                        
                        outpscan = clone(pscan);
                        outmscan = clone(mscan);

                        out = compareMorphs(pscan,mscan);
                        
                        ax = subplot(2,3,1);
                        newscan = clone(scan);
                        newscan.Value = out.Distances;
                        newscan.Axes = ax;
                        newscan.Visible = true;
                        camlight headlight;axis equal;colorbar;
                        title(['Module ' num2str(modindex(m))]);
                        xlabel('Magnitude');
                        ax.Visible = 'off';

                        ax = subplot(2,3,2);
                        newscan = clone(scan);
                        newscan.Value = out.NormalDistances;
                        maxval = max(abs(out.NormalDistances));
                        newscan.Axes = ax;
                        newscan.Visible = true;
                        camlight headlight;axis equal;colormap('jet');colorbar;
                        set(gca,'clim',[-1*maxval maxval]);
                        xlabel('Normal Distances');
                        title(obj.RS);
                        ax.Visible = 'off';

                        ax = subplot(2,3,3);
                        newscan = clone(scan);
                        newscan.Value = out.AreaRatios;
                        maxval = max(abs(out.AreaRatios));
                        newscan.Axes = ax;
                        newscan.Visible = true;
                        camlight headlight;axis equal;colormap('jet');colorbar;
                        set(gca,'clim',[-1*maxval maxval]);
                        xlabel('Area Ratios');
                        ax.Visible = 'off';


                        maxval = max([abs(out.signedCurvature1) abs(out.signedCurvature2)]);
                        maxval = max(diff(:));

                        ax = subplot(2,3,4);
                        newscan = clone(scan);
                        %newscan.Value = out.signedCurvature1;
                        newscan.Value = diff(1,:).*sign(scan.Vertices(1,:));
                        %maxval = max(abs(diff(1,:)));
                        newscan.Axes = ax;
                        newscan.Visible = true;
                        camlight headlight;axis equal;colormap('jet');colorbar;
                        set(gca,'clim',[-1*maxval maxval]);
                        xlabel('aa Curvature');
                        ax.Visible = 'off';

                        ax = subplot(2,3,5);
                        newscan = clone(scan);
                        %newscan.Value = out.signedCurvature2;
                        newscan.Value = diff(2,:);%.*sign(scan.Vertices(2,:));
                        %maxval = max(abs(diff(2,:)));
                        newscan.Axes = ax;
                        newscan.Visible = true;
                        camlight headlight;axis equal;colormap('jet');colorbar;
                        set(gca,'clim',[-1*maxval maxval]);
                        xlabel('AA Curvature');
                        ax.Visible = 'off';
                        
                        
                        ax = subplot(2,3,6);
                        newscan = clone(scan);
                        %newscan.Value = out.signedCurvatureDiff;
                        %maxval = max(abs(out.signedCurvatureDiff));
                        newscan.Value = diff(3,:);%.*sign(scan.Vertices(3,:));
                        %maxval = max(abs(diff(3,:)));
                        newscan.Axes = ax;
                        newscan.Visible = true;
                        camlight headlight;axis equal;colormap('jet');colorbar;
                        set(gca,'clim',[-1*maxval maxval]);
                        xlabel('Curvature Diff');
                        ax.Visible = 'off';
                end
                 
            end
            
            
        end
        function B = getMod11Effect(obj,DepVar,COV)
                Y = getCellData(obj.HI,DepVar,1);
%                 forspace = getCellData(obj.HI,SPACE,1);   
                X = PPMMV3.getResiduals(COV,obj.XREG);
                B =  PPMMV3.PLSR(X,Y,[]);     
        end
        function B = getModEffect(obj,DepVar,COV,modindex)
                Y = getCellData(obj.HI,DepVar,modindex);
%                 forspace = getCellData(obj.HI,SPACE,1);   
                X = PPMMV3.getResiduals(COV,obj.XREG);
                B =  PPMMV3.PLSR(X,Y,[]);
        end
    end
    methods % BEST CLASSIFIER
        function out = bestClassifierGA(obj,DepVar,RX,RDepVar,varargin)
           nvars = obj.nModules;
           popSize = 50;
           initPop = round(rand(popSize,nvars));
           options = gaoptimset('Generations',10,'InitialPopulation',initPop,'PopulationSize',popSize,'PopulationType','bitstring','Vectorized','on','Display','iter');
           minfun = @(z)popClassifier(obj,DepVar,z);
           obj.fCount = 0;
           [x,~,~,~,population] = ga(minfun,nvars,[],[],[],[],[],[],[],options);
           population = unique(population,'rows');
           obj.SELECT = x;
           modindex = find(x);
           FS = getCellData(obj.HI,DepVar,modindex);
           trainClassifier(obj,obj.GT,FS,varargin{:});
           RFS = getCellData(obj.HI,RDepVar,modindex);
           out = testClassifier(obj,RX,RFS,varargin{:});
        end
        function out = popClassifier(obj,DepVar,z)
                nZ = size(z,1);
                out = ones(nZ,1);
                obj.fCount = obj.fCount+1;
                [X,COST,W,ind,~,~] = prepTrainClassifier(obj,obj.X);
                rng(obj.fCount);K = cvpartition(length(ind),'KFold',3);
                parfor i=1:nZ
                    % i=1;
                    modindex = find(z(i,:));
                    FS = getCellData(obj.HI,DepVar,modindex);
                    tmp = CATMV3.trainSVMRBF(FS(ind,:),X(ind),COST,W(ind),K,obj.Optimize,obj.Standardize,obj.PosClass);
                    out(i) = nanmean(tmp.pAUC);
                end
        end
        function out = bestClassifierGAVal(obj,DepVar,RX,RDepVar,varargin)
           nvars = double(obj.nModules);
           popSize = 50;
           initPop = round(rand(popSize,nvars));
           options = gaoptimset('Generations',10,'InitialPopulation',initPop,'PopulationSize',popSize,'PopulationType','bitstring','Vectorized','on','Display','iter');
           minfun = @(z)popClassifierVal(obj,DepVar,RX,RDepVar,z);
           obj.fCount = 0;
           [x,~,~,~,population] = ga(minfun,nvars,[],[],[],[],[],[],[],options);
           population = unique(population,'rows');
           obj.SELECT = x;
           modindex = find(x);
           FS = getCellData(obj.HI,DepVar,modindex);
           trainClassifier(obj,obj.GT,FS,varargin{:});
           RFS = getCellData(obj.HI,RDepVar,modindex);
           out = testClassifier(obj,RX,RFS,varargin{:});
        end
        function out = popClassifierVal(obj,DepVar,ValX,ValDepVar,z)
                nZ = size(z,1);
                out = ones(nZ,1);
                obj.fCount = obj.fCount+1;
                parfor i=1:nZ
                    % i=1;
                    cobj = clone(obj);
                    modindex = find(z(i,:));
                    FS = getCellData(cobj.HI,DepVar,modindex);
                    trainClassifier(cobj,cobj.X,FS,'K',0,'optimize',false,'standardize',true);
                    FSVAL = getCellData(cobj.HI,ValDepVar,modindex);
                    tmpout = testClassifier(cobj,ValX,FSVAL,'display',false);
                    out(i) = tmpout.pAUC;
                end
        end
    end
    methods % REPORTING
        function out = reportAssociation(obj,varargin)
             switch lower(obj.Mode)
               case 'modular'
                   out = reportMod(obj,varargin{:});
               case 'hierarchical'
                   out = reportHier(obj,varargin{:});    
               case 'propagate'
                   out = reportProp(obj,varargin{:});
               case 'univariate'
                   out = reportUni(obj,varargin{:});    
                case 'multivariate'
                   out = reportMulti(obj,varargin{:});
             end
        end
        function out = reportUni(obj,varargin)
                 Input = find(strcmpi(varargin, 'pDis'));
                 if isempty(Input),pDis = 5*10^-8;else pDis = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'pRep'));
                 if isempty(Input),pRep = 5*10^-8;else pRep = varargin{Input+1};end
                 out.Discovery = sum(obj.pF<=pDis);
                 if ~isempty(obj.ReplicationTest)
                    out.Replication = sum(obj.ReplicationTest.pF<=pRep);
                 else
                    out.Replication = 0;
                 end
                 if ~isempty(obj.ClassifierTest)
                     tmp = zeros(1,length(obj.ClassifierTest));
                     for i=1:length(obj.ClassifierTest)
                        if isempty(obj.ClassifierTest{i}),continue;end
                        tmp(i) = obj.ClassifierTest{i}.pAUC<=pRep; 
                     end
                     out.Classifier = sum(tmp);
                 else
                     out.Classifier = 0;
                 end
                 out.ReplicationA = 0;
        end
        function out = reportMulti(obj,varargin)
                 Input = find(strcmpi(varargin, 'pDis'));
                 if isempty(Input),pDis = 5*10^-8;else pDis = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'pRep'));
                 if isempty(Input),pRep = 5*10^-8;else pRep = varargin{Input+1};end
                 out.Discovery = sum(obj.pF<=pDis);
                 if ~isempty(obj.ReplicationTest)
                    out.Replication = sum(obj.ReplicationTest.pF<=pRep);
                    out.ReplicationA = sum(obj.ReplicationTest.pA<=pRep);
                 else
                    out.Replication = 0;
                    out.ReplicationA = 0;
                 end
                 if ~isempty(obj.ClassifierTest)
                     out.Classifier = obj.ClassifierTest.pAUC<=pRep;
                 else
                     out.Classifier = 0;
                 end
        end
        function out = reportMod(obj,varargin)
                 Input = find(strcmpi(varargin, 'pDis'));
                 if isempty(Input),pDis = 5*10^-8;else pDis = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'pRep'));
                 if isempty(Input),pRep = 5*10^-8;else pRep = varargin{Input+1};end
                 out.Discovery = sum(obj.pF<=pDis);
                 if ~isempty(obj.ReplicationTest)
                    out.Replication = sum(obj.ReplicationTest.pF<=pRep);
                    out.ReplicationA = sum(obj.ReplicationTest.pA<=pRep);
                 else
                    out.Replication = 0;
                    out.ReplicationA = 0;
                 end
                 if ~isempty(obj.ClassifierTest)
                     tmp = zeros(1,length(obj.ClassifierTest));
                     for i=1:length(obj.ClassifierTest)
                        if isempty(obj.ClassifierTest{i}),continue;end
                        tmp(i) = obj.ClassifierTest{i}.pAUC<=pRep; 
                     end
                     out.Classifier = sum(tmp);
                 else
                     out.Classifier = 0;
                 end
        end
        function out = reportHier(obj,varargin)
                 Input = find(strcmpi(varargin, 'pDis'));
                 if isempty(Input),pDis = 5*10^-8;else pDis = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'pRep'));
                 if isempty(Input),pRep = 5*10^-8;else pRep = varargin{Input+1};end
                 out.Discovery = sum(obj.pF<=pDis);
                 if ~isempty(obj.ReplicationTest)
                    out.Replication = sum(obj.ReplicationTest.pF<=pRep);
                    out.ReplicationA = sum(obj.ReplicationTest.pA<=pRep);
                 else
                    out.Replication = 0;
                    out.ReplicationA = 0;
                 end
                 if ~isempty(obj.ClassifierTest)
                     tmp = zeros(1,length(obj.ClassifierTest));
                     for i=1:length(obj.ClassifierTest)
                        if isempty(obj.ClassifierTest{i}),continue;end
                        tmp(i) = obj.ClassifierTest{i}.pAUC<=pRep; 
                     end
                     out.Classifier = sum(tmp);
                 else
                     out.Classifier = 0;
                 end
        end
        function out = reportProp(obj,varargin)
                 Input = find(strcmpi(varargin, 'pDis'));
                 if isempty(Input),pDis = 5*10^-8;else pDis = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'pRep'));
                 if isempty(Input),pRep = 5*10^-8;else pRep = varargin{Input+1};end
                 out.Discovery = sum(obj.pF<=pDis);
                 if ~isempty(obj.ReplicationTest)
                    out.Replication = sum(obj.ReplicationTest.pF<=pRep);
                    out.ReplicationA = sum(obj.ReplicationTest.pA<=pRep);
                 else
                    out.Replication = 0;
                    out.ReplicationA = 0;
                 end
                 if ~isempty(obj.ClassifierTest)
                     out.Classifier = obj.ClassifierTest.pAUC<=pRep;
                 else
                     out.Classifier = 0;
                 end
        end
    end
    methods % INTERFACING
       function [out,Beta] = redDepVar(obj,in,COV)
           out = in;
           Beta = in;
           for c=1:1:obj.nModules
               [out{obj.Levels(c)}.DepVar{obj.Clusters(c)},Beta{obj.Levels(c)}.DepVar{obj.Clusters(c)}] = ...
                            PPMMV3.getResiduals(COV,in{obj.Levels(c)}.DepVar{obj.Clusters(c)});
           end
       end    
       function out = getTestX(obj,COV,GB,GT,RS,nrT)
           out = nan*zeros(nrT,1);
            switch obj.ID
                case {1 2 3 4 5 6 7 8 9 10}
                    if isempty(COV), return; end
                    out = COV(:,obj.ID);
                case 100
                    if isempty(GB), return; end
                    out = GB(:,obj.GBID);
                case 1000
                    ind = find(strcmp(obj.RS,RS));
                    if isempty(ind),return;end
                    out = GT(:,ind);
                otherwise
                    return;      
            end 
       end
       function out = getPosClassBal(obj,X)
           X(isnan(X)) = [];
           Tind = X==obj.PosClass;
           out = sum(Tind)/length(X);
       end
       function out = imageMatrix(obj,in,plotFrame)
                if nargin<3, plotFrame = true; end
                out = imageMatrix(obj.HI,in,plotFrame);
       end
       function reset(obj)
           obj.pF = [];
           obj.F = [];
           obj.ReplicationTest = [];
           obj.ClassifierTest = [];
           obj.Classifier = [];
           obj.CV= [];
           obj.OffsetPXC = [];
           obj.ScalePXC = [];
           obj.AvgFS = [];
           obj.StdFS = [];
           obj.SELECT = [];
           obj.ACTION = [];
           obj.MODSELECT = [];
       end  
    end
    methods % IMAGING
        function v = imageSELECT(obj,RefScan,label,varargin)
            %if nargin<4,sel = obj.SELECT;end
            scan = clone(RefScan);
            scan.Material = 'Dull';
            values = zeros(1,scan.nrV);
            index = find(obj.SELECT==1);
            for i=length(index):-1:1
               [lev,cl] = Ind2LC(obj.HI,index(i));
               values(label(lev,:)==cl) = -log10(obj.pF(index(i)));
            end
            scan.Value = values;
            scan.ColorMode = 'Indexed';
            v = PPMMV3.setupViewer(scan);
            title(v.RenderAxes,num2str(obj.pF(1)));
            %colormap(v.RenderAxes,'jet');
        end
        function v = imageSELECTHits(obj,RefScan,label,varargin)
            %if nargin<4,sel = obj.SELECT;end
            scan = clone(RefScan);
            scan.Material = 'Dull';
            values = zeros(1,scan.nrV);
            index = find(obj.SELECT==1);
            for i=length(index):-1:1
               [lev,cl] = Ind2LC(obj.HI,index(i));
               values(label(lev,:)==cl) = values(label(lev,:)==cl)+1;
            end
            scan.Value = values;
            scan.ColorMode = 'Indexed';
            v = PPMMV3.setupViewer(scan);
            title(v.RenderAxes,num2str(obj.pF(1)));
            %colormap(v.RenderAxes,'jet');
        end
        function f = imagePATHS(obj)
            for i=1:obj.nPATH
                f = figure;
                ax1 = subplot(1,2,1);
                ax2 = subplot(1,2,2);
                sel = obj.PATH(i,:);
                plotHierValues3Dv2(sel,'range',[0 1],'peer',ax1);
                modindex = getAllModulesFromSelection(obj.HI,sel);
                tmp = zeros(1,obj.HI.nLC);
                tmp(modindex) = 1;
                plotHierValues3Dv2(tmp,'range',[0 1],'peer',ax2);
                drawnow;
            end
        end
        function f = imageRosettes(obj)
                 f = figure;f.Color = [1 1 1];
                 %ax = subplot(3,3,1);
                 plotHierValues3Dv2(obj.R,'range',[0 0.5],'peer',gca,'title','Canonical Correlation');
                 ax = subplot(3,3,2);
                 plotHierValues3Dv2(-log10(obj.pF),'crit',-log10(5*10^-8),'range',[0 max(-log10(obj.pF))],'peer',ax,'title',num2str(min(obj.pF)));
                 ax = subplot(3,3,3);
                 f = figure;f.Color = [1 1 1];
                 plotHierValues3Dv2(obj.ReplicationTest.CCARIP_BD,'range',[-1 1],'peer',gca,'title','RIP Beta Disc.');colormap(gca,'jet')
            
        end
    end
    methods (Abstract = true)
    end
% STATIC METHODS    
    methods (Static = true) % ASSOCIATION ANALYSIS
        function [F,pF,B,DIM,R,DF1,DF2,STATS] = parCCA(X,Y,C,pcaperc,correctY)
           if nargin<5, correctY = true;end
           if ~isempty(C)
            index = intersect(PPMMV3.notNAN(C),PPMMV3.notNAN(X));
            if correctY
                Y = PPMMV3.getResiduals(C(index,:),Y(index,:));
                %disp('Y WAS CORRECTED');
            else
                Y = Y(index,:);
            end
            X = PPMMV3.getResiduals(C(index,:),X(index,:));
            %X = X(index,:);
           else
            index = PPMMV3.notNAN(X);   
            Y = Y(index,:);
            X = X(index,:);
           end
           Y = PPMMV3.pcaY(Y,pcaperc);
           %Y = PPMMV3.checkRankYCCA(Y);
           [~,B,R,~,~,STATS] = canoncorr(X,Y);
           F = STATS.F(end);
           pF = STATS.pF(end);
           DIM = size(Y,2);
           DF1 = STATS.df1;
           DF2 = STATS.df2;
        end
        function [F,pF,B] = parRIP(X,Y,C,correctY,method,K)
            if nargin<4, correctY = true;end
            if nargin<5, method = 'plsr';end
            if nargin<6, K = 1; end
            if ~isempty(C)
               index = intersect(PPMMV3.notNAN(C),PPMMV3.notNAN(X));
               if correctY
                  Y = PPMMV3.getResiduals(C(index,:),Y(index,:));
                  %disp('Y WAS CORRECTED');
               else
                  Y = Y(index,:);
               end
               X = PPMMV3.getResiduals(C(index,:),X(index,:));
               %X = X(index,:);
            else
              index = PPMMV3.notNAN(X);   
              Y = Y(index,:);
              X = X(index,:);
            end
            if K>0
                switch method
                    case 'cca'
                        RIP = PPMMV3.getLKORIP(X,Y,K,'cca');
                    case 'plsr'
                        RIP = PPMMV3.getLKORIP(X,Y,K,'plsr');
                end
                STATS = regstats(RIP',X,'linear','tstat');
                F = STATS.tstat.t(end);
                pF = STATS.tstat.pval(end);
                B = STATS.tstat.beta(end);
            else
                STATS = PPMMV3.getPERMRIP(X,Y,method,10000);
                F = STATS.T;
                pF = STATS.permpT;
                B = STATS.B;
            end
        end
        function Y = checkRankYCCA(Y)
                 [nObs,nV] = size(Y);
                 nVfrac = 1.2;
                 [Q2,T22,perm2] = qr(Y,0);
                 rankY = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(nObs,nV));
%                  if rankY == 0
%                         error(message('stats:canoncorr:BadData', 'Y'));
%                  elseif rankY < nV
%                         disp('ADJUSTING RANK Y');
%                         adjust = true;
%                         counter = 1;
%                         Y = PPMMV3.pcaY(Y,100);
%                         while adjust
%                             counter = counter+1;
%                             [~,nV] = size(Y);
%                             [~,T22,~] = qr(Y,0);
%                             rankY = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(nObs,nV));
%                             if ~(rankY<nV), 
%                                adjust = false;
%                             else
%                                Y = Y(:,1:end-1);
%                             end
%                         end
%                  end
                 if rankY == 0
                        error(message('stats:canoncorr:BadData', 'Y'));
                 elseif rankY < nV
                        disp('ADJUSTING RANK Y');
                        adjust = true;
                        counter = 0;
                        while adjust
                            Y = PPMMV3.pcaY(Y,99.9);
                            counter = counter+1;
                            [~,nV] = size(Y);
                            [Q2,T22,perm2] = qr(Y,0);
                            rankY = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(nObs,nV));
                            if ~(rankY<nV), adjust = false;end
                        end
                 end
        end
        function [F,pF,B] = parMANCOVAN(X,Y,C,pcaperc)
            if ~isempty(C)
               index = intersect(PPMMV3.notNAN(C),PPMMV3.notNAN(X));   
            else
               index = PPMMV3.notNAN(X);   
            end
            X = X(index);
            Y = Y(index,:);
            if ~isempty(C), C = C(index,:); end
            Y = PPMMV3.pcaY(Y,pcaperc);
            Y = PPMMV3.checkRankYMANCOVAN(X,Y,C);
            if nargout < 3
               [F,pF] = mySNPmancovan(Y, X, C);
            else
               [F,pF,tmp] = mySNPmancovan(Y, X, C);
               B = tmp.B(2,:);
            end
        end
        function Y = checkRankYMANCOVAN(X,Y,C)
                 X = mX(X, C, cell(0));
                 if size(Y, 2) <= size(Y, 1) - rank(X) - 2, return;end
                 maxDim = (size(Y, 1) - rank(X) - 2);
                 Y = Y(:,1:maxDim);
        end
        function Y = pcaY(Y,pcaperc)
            if nargin < 2,pcaperc=100;end
            if isempty(pcaperc), return;end
            tmp = plainPCA;
            tmp.Centering = 0;
            getAverage(tmp,Y');
            getModel(tmp,Y');
            if pcaperc==100, Y = tmp.Tcoeff;return;end
            stripPercVar(tmp,pcaperc);
            Y = tmp.Tcoeff;
        end
        function B = parPLSR(X,Y,C)
           if ~isempty(C)
            index = intersect(PPMMV3.notNAN(C),PPMMV3.notNAN(X));   
            Y = PPMMV3.getResiduals(C(index,:),Y(index,:));
            X = PPMMV3.getResiduals(C(index,:),X(index,:));
            %X = X(index,:);
           else
            index = PPMMV3.notNAN(X);   
            Y = Y(index,:);
            X = X(index,:);
           end
           [~,~,~,~,betha] = plsregress(X,Y,1);
           B = betha(end,:);
        end
        function [T,pT,pFastT,avgT] = parRegUnivariate(X,Y,C)
           if ~isempty(C)
                index = intersect(PPMMV2.notNAN(C),PPMMV2.notNAN(X));
                CX = [C(index,:) X(index)];
           else
               index = PPMMV3.notNAN(X);
               CX = X(index);
           end
           Y = Y(index,:);
           nVar = size(Y,2);
           T = zeros(1,nVar);
           pT = zeros(1,nVar);
           for i=1:1:nVar
               STATS = regstats(Y(:,i),CX,'linear','tstat');
               T(i) = STATS.tstat.t(end);
               pT(i) = STATS.tstat.pval(end);
           end
           pFastT = pfast(pT);
           W = -log10(pT);
           avgT = sum(abs(T).*W)/sum(W);
        end
        function [out,M] = getResiduals(X,Y)
                [~,~,~,~,M] = plsregress(X,Y,min(size(X,2),size(Y,2)));
                Y_est = [ones(size(X,1),1) X]*M;
                out = Y-Y_est;
        end
        function out = lookUppRA(A,Dim,Table)
           out = interp1(Table.Angles,Table.pAngles(Dim,:),A);
        end
        function [A,pA] = angleTest(X1,Y1,C1,X2,Y2,C2,t)
            if nargin<7, t = 0;end
            %B1 =  PPMMV3.parPLSR(X1,Y1,C1);
            %B2 = PPMMV3.parPLSR(X2,Y2,C2);
            B1 =  PPMMV3.PLSR(X1,Y1,C1);
            B2 = PPMMV3.PLSR(X2,Y2,C2);
            A = angle(B1',B2');
            if t==0, pA = nan; return; end
            Acount = zeros(1,t);
            N = length(X1);
            for p=1:1:t
                rng(p);
                ind = randperm(N);
                %forB = PPMMV3.parPLSR(X1,Y1(ind,:),C1);
                forB = PPMMV3.PLSR(X1,Y1(ind,:),C1);
                Acount(p) = angle(forB',B2');
                if mod(p,50)
                   tmppA = (length(find(Acount(1:p)>=A))+1)/(p+1); 
                   acc = 10/p;
                   if tmppA>acc, break;end
                end
            end
            pA = (length(find(Acount(1:p)>=A))+1)/(p+1);
        end
        function [F,pF,T,pT] = ripTest(X1,Y1,C1,X2,Y2,C2)
                 if ~isempty(C1),X1 = PPMMV3.getResiduals(C1,X1);end
                 B1 =  PPMMV3.PLSR(X1,Y1,[]);
                 %B1 =  PPMMV3.BRIM(X1,Y1,C1);
                 %[~,~,B1] = PPMMV3.parMANCOVAN(X1,Y1,C1,[]);
                 RIP = PPMMV3.getRIP(Y2,B1);
                 % ANOVA TESTs
                 [F,pF] = PPMMV3.parMANCOVAN(X2,RIP',[],[]);
                 % REGRESSION TEST
                 %[T,pT] = PPMMV3.parRegUnivariate(X2,RIP',C2);
                 if ~isempty(C2), X2 = PPMMV3.getResiduals(C2,X2);end
                 [T,pT] = PPMMV3.parRegUnivariate(X2,RIP',[]); 
        end
        function [F,pF,T,pT] = ripTestv1(X1,Y1,C1,X2,Y2,C2)
                 B1 =  PPMMV3.PLSR(X1,Y1,[]);
                 %B1 =  PPMMV3.BRIM(X1,Y1,C1);
                 %[~,~,B1] = PPMMV3.parMANCOVAN(X1,Y1,C1,[]);
                 RIP = PPMMV3.getRIP(Y2,B1);
                 % ANOVA TEST
                 [F,pF] = PPMMV3.parMANCOVAN(X2,RIP',[],[]);
                 % REGRESSION TEST
                 %[T,pT] = PPMMV3.parRegUnivariate(X2,RIP',C2);
                 [T,pT] = PPMMV3.parRegUnivariate(X2,RIP',[]); 
        end
        function [out] = getRIP(in,M)
            out = dot(in',repmat(M'/norm(M'),1,size(in,1)));
        end
        function out = combinedMVReplicationTest(X1,Y1,C1,X2,Y2,C2,t,pcaperc,angletable)
                [out.CCAF,out.CCApF] =  PPMMV3.parCCA(X2,Y2,C2,pcaperc);
                [out.RIPF,out.RIPpF,out.RIPT,out.RIPpT] = PPMMV3.ripTest(X1,Y1,C1,X2,Y2,C2);
                if isempty(angletable)
                    [out.A,out.pA] = PPMMV3.angleTest(X1,Y1,C1,X2,Y2,C2,t);
                else
                    [out.A,~] = PPMMV3.angleTest(X1,Y1,C1,X2,Y2,C2,0);
                    out.pA = PPMMV3.lookUppRA(out.A,size(Y1,2),angletable);
                end
                
        end
        function [out,ripout] = FullReplicationTest(X1,Y1,C1,X2,Y2,C2,pcaperc,correctY,angletable,t)
            
            index = intersect(PPMMV3.notNAN(C1),PPMMV3.notNAN(X1));
            C1 = C1(index,:);
            X1 = X1(index,:);
            Y1 = Y1(index,:);
            
            index = intersect(PPMMV3.notNAN(C2),PPMMV3.notNAN(X2));
            C2 = C2(index,:);
            X2 = X2(index,:);
            Y2 = Y2(index,:);
            
            if ~isempty(C1)
               if correctY
                  [Y1,YM1] = PPMMV3.getResiduals(C1,Y1);
                  %Y_est = [ones(size(C2,1),1) C2]*YM1;
                  [Y2,YM2] = PPMMV3.getResiduals(C2,Y2);
                  %Y2 = Y2-Y_est;
                  %disp('Y WAS CORRECTED');
               else
                  %Y = Y(index,:);
               end
               [X1,XM1] = PPMMV3.getResiduals(C1,X1);
               [X2,XM2] = PPMMV3.getResiduals(C2,X2);
               %X_est = [ones(size(C2,1),1) C2]*XM1;
               %X2 = X2-X_est;
               %X = X(index,:);
            else
              %index = PPMMV3.notNAN(X);   
              %Y = Y(index,:);
              %X = X(index,:);
            end
            [out.CCA_FD,out.CCA_pFD,CCAB1,~,out.CCARD] =  PPMMV3.parCCA(X1,Y1,[],pcaperc);
            [out.CCA_FR,out.CCA_pFR,CCAB2,~,out.CCARR] =  PPMMV3.parCCA(X2,Y2,[],pcaperc);
            PLSRB1 =  PPMMV3.PLSR(X1,Y1,[])';
            PLSRB2 =  PPMMV3.PLSR(X2,Y2,[])';
            
            PLSRB2 = sign(angle(PLSRB1,PLSRB2))*PLSRB2;
            %X2 = sign(angle(PLSRB1,PLSRB2))*X2;
            
            CCAB1 = sign(angle(PLSRB1,CCAB1))*CCAB1;
            CCAB2 = sign(angle(PLSRB2,CCAB2))*CCAB2;
            
            out.CCA_A = angle(CCAB1,CCAB2);
            out.CCA_pA = PPMMV3.lookUppRA(abs(out.CCA_A),size(Y1,2),angletable);
            
            %RIP = PPMMV3.getRIP(Y1,CCAB1');
%             RIP = PPMMV3.getLKORIP(X1,Y1,1,'cca',CCAB1);
%             ripout.CCARIPD = RIP;
%             STATS = regstats(RIP',X1,'linear','tstat');
%             out.CCARIP_TD = STATS.tstat.t(end);
%             out.CCARIP_pTD = STATS.tstat.pval(end);
%             out.CCARIP_BD = STATS.tstat.beta(end);

%             STATS = PPMMV3.get2FOLDRIP(X1,Y1,'cca',CCAB1);
%             out.CCARIP_TD1 = STATS.FoldT(1);
%             out.CCARIP_pTD1 = STATS.FoldpT(1);
%             out.CCARIP_BD1 = STATS.FoldB(1);
%             out.CCARIP_TD2 = STATS.FoldT(2);
%             out.CCARIP_pTD2 = STATS.FoldpT(2);
%             out.CCARIP_BD2 = STATS.FoldB(2);
%             out.CCARIP_12A = STATS.A;
%             out.CCARIP_12pA = PPMMV3.lookUppRA(out.CCARIP_12A,size(Y1,2),angletable);
            
            
            RIP = PPMMV3.getRIP(Y2,CCAB1');
            ripout.CCARIPR = RIP;
            STATS = regstats(RIP',X2,'linear','tstat');
            out.CCARIP_TR = STATS.tstat.t(end);
            out.CCARIP_pTR = STATS.tstat.pval(end);
            out.CCARIP_BR = abs(STATS.tstat.beta(end));
            out.CCARIP_SE = STATS.tstat.se(end);
            
            
            out.PLSR_A = angle(PLSRB1,PLSRB2);
            out.PLSR_pA = PPMMV3.lookUppRA(out.PLSR_A,size(Y1,2),angletable);
            
%             RIP = PPMMV3.getLKORIP(X1,Y1,1,'plsr',PLSRB1);
%             ripout.PLSRRIPD = RIP;
%             STATS = regstats(RIP',X1,'linear','tstat');
%             out.PLSRRIP_TD = STATS.tstat.t(end);
%             out.PLSRRIP_pTD = STATS.tstat.pval(end);
%             out.PLSRRIP_BD = STATS.tstat.beta(end);

%             STATS = PPMMV3.get2FOLDRIP(X1,Y1,'plsr',PLSRB1);
%             out.PLSRRIP_TD1 = STATS.FoldT(1);
%             out.PLSRRIP_pTD1 = STATS.FoldpT(1);
%             out.PLSRRIP_BD1 = STATS.FoldB(1);
%             out.PLSRRIP_TD2 = STATS.FoldT(2);
%             out.PLSRRIP_pTD2 = STATS.FoldpT(2);
%             out.PLSRRIP_BD2 = STATS.FoldB(2);
%             out.PLSRRIP_12A = STATS.A;
%             out.PLSRRIP_12pA = PPMMV3.lookUppRA(out.PLSRRIP_12A,size(Y1,2),angletable);
            
            RIP = PPMMV3.getRIP(Y2,PLSRB1');
            ripout.PLSRRIPR = RIP;
            STATS = regstats(RIP',X2,'linear','tstat');
            out.PLSRRIP_TR = STATS.tstat.t(end);
            out.PLSRRIP_pTR = STATS.tstat.pval(end);
            out.PLSRRIP_BR = abs(STATS.tstat.beta(end));
            
            
%             if t==0; return; end
%             nS1 = length(X1);
%             Acount = zeros(1,t);
%             T1Count = zeros(1,t);
%             T2Count = zeros(1,t);
%             T1obs = out.PLSRRIP_TD;
%             T2obs = out.PLSRRIP_TR;
%             obsAngle = out.PLSR_A;
%             parfor i=1:1:t
%                 % i=1;
%                 rng(i);
%                 ind = randperm(nS1);
%                 forB1 =  PPMMV3.PLSR(X1,Y1(ind,:),[]);
%                 forAngle = angle(forB1',B2');
%                 Acount(i) = abs(forAngle)>=abs(obsAngle);
%                 forRIP = PPMMV3.getRIP(Y1,forB1);
%                 STATS = regstats(forRIP',X1,'linear','tstat');
%                 T1Count(i) = STATS.tstat.t(end)>=T1obs;
%                 forRIP = PPMMV3.getRIP(Y2,forB1);
%                 STATS = regstats(forRIP',X2,'linear','tstat');
%                 T2Count(i) = STATS.tstat.t(end)>=T2obs;
%             end     
%             out.permPLSRpA = (length(find(Acount))+1)/(t+1);
%             out.permPLSRpTD = (length(find(T1Count))+1)/(t+1);
%             out.permPLSRpTR = (length(find(T2Count))+1)/(t+1);

        end
        function [out,ripout] = FullReplicationTestINTERNAL(X1,Y1,C1,X2,Y2,C2,pcaperc,correctY,angletable,t)
            if ~isempty(C1)
               if correctY
                  [Y1,YM1] = PPMMV3.getResiduals(C1,Y1);
                  Y_est = [ones(size(C2,1),1) C2]*YM1;
                  Y2 = Y2-Y_est;
                  %disp('Y WAS CORRECTED');
               else
                  %Y = Y(index,:);
               end
               [X1,XM1] = PPMMV3.getResiduals(C1,X1);
               X_est = [ones(size(C2,1),1) C2]*XM1;
               X2 = X2-X_est;
               %X = X(index,:);
            else
              %index = PPMMV3.notNAN(X);   
              %Y = Y(index,:);
              %X = X(index,:);
            end
            [out.CCA_FD,out.CCA_pFD,CCAB1,~,out.CCARD] =  PPMMV3.parCCA(X1,Y1,[],pcaperc);
            [out.CCA_FR,out.CCA_pFR,CCAB2,~,out.CCARR] =  PPMMV3.parCCA(X2,Y2,[],pcaperc);
            PLSRB1 =  PPMMV3.PLSR(X1,Y1,[])';
            PLSRB2 =  PPMMV3.PLSR(X2,Y2,[])';
            CCAB1 = sign(angle(PLSRB1,CCAB1))*CCAB1;
            CCAB2 = sign(angle(PLSRB2,CCAB2))*CCAB2;
            out.CCA_A = angle(CCAB1,CCAB2);
            out.CCA_pA = PPMMV3.lookUppRA(abs(out.CCA_A),size(Y1,2),angletable);
            
            %RIP = PPMMV3.getRIP(Y1,CCAB1');
%             RIP = PPMMV3.getLKORIP(X1,Y1,1,'cca',CCAB1);
%             ripout.CCARIPD = RIP;
%             STATS = regstats(RIP',X1,'linear','tstat');
%             out.CCARIP_TD = STATS.tstat.t(end);
%             out.CCARIP_pTD = STATS.tstat.pval(end);
%             out.CCARIP_BD = STATS.tstat.beta(end);

            STATS = PPMMV3.get2FOLDRIP(X1,Y1,'cca',CCAB1);
            out.CCARIP_TD1 = STATS.FoldT(1);
            out.CCARIP_pTD1 = STATS.FoldpT(1);
            out.CCARIP_BD1 = STATS.FoldB(1);
            out.CCARIP_TD2 = STATS.FoldT(2);
            out.CCARIP_pTD2 = STATS.FoldpT(2);
            out.CCARIP_BD2 = STATS.FoldB(2);
            out.CCARIP_12A = STATS.A;
            out.CCARIP_12pA = PPMMV3.lookUppRA(out.CCARIP_12A,size(Y1,2),angletable);
            
            
            RIP = PPMMV3.getRIP(Y2,CCAB1');
            ripout.CCARIPR = RIP;
            STATS = regstats(RIP',X2,'linear','tstat');
            out.CCARIP_TR = STATS.tstat.t(end);
            out.CCARIP_pTR = STATS.tstat.pval(end);
            out.CCARIP_BR = STATS.tstat.beta(end);
            
            
            out.PLSR_A = angle(PLSRB1,PLSRB2);
            out.PLSR_pA = PPMMV3.lookUppRA(out.PLSR_A,size(Y1,2),angletable);
            
%             RIP = PPMMV3.getLKORIP(X1,Y1,1,'plsr',PLSRB1);
%             ripout.PLSRRIPD = RIP;
%             STATS = regstats(RIP',X1,'linear','tstat');
%             out.PLSRRIP_TD = STATS.tstat.t(end);
%             out.PLSRRIP_pTD = STATS.tstat.pval(end);
%             out.PLSRRIP_BD = STATS.tstat.beta(end);

            STATS = PPMMV3.get2FOLDRIP(X1,Y1,'plsr',PLSRB1);
            out.PLSRRIP_TD1 = STATS.FoldT(1);
            out.PLSRRIP_pTD1 = STATS.FoldpT(1);
            out.PLSRRIP_BD1 = STATS.FoldB(1);
            out.PLSRRIP_TD2 = STATS.FoldT(2);
            out.PLSRRIP_pTD2 = STATS.FoldpT(2);
            out.PLSRRIP_BD2 = STATS.FoldB(2);
            out.PLSRRIP_12A = STATS.A;
            out.PLSRRIP_12pA = PPMMV3.lookUppRA(out.PLSRRIP_12A,size(Y1,2),angletable);
            
            RIP = PPMMV3.getRIP(Y2,PLSRB1');
            ripout.PLSRRIPR = RIP;
            STATS = regstats(RIP',X2,'linear','tstat');
            out.PLSRRIP_TR = STATS.tstat.t(end);
            out.PLSRRIP_pTR = STATS.tstat.pval(end);
            out.PLSRRIP_BR = STATS.tstat.beta(end);
            
            
%             if t==0; return; end
%             nS1 = length(X1);
%             Acount = zeros(1,t);
%             T1Count = zeros(1,t);
%             T2Count = zeros(1,t);
%             T1obs = out.PLSRRIP_TD;
%             T2obs = out.PLSRRIP_TR;
%             obsAngle = out.PLSR_A;
%             parfor i=1:1:t
%                 % i=1;
%                 rng(i);
%                 ind = randperm(nS1);
%                 forB1 =  PPMMV3.PLSR(X1,Y1(ind,:),[]);
%                 forAngle = angle(forB1',B2');
%                 Acount(i) = abs(forAngle)>=abs(obsAngle);
%                 forRIP = PPMMV3.getRIP(Y1,forB1);
%                 STATS = regstats(forRIP',X1,'linear','tstat');
%                 T1Count(i) = STATS.tstat.t(end)>=T1obs;
%                 forRIP = PPMMV3.getRIP(Y2,forB1);
%                 STATS = regstats(forRIP',X2,'linear','tstat');
%                 T2Count(i) = STATS.tstat.t(end)>=T2obs;
%             end     
%             out.permPLSRpA = (length(find(Acount))+1)/(t+1);
%             out.permPLSRpTD = (length(find(T1Count))+1)/(t+1);
%             out.permPLSRpTR = (length(find(T2Count))+1)/(t+1);

        end
        function out = FullReplicationTestv2(X1,Y1,C1,X2,Y2,C2,pcaperc,correctY,angletable,t)
            if ~isempty(C1)
               if correctY
                  [Y1,YM1] = PPMMV3.getResiduals(C1,Y1);
                  Y_est = [ones(size(C2,1),1) C2]*YM1;
                  Y2 = Y2-Y_est;
                  %disp('Y WAS CORRECTED');
               else
                  %Y = Y(index,:);
               end
               [X1,XM1] = PPMMV3.getResiduals(C1,X1);
               X_est = [ones(size(C2,1),1) C2]*XM1;
               X2 = X2-X_est;
               %X = X(index,:);
            else
              %index = PPMMV3.notNAN(X);   
              %Y = Y(index,:);
              %X = X(index,:);
            end
            [out.CCA_FD,out.CCA_pFD,CCAB1] =  PPMMV3.parCCA(X1,Y1,[],pcaperc);
            [out.CCA_FR,out.CCA_pFR,CCAB2] =  PPMMV3.parCCA(X2,Y2,[],pcaperc);
            out.CCA_A = angle(CCAB1,CCAB2);
            out.CCA_pA = PPMMV3.lookUppRA(abs(out.CCA_A),size(Y1,2),angletable);
            
            RIP = PPMMV3.getRIP(Y1,CCAB1');
            STATS = regstats(RIP',X1,'linear','tstat');
            out.CCARIP_TD = STATS.tstat.t(end);
            out.CCARIP_pTD = STATS.tstat.pval(end);
            out.CCARIP_BD = STATS.tstat.beta(end);
            
            RIP = PPMMV3.getRIP(Y2,CCAB1');
            STATS = regstats(RIP',X2,'linear','tstat');
            out.CCARIP_TR = STATS.tstat.t(end);
            out.CCARIP_pTR = STATS.tstat.pval(end);
            out.CCARIP_BR = STATS.tstat.beta(end);        
            
            B1 =  PPMMV3.PLSR(X1,Y1,[]);
            B2 =  PPMMV3.PLSR(X2,Y2,[]);
            out.PLSR_A = angle(B1',B2');
            out.PLSR_pA = PPMMV3.lookUppRA(out.PLSR_A,size(Y1,2),angletable);
            
            out.PLSRCCA_A = angle(B1',CCAB1);
            
            RIP = PPMMV3.getRIP(Y1,B1);
            STATS = regstats(RIP',X1,'linear','tstat');
            out.PLSRRIP_TD = STATS.tstat.t(end);
            out.PLSRRIP_pTD = STATS.tstat.pval(end);
            out.PLSRRIP_BD = STATS.tstat.beta(end);
            
            RIP = PPMMV3.getRIP(Y2,B1);
            STATS = regstats(RIP',X2,'linear','tstat');
            out.PLSRRIP_TR = STATS.tstat.t(end);
            out.PLSRRIP_pTR = STATS.tstat.pval(end);
            out.PLSRRIP_BR = STATS.tstat.beta(end);
            
            
            if t==0; return; end
            nS1 = length(X1);
            Acount = zeros(1,t);
            T1Count = zeros(1,t);
            T2Count = zeros(1,t);
            T1obs = out.PLSRRIP_TD;
            T2obs = out.PLSRRIP_TR;
            obsAngle = out.PLSR_A;
            parfor i=1:1:t
                % i=1;
                rng(i);
                ind = randperm(nS1);
                forB1 =  PPMMV3.PLSR(X1,Y1(ind,:),[]);
                forAngle = angle(forB1',B2');
                Acount(i) = abs(forAngle)>=abs(obsAngle);
                forRIP = PPMMV3.getRIP(Y1,forB1);
                STATS = regstats(forRIP',X1,'linear','tstat');
                T1Count(i) = STATS.tstat.t(end)>=T1obs;
                forRIP = PPMMV3.getRIP(Y2,forB1);
                STATS = regstats(forRIP',X2,'linear','tstat');
                T2Count(i) = STATS.tstat.t(end)>=T2obs;
            end     
            out.permPLSRpA = (length(find(Acount))+1)/(t+1);
            out.permPLSRpTD = (length(find(T1Count))+1)/(t+1);
            out.permPLSRpTR = (length(find(T2Count))+1)/(t+1);

        end
        function B = PLSR(X,Y,C)
            if ~isempty(C)
               index = intersect(PPMMV3.notNAN(C),PPMMV3.notNAN(X));   
               CX = [C(index,:) X(index)];
            else
               index = PPMMV3.notNAN(X);   
               CX = X(index);
            end
            [~,~,~,~,betha] = plsregress(CX,Y(index,:),1);
            B = betha(end,:); 
        end
        function B = BRIM(X,Y,C)
           if ~isempty(C)
              index = intersect(PPMMV3.notNAN(C),PPMMV3.notNAN(X));   
              CX = [C(index,:) X(index)];
           else
              index = PPMMV3.notNAN(X);   
              CX = X(index);
           end
           K = 10;
           nrBoot = 3;
           nrSamples = length(index);
           Y = Y(index,:);
           for fo = 1:1:nrBoot
               Finner = crossvalind('Kfold',nrSamples,K);
               TMP = CX(:,end);% Allocate Memory, only last collumn
               for fi=1:1:K
                   % fi=1;
                   FiTestInd = find(Finner==fi);
                   FiTrInd = setdiff(1:nrSamples,FiTestInd);
                   [~,~,~,~,betha] = plsregress(CX(FiTrInd,:),Y(FiTrInd,:));
                   TMP(FiTestInd) = PPMMV3.getRIP(Y(FiTestInd,:),betha(end,:));
               end
               CX(:,end) = TMP;
           end
           [~,~,~,~,betha] = plsregress(CX,Y,1);
           B = betha(end,:);
        end
        function out = getLKORIP(X,Y,K,method,REFB)
                 nS = length(X);
                 K = round(nS/K);
                 rng(K);
                 F = crossvalind('Kfold',nS,K);
                 out = zeros(1,nS);
                 allind = 1:nS;
                 if nargin<5 % establishing reference direction to cope with sign flip
                    switch method
                        case 'cca'
                           [~,B] = canoncorr(X,Y);
                           [~,~,~,~,B1] = plsregress(X,Y);
                           B1 = B1(end,:)';
                           REFB = sign(angle(B1,B))*B;
                        case 'plsr'
                           [~,~,~,~,REFB] = plsregress(X,Y);
                           REFB = REFB(end,:)';
                    end
                 end
                 for k=1:1:K
                    TestInd = find(F==k);
                    TrInd = setdiff(allind,TestInd);
                    switch method
                        case 'cca'
                           [~,B] = canoncorr(X(TrInd),Y(TrInd,:));
                        case 'plsr'
                           [~,~,~,~,B] = plsregress(X(TrInd),Y(TrInd,:));
                           B = B(end,:)';
                    end
                    B = sign(angle(REFB,B))*B;
                    out(TestInd) = PPMMV3.getRIP(Y(TestInd,:),B');
                 end
        end
        function [out,STATS] = get3FOLDRIP(X,Y,method,REFB)
                 nS = length(X);
                 rng(3);F = crossvalind('Kfold',nS,3);
                 out = zeros(1,nS);
                 Bs = nan*zeros(size(Y,2),3);
                 if nargin<5 % establishing reference direction to cope with sign flip
                    switch method
                        case 'cca'
                           [~,B] = canoncorr(X,Y);
                           [~,~,~,~,B1] = plsregress(X,Y);
                           B1 = B1(end,:)';
                           REFB = sign(angle(B1,B))*B;
                        case 'plsr'
                           [~,~,~,~,REFB] = plsregress(X,Y);
                           REFB = REFB(end,:)';
                    end
                 end
                 orders = [1 2 3;3 1 2;2 3 1];
                 T = nan*zeros(1,3);
                 pT = nan*zeros(1,3);
                 for k=1:1:3
                    TestInd = find(F==orders(k,1));
                    TrInd = find(F==orders(k,2));
                    switch method
                        case 'cca'
                           [~,B] = canoncorr(X(TrInd),Y(TrInd,:));
                        case 'plsr'
                           [~,~,~,~,B] = plsregress(X(TrInd),Y(TrInd,:));
                           B = B(end,:)';
                    end
                    B = sign(angle(REFB,B))*B;
                    out(TestInd) = PPMMV3.getRIP(Y(TestInd,:),B');
                    Bs(:,k) = B;
                    tmp = regstats(out(TestInd)',X(TestInd),'linear','tstat');
                    T(k) = tmp.tstat.t(end);
                    pT(k) = tmp.tstat.pval(end);
                 end
                 tmp = regstats(out',X,'linear','tstat');
                 STATS.T = tmp.tstat.t(end);
                 STATS.pT = tmp.tstat.pval(end);
                 STATS.B = tmp.tstat.beta(end);
                 STATS.pfastT = pfast(pT);
                 STATS.A = -1*(squareform(pdist(Bs','cosine'))-1);
                 STATS.FoldT = T;
                 STATS.FoldpT = pT;
        end
        function STATS = get2FOLDRIP(X,Y,method,REFB)
                 nS = length(X);
                 K = 2;
                 rng(K);F = crossvalind('Kfold',nS,K);
                 Bs = nan*zeros(size(Y,2),K);
                 allind = 1:nS;
                 if nargin<5 % establishing reference direction to cope with sign flip
                    switch method
                        case 'cca'
                           [~,B] = canoncorr(X,Y);
                           [~,~,~,~,B1] = plsregress(X,Y);
                           B1 = B1(end,:)';
                           REFB = sign(angle(B1,B))*B;
                        case 'plsr'
                           [~,~,~,~,REFB] = plsregress(X,Y);
                           REFB = REFB(end,:)';
                    end
                 end
                 T = nan*zeros(1,K);
                 pT = nan*zeros(1,K);
                 Beta = nan*zeros(1,K);
                 for k=1:1:K
                    TestInd = find(F==k);
                    TrInd = setdiff(allind,TestInd);
                    switch method
                        case 'cca'
                           [~,B] = canoncorr(X(TrInd),Y(TrInd,:));
                        case 'plsr'
                           [~,~,~,~,B] = plsregress(X(TrInd),Y(TrInd,:));
                           B = B(end,:)';
                    end
                    B = sign(angle(REFB,B))*B;
                    RIP = PPMMV3.getRIP(Y(TestInd,:),B');
                    Bs(:,k) = B;
                    tmp = regstats(RIP',X(TestInd),'linear','tstat');
                    T(k) = tmp.tstat.t(end);
                    pT(k) = tmp.tstat.pval(end);
                    Beta(k) = tmp.tstat.beta(end);
                 end
                 [~,index] = max(pT);
                 STATS.T = T(index);
                 STATS.pT = pT(index);
                 STATS.B = Beta(index);
                 STATS.pfastT = pfast(pT);
                 STATS.A = angle(Bs(:,1),Bs(:,2));
                 STATS.FoldT = T;
                 STATS.FoldpT = pT;
                 STATS.FoldB = Beta;
        end
        function STATS = getPERMRIP(X,Y,method,t)
                 if nargin < 5, t = 10000;end
                 nS = length(X);
                 switch method
                        case 'cca'
                           [~,B] = canoncorr(X,Y);
                           [~,~,~,~,B1] = plsregress(X,Y);
                           B1 = B1(end,:)';
                           REFB = sign(angle(B1,B))*B;
                        case 'plsr'
                           [~,~,~,~,REFB] = plsregress(X,Y);
                           REFB = REFB(end,:)';
                 end
                 RIP = PPMMV3.getRIP(Y,REFB');
                 tmp = regstats(RIP',X,'linear','tstat');
                 STATS.T = tmp.tstat.t(end);
                 STATS.pT = tmp.tstat.pval(end);
                 STATS.B = tmp.tstat.beta(end);
                 clear B;
                 if t==0, STATS.permpT = STATS.pT;return;end
                 Tcount = zeros(1,t);
                 rng('shuffle');
                 parfor i=1:t
                    permind = randperm(nS);
                    permB = [];
                    switch method
                        case 'cca'
                           [~,permB] = canoncorr(X,Y(permind,:));
                        case 'plsr'
                           [~,~,~,~,permB] = plsregress(X,Y(permind,:));
                           permB = permB(end,:)';
                    end
                    permRIP = PPMMV3.getRIP(Y(permind,:),permB');
                    tmp = regstats(permRIP',X,'linear','tstat'); 
                    Tcount(i) = tmp.tstat.t(end)>=STATS.T; 
                 end
                 STATS.permpT = (length(find(Tcount))+1)/(t+1);            
%                  Tcount = zeros(1,t);
%                  rng('shuffle');
%                  for p=1:t
%                     permind = randperm(nS);
%                     permB = [];
%                     switch method
%                         case 'cca'
%                            [~,permB] = canoncorr(X,Y(permind,:));
%                         case 'plsr'
%                            [~,~,~,~,permB] = plsregress(X,Y(permind,:));
%                            permB = permB(end,:)';
%                     end
%                     permRIP = PPMMV3.getRIP(Y(permind,:),permB');
%                     tmp = regstats(permRIP',X,'linear','tstat'); 
%                     Tcount(p) = tmp.tstat.t(end)>=STATS.T;
%                     if mod(p,50)
%                        tmpp = (sum(Tcount(1:p))+1)/(p+1); 
%                        acc = 2/p;
%                        if tmpp>acc, break;end
%                     end
%                  end
%                 STATS.permpT = (sum(Tcount(1:p))+1)/(p+1);
        end
    end
    methods (Static = true) % INTERFACING
       function out = notNAN(in)
           index = (1:size(in,1));
           [i,~] = find(isnan(in));
           out = setdiff(index,unique(i));
       end
       function out = getL1CL1DepVar(in)
           if ~iscell(in), out = in;return;end
           out = in{1}.DepVar{1};
       end 
       function out = convertUInt(in)
           if isempty(in), out = []; return; end
           M = max(in(:));
           if M<=intmax('uint8'),out = uint8(in);return;end
           if M<=intmax('uint16'), out = uint16(in);return;end
           if M<=intmax('uint32'), out = uint32(in);return;end
           out = uint64(in);
       end
       function out = selectDepVar(DepVar,index)
                 out = cell(size(DepVar));   
                 for i=1:1:length(DepVar);
                     out{i}.DepVar = cell(size(DepVar{i}.DepVar));
                     for j=1:1:length(DepVar{i}.DepVar)
                        tmp = DepVar{i}.DepVar{j};
                        out{i}.DepVar{j} = tmp(index,:); 
                     end
                 end 
       end
       function [out,inavg,instd] = standardize(in)
                [nO,nV] = size(in); 
                inavg = mean(in,1);
                instd = std(in,[],1);
                out = in-repmat(inavg,nO,1);
                out = out./repmat(instd,nO,1);
       end
       function out = CellList2Structure(in)
                out = [];
                if isempty(in), return; end
                nrCells = length(in);
                succes = 0;
                for i=1:1:nrCells
                   if isempty(in{i}),continue;end 
                   names = fieldnames(in{i});
                   for n=1:1:length(names)
                       eval(['out.' names{n} ' = [];']);
                   end
                   succes = 1;
                   break;
                end
                if ~succes, return;end
                for i=1:1:nrCells
                   if isempty(in{i})
                      for n=1:1:length(names)
                        eval(['out.' names{n} ' = [out.' names{n} ' nan];']);
                      end
                   else
                      for n=1:1:length(names)
                        eval(['out.' names{n} ' = [out.' names{n} ' in{i}.' names{n} '];']);
                      end
                   end
                end
       end
    end
    methods (Static = true) % IMAGING
        function v = setupViewer(scan)
            if nargin<1
               v = viewer3DObj;
            else
               v = viewer(scan);
            end
            v.BackgroundColor = [1 1 1];v.AxesVisible = false;v.AxesWallColor = [1 1 1];
            v.AxesXColor = [0 0 0];v.AxesYColor = [0 0 0];v.AxesZColor = [0 0 0];
            v.SceneLightVisible = true;
            v.SceneLightLinked = true;
        end
    end
    methods (Static = true) % GA
        
        
    end
end
