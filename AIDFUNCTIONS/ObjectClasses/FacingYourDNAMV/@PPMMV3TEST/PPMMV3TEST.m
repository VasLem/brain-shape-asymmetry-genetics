classdef PPMMV3TEST <superClassLight
% PROPERTIES
    properties % GENERAL INTERFACING
        Mode = [];
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
    end
    properties % ASSOCIATION TESTING
       PCAPerc = 99.9;
       Center = false;
       F;
       pF;
       ReplicationTest;
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
    end
    properties (Dependent = true)
    end
    properties (Hidden = true)
    end
    properties(Hidden = true, Dependent = true)  
    end
    properties(Abstract = true) % ABSTRACT PROPERTIES
       X;
       XREG;
       %Predictor;
    end
% OBJECT METHODS    
    methods % CONSTRUCTOR
       function obj = PPMMV3TEST(varargin)
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
           end
       end
       function out = discoveryUniMulti(obj,DepVar,COV,mode)
                if iscell(DepVar), DepVar = getCellData(obj.HI,DepVar,1);end% FULL FACE PCs
                switch lower(mode)
                    case 'multivariate'
                       [out.F,out.pF] =  PPMMV3TEST.parCCA(obj.XREG,DepVar,COV,[]);% MULTIVARIATE
                    case 'univariate'
                       [out.F,out.pF] = PPMMV3TEST.parRegUnivariate(obj.XREG,DepVar,COV);% UNIVARIATE
                end
                obj.pF = out.pF;obj.F = out.F;
       end
       function out = discoveryModHier(obj,DepVar,COV,mode,varargin)
           Input = find(strcmpi(varargin, 'display'));
           if isempty(Input),display = false;else display = varargin{Input+1};end
           F = nan*zeros(1,obj.nModules);
           pF = nan*zeros(1,obj.nModules);
           B = cell(1,obj.nModules);
           DIM = zeros(1,obj.nModules);
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
              [F(i),pF(i),B{i},DIM(i)] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,pcaperc); %#ok<*PFBNS>
           end
           out.F = F;out.pF = pF;out.B = B;out.DIM = DIM;
           if display
              plotHierValues3Dv2(-log10(out.pF),'crit',-log10(5*10^-8),'range',[0 max([-log10(5*10^-8) -log10(out.pF)])],'title',num2str(min(out.pF)));
           end
           obj.F = F;obj.pF = pF;obj.CCADIM = DIM;
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
               [F(modindex(i)),pF(modindex(i))] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,[]);
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
                           [combF,combpF] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc);
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
                          [tmpF,tmppF] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc); 
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
                       [modF,modpF] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,[]); 
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
                           [tmpF,tmppF] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc);
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
                           [tmpF,tmppF] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc);
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
    end
    methods % REPLICATION ANALYSIS
        function out = runReplicationTest(obj,DepVar,COV,RX,RDepVar,RCOV,varargin)
            switch lower(obj.Mode)
               case 'modular'
                   out = replicationMod(obj,DepVar,COV,RX,RDepVar,RCOV,varargin{:});
               case 'hierarchical'
                   out = replicationHier(obj,DepVar,COV,RX,RDepVar,RCOV,varargin{:});    
               case 'propagate'
                   out = replicationProp(obj,DepVar,COV,RX,RDepVar,RCOV,varargin{:});
               case 'univariate'
                   out = replicationUni(obj,RX,RDepVar,RCOV,varargin{:});    
                case 'multivariate'
                   out = replicationMulti(obj,DepVar,COV,RX,RDepVar,RCOV,varargin{:});
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
                     [F(axindex(i)),pF(axindex(i))] = PPMMV3TEST.parRegUnivariate(Robj.XREG,RDepVar(:,axindex(i)),RCOV);
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
                [out.F,out.pF] =  PPMMV3TEST.parCCA(Robj.XREG,RDepVar,RCOV,[]);
                [out.A,out.pA] = PPMMV3TEST.angleTest(obj.XREG,DepVar,COV,Robj.XREG,RDepVar,RCOV,t);
        end
        function out = replicationMod(obj,DepVar,COV,RX,RDepVar,RCOV,varargin)
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
                   forDepVar = getCellData(obj.HI,DepVar,modindex(i));
                   forRDepVar = getCellData(obj.HI,RDepVar,modindex(i));
                   [tmpF(i),tmppF(i)] =  PPMMV3TEST.parCCA(Robj.XREG,forRDepVar,RCOV,[]); %#ok<*PFBNS>
                   [tmpA(i),tmppA(i)] = PPMMV3TEST.angleTest(obj.XREG,forDepVar,COV,Robj.XREG,forRDepVar,RCOV,t);
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
                   forDepVar = getCellData(obj.HI,DepVar,useindex);
                   forRDepVar = getCellData(obj.HI,RDepVar,useindex);
                   [tmpF(i),tmppF(i)] =  PPMMV3TEST.parCCA(Robj.XREG,forRDepVar,RCOV,pcaperc); %#ok<*PFBNS>
                   [tmpA(i),tmppA(i)] = PPMMV3TEST.angleTest(obj.XREG,forDepVar,COV,Robj.XREG,forRDepVar,RCOV,t);
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
            [out.F,out.pF] =  PPMMV3TEST.parCCA(Robj.XREG,RDepVar,RCOV,obj.PCAPerc);
            [out.A,out.pA] = PPMMV3TEST.angleTest(obj.XREG,DepVar,COV,Robj.XREG,RDepVar,RCOV,t);
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
                 switch combined
                     case true
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
                            PPMMV3TEST.getResiduals(COV,in{obj.Levels(c)}.DepVar{obj.Clusters(c)});
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
            v = PPMMV3TEST.setupViewer(scan);
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
               values(label(lev,:)==cl) = values(label(lev,:)==cl)+1;;
            end
            scan.Value = values;
            scan.ColorMode = 'Indexed';
            v = PPMMV3TEST.setupViewer(scan);
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
    end
    methods (Abstract = true)
    end
% STATIC METHODS    
    methods (Static = true) % ASSOCIATION ANALYSIS
        function [F,pF,B,DIM] = parCCA(X,Y,C,pcaperc)
           if ~isempty(C)
            index = intersect(PPMMV3TEST.notNAN(C),PPMMV3TEST.notNAN(X));   
            Y = PPMMV3TEST.getResiduals(C(index,:),Y(index,:));
            X = PPMMV3TEST.getResiduals(C(index,:),X(index,:));
            %X = X(index,:);
           else
            index = PPMMV3TEST.notNAN(X);   
            Y = Y(index,:);
            X = X(index,:);
           end
           Y = PPMMV3TEST.pcaY(Y,pcaperc);
           Y = PPMMV3TEST.checkRankY(Y);
           [~,B,~,~,~,STATS] = canoncorr(X,Y);
           F = STATS.F(end);
           pF = STATS.pF(end);
           DIM = size(Y,2);
        end
        function Y = checkRankY(Y)
                 [nObs,nV] = size(Y);
                 nVfrac = 1.2;
                 [Q2,T22,perm2] = qr(Y,0);
                 rankY = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(nObs,nV));
                 if rankY == 0
                        error(message('stats:canoncorr:BadData', 'Y'));
                 elseif rankY < nV
                        %disp('ADJUSTING RANK Y');
                        adjust = true;
                        counter = 0;
                        while adjust
                            Y = PPMMV3TEST.pcaY(Y,99.9);
                            counter = counter+1;
                            [~,nV] = size(Y);
                            [Q2,T22,perm2] = qr(Y,0);
                            rankY = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(nObs,nV));
                            if ~(rankY<nV), adjust = false;end
                        end
                 end
        end
        function Y = pcaY(Y,pcaperc)
            if nargin < 2,pcaperc=100;end
            if isempty(pcaperc), return;end
            tmp = plainPCA;
            tmp.Centering = 1;
            getAverage(tmp,Y');
            getModel(tmp,Y');
            if pcaperc==100, Y = tmp.Tcoeff;return;end
            stripPercVar(tmp,pcaperc);
            Y = tmp.Tcoeff;
        end
        function B = parPLSR(X,Y,C)
           if ~isempty(C)
            index = intersect(PPMMV3TEST.notNAN(C),PPMMV3TEST.notNAN(X));   
            Y = PPMMV3TEST.getResiduals(C(index,:),Y(index,:));
            X = PPMMV3TEST.getResiduals(C(index,:),X(index,:));
            %X = X(index,:);
           else
            index = PPMMV3TEST.notNAN(X);   
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
               index = PPMMV3TEST.notNAN(X);
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
            B1 =  PPMMV3TEST.parPLSR(X1,Y1,C1);
            B2 = PPMMV3TEST.parPLSR(X2,Y2,C2);
            A = angle(B1',B2');
            if t==0, pA = nan; return; end
            Acount = zeros(1,t);
            N = length(X1);
            for p=1:1:t
                rng(p);
                ind = randperm(N);
                forB = PPMMV3TEST.parPLSR(X1,Y1(ind,:),C1);
                Acount(p) = angle(forB',B2');
                if mod(p,50)
                   tmppA = (length(find(Acount(1:p)>=A))+1)/(p+1); 
                   acc = 10/p;
                   if tmppA>acc, break;end
                end
            end
            pA = (length(find(Acount(1:p)>=A))+1)/(p+1);
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
end




%        function out = runPropParCCA(obj,DepVar,COV,pStart,display)
%            SELECT = zeros(1,obj.nModules);
%            ACTION = zeros(1,obj.nModules);
%            % POSSIBLE ACTION CODES
%            % 3 = Fully Integrated
%            % -3 = Children only Integrated
%            % 2 = best child Integrated
%            % -2 = worst child Integrated
%            % 1 = best child propagated
%            % -1 = Modular
%            % 0 = lowest level modules
%            actions = -3:1:3;
%            MODSELECT = zeros(obj.nModules,obj.nModules);
%            F = zeros(1,obj.nModules);
%            pF = ones(1,obj.nModules);
%            % START WITH LOWEST LEVEL
%            L = obj.nLev;nC = 2^(L-1);
%            modindex = LC2Ind(obj.HI,L*ones(1,nC),1:nC);
%            for i=1:nC
%                forDepVar = getCellData(obj.HI,DepVar,modindex(i));
%                [F(modindex(i)),pF(modindex(i))] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,[]);
%                SELECT(modindex(i)) = pF(modindex(i)) <= pStart;
%                MODSELECT(modindex(i),modindex(i)) = SELECT(modindex(i));
%                ACTION(modindex(i)) = actions(3);
%            end
%            if display
%                plotrow = 1;
%                f = figure;f.Color = [1 1 1];
%                %f.Position = [1          41        1920        1084];
%                plotaxes = cell(2,6);
%                for pp=1:1:3
%                    plotaxes{1,pp} = subplot(2,3,pp); 
%                end
%                for pp=1:1:3
%                    plotaxes{2,pp} = subplot(2,3,pp+3); 
%                end
%                plotHierValues3Dv2(-log10(pF),'crit',-log10(5*10^-8),'peer',plotaxes{plotrow,1});
%                %title(axpvalue,['LEVEL: ' num2str(L) ' PVALUE']);
%                plotHierValues3Dv2(SELECT,'range',[0 1],'peer',plotaxes{plotrow,2});
%                %title(axsel,['LEVEL: ' num2str(L) ' SELECT']);
%                plotHierValues3Dv2(ACTION,'range',[-3 3],'peer',plotaxes{plotrow,3});
%                %title(axaction,['LEVEL: ' num2str(L) ' SELECT']);
%                drawnow;
%                plotrow= plotrow+1;
%                if plotrow==3, plotrow = 1;end
%            end
%            % WORK YOUR WAY UP
%            for L=obj.nLev-1:-1:1
%                % L = 3;
%                nC = 2^(L-1);
%                modindex = LC2Ind(obj.HI,L*ones(1,nC),1:nC);
%                for i=1:nC
%                    % i=1
%                   % INITIALIZE 
%                        ind = modindex(i);
%                        children = getChildren(obj.HI,ind);
%                        pChildren = pF(children);
%                        [bestp,tmpind] = min(pChildren);
%                        bestChild = children(tmpind);
%                        worstChild = children(setdiff(1:2,tmpind));
%                        childrencombined = false;
%                   % COMBINED CHILDREN ANALYSIS
%                        allchildren = getAllChildren(obj.HI,ind);
%                        index = find(SELECT(allchildren)==1);
%                        if ~isempty(index)
%                           setindex = allchildren(index);
%                           forDepVar = getCellData(obj.HI,DepVar,setindex);
%                           [combF,combpF] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc);
%                           if combpF<=bestp
%                              bestp = combpF;
%                              childrencombined = true;
%                              %disp([num2str(i) ': CHILDREN COMBINED']);
%                           end
%                        end
%                   % HIERARCHICAL ANALYSIS 
%                        if isempty(index)% no children left
%                           tmppF = 1;
%                           tmpF = 0;
%                        else
%                           setindex = [ind allchildren(index)]; 
%                           forDepVar = getCellData(obj.HI,DepVar,setindex);  
%                           [tmpF,tmppF] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc); 
%                        end
%                        if tmppF<=bestp
%                           SELECT(ind) = 1;
%                           F(ind) = tmpF;
%                           pF(ind) = tmppF;
%                           MODSELECT(ind,setindex) = 1;
%                           ACTION(ind) = actions(7);
%                           if display, disp([num2str(i) ': FULLY INTEGRATED']);end
%                           continue; % DONE YOU HAVE IMPROVED THE P-VALUE
%                        end
%                   % MODULAR ANALYSIS 
%                        forDepVar = getCellData(obj.HI,DepVar,ind);
%                        [modF,modpF] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,[]); 
%                        if modpF<=bestp
%                           if modpF<=pStart; 
%                              SELECT(ind) = 1;
%                           else
%                              SELECT(ind) = 0;                          
%                           end
%                           F(ind) = modF;
%                           pF(ind) = modpF;
%                           MODSELECT(ind,ind) = SELECT(ind);
%                           ACTION(ind) = actions(3);
%                           index = getAllChildren(obj.HI,ind);
%                           SELECT(index) = 0;
%                           if display, disp([num2str(i) ': MODULAR']);end
%                           continue;
%                        end
%                   % BRANCH WITH THE BEST CHILD;
%                        allchildren = [bestChild getAllChildren(obj.HI,bestChild)];
%                        index = find(SELECT(allchildren)==1);
%                        setindex = [ind allchildren(index)];
%                        forDepVar = getCellData(obj.HI,DepVar,setindex);  
%                        [tmpF,tmppF] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc);
%                        if tmppF<=bestp
%                           SELECT(ind) = 1;
%                           F(ind) = tmpF;
%                           pF(ind) = tmppF;
%                           MODSELECT(ind,setindex) = 1;
%                           ACTION(ind) = actions(6);
%                           if display, disp([num2str(i) ': BEST CHILD INTEGRATED']);end
%                           % remove other branch from select
%                           allchildren = [worstChild getAllChildren(obj.HI,worstChild)];
%                           SELECT(allchildren) = 0;
%                           continue; % DONE YOU HAVE IMPROVED THE P-VALUE
%                        end
%                    % BRANCH WITH THE BEST CHILD;
%                        allchildren = [worstChild getAllChildren(obj.HI,worstChild)];
%                        index = find(SELECT(allchildren)==1);
%                        setindex = [ind allchildren(index)];
%                        forDepVar = getCellData(obj.HI,DepVar,setindex);  
%                        [tmpF,tmppF] =  PPMMV3TEST.parCCA(obj.XREG,forDepVar,COV,obj.PCAPerc);
%                        if tmppF<=bestp
%                           SELECT(ind) = 1;
%                           F(ind) = tmpF;
%                           pF(ind) = tmppF;
%                           MODSELECT(ind,setindex) = 1;
%                           ACTION(ind) = actions(2);
%                           if display, disp([num2str(i) ': WORST CHILD INTEGRATED']);end
%                           % remove other branch from select
%                           allchildren = [bestChild getAllChildren(obj.HI,bestChild)];
%                           SELECT(allchildren) = 0;
%                           continue; % DONE YOU HAVE IMPROVED THE P-VALUE
%                        end    
%                    % REMOVE MODULE
%                        SELECT(ind) = 0;
%                        if childrencombined
%                           F(ind) = combF;
%                           pF(ind) = combpF;
%                           ACTION(ind) = actions(1);
%                           allchildren = getAllChildren(obj.HI,ind);
%                           MODSELECT(ind,allchildren(SELECT(allchildren)==1)) = 1;
%                           if display, disp([num2str(i) ': CHILDREN INTEGRATED']); end
%                        else
%                           F(ind) = F(bestChild); 
%                           pF(ind) = pF(bestChild);
%                           ACTION(ind) = actions(5);
%                           SELECT(worstChild) = 0;
%                           allchildren = [bestChild getAllChildren(obj.HI,bestChild)];
%                           MODSELECT(ind,allchildren(SELECT(allchildren)==1)) = 1;
%                           if display, disp([num2str(i) ': BEST CHILD PROPAGATED']);end
%                           allchildren = [worstChild getAllChildren(obj.HI,worstChild)];
%                           SELECT(allchildren) = 0;
%                        end
%                end
%                if display
%                    plotHierValues3Dv2(-log10(pF),'crit',-log10(5*10^-8),'peer',plotaxes{plotrow,1});
%                    %title(axpvalue,['LEVEL: ' num2str(L) ' PVALUE']);
%                    plotHierValues3Dv2(SELECT,'range',[0 1],'peer',plotaxes{plotrow,2});
%                    %title(axsel,['LEVEL: ' num2str(L) ' SELECT']);
%                    plotHierValues3Dv2(ACTION,'range',[-3 3],'peer',plotaxes{plotrow,3});
%                    %title(axaction,['LEVEL: ' num2str(L) ' SELECT']);
%                    drawnow;
%                    plotrow= plotrow+1;
%                    if plotrow==3, plotrow = 1;end
%                    pause;
%                end
%            end
%            out.F = F;out.pF = pF;out.SELECT = SELECT;out.MODSELECT = MODSELECT;out.ACTION = ACTION;
%            if nargout==1, return; end
%            obj.F = F;obj.pF = pF;
%            obj.SELECT = SELECT;obj.MODSELECT = MODSELECT;obj.ACTION = ACTION;
%        end
