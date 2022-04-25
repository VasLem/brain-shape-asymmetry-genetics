function out = normEquivalent(obj,model,Kappa,vis,priorvalues,record)
                if nargin < 6, record = false; end
            % building the surfaces
                TS = clone(obj);
                FS = model.Average;
%                 if isempty(FS.PoseLM), indicatePoseLM(FS); end
%                 normalisePose(TS,FS.PoseLM);
            % prior knowledge?
                if isempty(priorvalues)
                   priorvalues = ones(1,FS.nrV);
                end
            % Building similarity measure     
                S = FixedPtsDistanceSM;
                I = gaussianIP;
                O = uniformOP;
                L1 = randomBerLV;
                L2 = detInputLV;
                L2.PriorValues = priorvalues;
                C = singleCP('InlierP',I,'OutlierP',O,'LatentV',combinedLV('LatentV',{L1 L2}),'Smeasure',S);
            % Buidling Tranformation Model
                T = RigidPca;
                T.Fast = false;
                T.Model = clone(model);
                T.Nu = 1;
                T.Pca.Trim = true;
                T.Pca.TrimKappa = 2;
            % Building MAP
                M = singleMAP('CompleteP',C,'Tmodel',T);%
                M.Floating = FS;%OK
                M.Target = TS;%OK
            % extract Normal equivalent PART 3    
                objICP = ICP;
                objICP.ChangeTol = 1;
                objICP.MaxIter = 100;
                DA = DeterministicAnnealing;
                DA.TempFraction = 1e6;
                DA.Rate = 0.6;
                DA.MaxIter = 10;
                DA.FinalZero = true;
                if vis
                   v1 = viewer(TS);
                   v2 = viewer(FS);
                   pos1 = get(v1.Figure,'Position');
                   pos2 = get(v2.Figure,'Position');
                   pos2(1) = pos2(1) + pos1(3);
                   set(v2.Figure,'Position',pos2);
                   pause(1);
                   %FS.Visible = false;
                end
                if ~record
                    solve(objICP,'ObjFun',M,'DA',DA,'ML',false,'nrPoints',FS.nrV,'Kappa',Kappa);
                else
                    solveRecord(objICP,'ObjFun',M,'DA',DA,'ML',false,'nrPoints',FS.nrV,'Kappa',Kappa);
                end
                out = clone(M.Solution);
                T = rigidTM;
                match(T,obj,TS);
                transform(T,out);
                if ~isempty(obj.TextureColor), out.TextureColor = obj.TextureColor; end
                if ~isempty(obj.UV), out.UV = obj.UV; end
                if ~isempty(obj.TextureMap), out.TextureMap = clone(obj.TextureMap); end
                delete(M);
                delete(TS);
                if vis
                    delete(v1);
                    delete(v2);
                end
        end