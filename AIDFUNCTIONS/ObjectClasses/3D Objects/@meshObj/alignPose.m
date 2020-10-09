function out =  alignPose(obj,refobj,kappa)
    if nargin < 3
       kappa = 0;
    end
    if nargout == 1
       obj = clone(obj);
       out = obj;
    end
    T = rigidTM;
    %T = scaledRigidTM;
    if kappa == 0
        w = ones(3,obj.nrV);
        w(isnan(obj.Vertices))= 0;
        match(T,refobj,obj,w);
        transform(T,obj);
        obj.Value = sum(w)/3;
        return;
    end
    if isnan(kappa)
       C = singleCP('InlierP',gaussianIP,...
                    'OutlierP',dummyOP,...
                    'LatentV',dummyLV,...
                    'Smeasure',FixedPtsDistanceSM);
    else
        C = singleCP('InlierP',gaussianIP,...
                    'OutlierP',uniformOP,...
                    'LatentV',randomBerLV,...
                    'Smeasure',FixedPtsDistanceSM);
    end
        
    % MAP definition
    M = singleMAP('CompleteP',C,...
                  'Tmodel',T);
    M.Floating = obj;%OK
    M.Target = refobj;%OK
    % Optimizer
    objICP = ICP;
    objICP.ChangeTol = 0.01;
    objICP.MaxIter = 400;
    
    solve(objICP,'ObjFun',M,'DA',false,'ML',false,'nrPoints',obj.nrV,'Kappa',kappa,'DefaultOptions');
%     viewer(M.Target);
%     viewer(M.Floating);
%     solveRecord(objICP,'ObjFun',M,'DA',false,'ML',false,'nrPoints',obj.nrV,'Kappa',kappa,'DefaultOptions');
    transform(T,obj);
    obj.Value = M.LatentV.Value;
    try
        addprop(obj,'NoiseLevel');
    catch
        % do nothing
    end
    obj.NoiseLevel = M.InlierP.Sigma;
end

% function out =  alignPose(obj,refobj,kappa)
%     if nargin < 3
%        kappa = 0;
%     end
%     if nargout == 1
%        obj = clone(obj);
%        out = obj;
%     end
%     %T = rigidTM;
%     T = scaledRigidTM;
%     if kappa == 0
%         w = ones(3,obj.nrV);
%         w(isnan(obj.Vertices))= 0;
%         match(T,refobj,obj,w);
%         transform(T,obj);
%         obj.Value = sum(w)/3;
%         return;
%     end
%     C = singleCP('InlierP',gaussianIP,...
%                  'OutlierP',uniformOP,...
%                  'LatentV',randomBerLV,...
%                  'Smeasure',FixedPtsDistanceSM);
%     % MAP definition
%     M = singleMAP('CompleteP',C,...
%                   'Tmodel',T);
%     M.Floating = obj;%OK
%     M.Target = refobj;%OK
%     % Optimizer
%     objICP = ICP;
%     objICP.ChangeTol = 0.01;
%     objICP.MaxIter = 400;
%     
%     solve(objICP,'ObjFun',M,'DA',false,'ML',false,'nrPoints',obj.nrV,'Kappa',kappa,'DefaultOptions');
% 
% %     viewer(M.Target);
% %     viewer(M.Floating);
% %     solveRecord(objICP,'ObjFun',M,'DA',false,'ML',false,'nrPoints',obj.nrV,'Kappa',kappa,'DefaultOptions');
%     transform(T,obj);
%     obj.Value = M.LatentV.Value;
% end