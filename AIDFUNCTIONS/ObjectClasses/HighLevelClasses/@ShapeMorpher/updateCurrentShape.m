function updateCurrentShape(obj,what)
         switch what
             case 'PC'
                 Data = get(obj.Handles.PC.Table,'Data');
                 obj.CurrentCoeff = cell2mat(Data(:,3))';
                 updateShape(obj);
%                  updateShapeSlider(obj);
                 if isempty(obj.PredModel), return; end
                 obj.CurrentPred = XInvfromY(obj.PredModel,obj.CurrentCoeff);
                 updatePredTable(obj);  
                 updatePredSlider(obj);
             case 'PRED'
                 Data = get(obj.Handles.Pred.Table,'Data');
                 oldcoeff = obj.CurrentCoeff;
                 oldpred = obj.CurrentPred;
                 obj.CurrentPred = cell2mat(Data(:,3))';
                 obj.CurrentCoeff = changeX(obj.PredModel,obj.CurrentPred,oldpred,oldcoeff);
                 updateShape(obj);
                 updateShapeTable(obj);
                 updateShapeSlider(obj);
                 %updatePredSlider(obj);
             case 'RESET'
                 obj.CurrentCoeff = obj.BaseCoeff;
                 obj.CurrentPred = obj.BasePred;
                 updateShape(obj);
                 updateShapeTable(obj);
                 updateShapeSlider(obj);
                 updatePredTable(obj);
                 updatePredSlider(obj);
             otherwise 
                 return;
         end
end

