function updateCurrentShape(obj,what)
         switch what
             case 'PC'
                 Data = get(obj.Handles.PC.Table,'Data');
                 obj.CurrentCoeff = cell2mat(Data(:,3))';
                 updateShape(obj);
                 obj.CurrentC = CfromY(obj.PredModel,obj.CurrentCoeff);
                 obj.CurrentX = XfromYC(obj.PredModel,obj.CurrentCoeff,obj.CurrentC);
                 updateCTable(obj);  
                 updateCSlider(obj);
                 updateXTable(obj);  
                 updateXSlider(obj);
                 updatePCBar(obj);
             case 'C'
                 Data = get(obj.Handles.C.Table,'Data');
                 oldcoeff = obj.CurrentCoeff;
                 oldC = obj.CurrentC;
                 obj.CurrentC = cell2mat(Data(:,3))';
                 obj.CurrentCoeff = changeCY(obj.PredModel,obj.CurrentC,oldC,oldcoeff);
                 updateShape(obj);
                 updateShapeTable(obj);
                 updateShapeSlider(obj);
                 updatePCBar(obj);
                 obj.CurrentX = XfromYC(obj.PredModel,obj.CurrentCoeff,obj.CurrentC);
                 updateXTable(obj);  
                 updateXSlider(obj);
             case 'X'
                 Data = get(obj.Handles.X.Table,'Data');
                 oldcoeff = obj.CurrentCoeff;
                 oldX = obj.CurrentX;
                 obj.CurrentX = cell2mat(Data(:,3))';
                 obj.CurrentCoeff = changeSingleX(obj.PredModel,obj.CurrentX(obj.ActiveX),oldX(obj.ActiveX),oldcoeff,obj.CurrentC);
                 updateShape(obj);
                 updateShapeTable(obj);
                 updateShapeSlider(obj);
                 updatePCBar(obj);
                 obj.CurrentX = XfromYC(obj.PredModel,obj.CurrentCoeff,obj.CurrentC);
                 updateXTable(obj);  
                 updateXSlider(obj);
             case 'RESET'
                 obj.CurrentCoeff = obj.BaseCoeff;
                 obj.CurrentC = obj.BaseC;
                 obj.CurrentX = obj.BaseX;
                 updateShape(obj);
                 updateShapeTable(obj);
                 updateShapeSlider(obj);
                 updatePCBar(obj);
                 updateCTable(obj);
                 updateCSlider(obj);
                 updateXTable(obj);
                 updateXSlider(obj);
             case 'SET'
                 obj.BaseCoeff = obj.CurrentCoeff;
                 obj.BaseC = obj.CurrentC;
                 obj.BaseX = obj.CurrentX;
              case 'AVERAGE'
                 obj.CurrentCoeff = obj.PredModel.YAvg;
                 obj.CurrentC = CfromY(obj.PredModel,obj.CurrentCoeff);
                 obj.CurrentX = XfromYC(obj.PredModel,obj.CurrentCoeff,obj.CurrentC);
                 updateShape(obj);
                 updateShapeTable(obj);
                 updateShapeSlider(obj);
                 updateCTable(obj);  
                 updateCSlider(obj);
                 updateXTable(obj);  
                 updateXSlider(obj);
                 updatePCBar(obj); 
             case 'IMPORT'
                 [filename, pathname] = uigetfile({'*.xls','Excell file'},'Select Face','MultiSelect','off');
                 if isequal([filename,pathname],[0,0]), return; end
                 cd(pathname);
                 [~,~,raw] = xlsread(filename);
                 obj.CurrentCoeff = cell2mat(raw(2:end,1))';
                 if size(raw,2)==3
                     obj.CurrentC = cell2mat(raw(2:end,2))';
                     obj.CurrentX = cell2mat(raw(2:end,3))';
                 else
                     obj.CurrentC = CfromY(obj.PredModel,obj.CurrentCoeff);
                     obj.CurrentX = XfromYC(obj.PredModel,obj.CurrentCoeff,obj.CurrentC);
                 end
                 updateShape(obj);
                 updateShapeTable(obj);
                 updateShapeSlider(obj);
                 updateCTable(obj);  
                 updateCSlider(obj);
                 updateXTable(obj);  
                 updateXSlider(obj);
                 updatePCBar(obj);    
             otherwise 
                 return;
         end
         %updateLMSelection(obj.ShapeViewer,'Updated Shape');
end

% Data = get(obj.Handles.X.Table,'Data');
%                  oldcoeff = obj.CurrentCoeff;
%                  oldX = obj.CurrentX;
%                  obj.CurrentX = cell2mat(Data(:,3))';
%                  obj.PredModel.ActiveX = obj.ActiveX;
%                  obj.CurrentCoeff = changeX(obj.PredModel,obj.CurrentX,oldX,oldcoeff,obj.CurrentC);
%                  updateShape(obj);
%                  updateShapeTable(obj);
%                  updateShapeSlider(obj);
%                  obj.CurrentX = XfromYC(obj.PredModel,obj.CurrentCoeff,obj.CurrentC);
%                  updateXTable(obj);  
%                  updateXSlider(obj);