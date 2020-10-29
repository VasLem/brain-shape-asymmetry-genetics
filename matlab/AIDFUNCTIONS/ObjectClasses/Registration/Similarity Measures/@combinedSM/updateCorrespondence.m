function out = updateCorrespondence(obj,Tmodel)
         % weighted combination of Smeasure correspondences
         if isempty(obj.Smeasure), if nargout == 1, out = [];end; return; end
         corresp = fastClone(Tmodel.Evaluation);
         if obj.nrSM == 1
            updateCorrespondence(obj.Smeasure{1},Tmodel);
            corresp = fastCopy(obj.Smeasure{1}.Correspondence);
         else % weighted combination of correspondences in all Smeasures
            for i=1:1:obj.nrSM % for every Smeaure update correspondence
                updateCorrespondence(obj.Smeasure{i},Tmodel);
            end
            fields = obj.Fields;%execute only once
            % Summing correspondeces per Field over all Smeasures
            % having an affect on that Field;
            for i=1:1:obj.nrFields% for every affected field by IP
                corresobj = {};w = [];index = 0;
                for j=1:1:obj.nrSM% for every Smeasure in combined Smeasure
                    if isField(obj.Smeasure{j},fields{i})% has effect on the field
                       index = index+1;
                       corresobj{index} = obj.Smeasure{j}.Correspondence;  %#ok<AGROW>
                       w(index,1) = obj.SmeasureW(j); %#ok<AGROW>
                    end
                end
                sumobject = patchObj.sum(corresobj,fields{i},w);
                corresp.(fields{i}) = sumobject.(fields{i});
                delete(sumobject);
             end
         end
         if nargout == 1, out = corresp; return; end
         obj.Correspondence = corresp;% performs a copy;
         delete(corresp);% hence corresp needs to be deleted           
end