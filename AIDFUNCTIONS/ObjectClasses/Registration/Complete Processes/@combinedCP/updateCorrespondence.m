function out = updateCorrespondence(obj,Tmodel)
         % weighted combination of CompleteP correspondences, similar to
         % combined Smeasure, with the difference that individual points
         % are weigthed by the LatentV values B
         if isempty(obj.CompleteP), if nargout == 1, out = [];end; return; end
         corresp = fastClone(Tmodel.Evaluation);
         if obj.nrCP == 1
            updateCorrespondence(obj.CompleteP{1},Tmodel);
            corresp = fastCopy(obj.CompleteP{1}.Correspondence);
         else % weighted combination of correspondences in all Complete Processes
            for i=1:1:obj.nrCP % for every Complete Process update correspondence
                updateCorrespondence(obj.CompleteP{i},Tmodel);
            end
            fields = obj.Fields;%execute only once
            % Summing correspondeces per Field over all Smeasures
            % having an affect on that Field;
            for i=1:1:obj.nrFields% for every affected field by correspondence in combined CP
                corresobj = {};w = [];index = 0;
                for j=1:1:obj.nrCP% for every CompleteP in combinedCompleteP
                    if isField(obj.CompleteP{j},fields{i})% has effect on the field
                       index = index+1;
                       corresobj{index} = obj.CompleteP{j}.Correspondence;  %#ok<AGROW>
                       w(index,:) = obj.CompletePW(j)*obj.CompleteP{j}.B; %#ok<AGROW>
                    end
                end
                sumobject = sum(corresobj,fields{i},w);
                corresp.(fields{i}) = sumobject.(fields{i});
                delete(sumobject);
             end
         end
         if nargout == 1, out = corresp; return; end
         obj.Correspondence = corresp;% performs a copy;
         delete(corresp);% hence corresp needs to be deleted           
end