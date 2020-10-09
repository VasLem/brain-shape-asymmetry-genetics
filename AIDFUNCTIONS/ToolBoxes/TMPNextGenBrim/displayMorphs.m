function displayMorphs(morphs,index)
    if nargin < 2, index = 1:size(morphs,2); end
    if size(morphs,1)>1
        nr = length(index);
        for i=1:1:nr
           for j=1:1:2 
               m = morphs{j,index(i)};
               m.SingleColor = [0.8 0.8 0.8];
               m.Material = 'Dull';
               v = viewer(m);
               v. SceneLightVisible = true;
               v.SceneLightLinked = true;
           end
        end
    else
        for j=1:1:2 
               m = morphs{1,j};
               m.SingleColor = [0.8 0.8 0.8];
               m.Material = 'Dull';
               v = viewer(m);
               v. SceneLightVisible = true;
               v.SceneLightLinked = true;
        end
    end
end