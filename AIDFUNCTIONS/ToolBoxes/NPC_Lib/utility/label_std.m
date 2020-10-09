function[labels]=label_std(kind,stats,stats2)  
%[labels]=label_std(kind,stats,stats2) 
%
%kind =  'partial'    (used by NP_1s, NP_2s*, NP_ANOVA1*),
%       'FWE', 'NPC', 'StOrd', (used by NPC, NPC_FWE, NP_StOrd*),
%       'ReM'        (used by NP_ReM*).
%stats = it is the vector/matrix of observed (partial) p-values
%stats2 = it only used with kind='ReM'; it is the DES matrix.
%
%Livio Finos
%e-mail: livio@stat.unipd.it

if not(ischar(kind)),
    stats=kind;
    kind='partial';
end

switch kind
    case {'Ranks'}
        for i=1:size(stats,2)
  labels.dims{2}{i}=['X' num2str(i)];
        end
        labels.dimslabel{2}={'VARIABLES'};
        for i=1:size(stats,1)
  labels.dims{3}{i}={' '};
        end
        labels.dimslabel{3}={' '};
    case {'partial','FWE','NPC','StOrd','BF','1s'}
        if size(stats,1)>1
  labels.dimslabel{1}={'Random Permutations'}; 
  labels.dims{1}{size(stats,1)}={'Observed'};
  for i=1:(size(stats,1)-1)
      labels.dims{1}{i}=[ 'Rand ' num2str(i)];
  end
  
%   labels.dims{1}(:)=[num2str((1:size(stats,1))')];
        else
  labels.dimslabel{1}={' '};
  labels.dims{1}{1}={' '};
        end

        for i=1:size(stats,2)
  labels.dims{2}{i}=['Y' num2str(i)];
        end
        labels.dimslabel{2}={'VARIABLES Y'};
        if size(stats,3)>1
  labels.dimslabel{3}={'VARIABLES X'}; 
  for i=1:size(stats,3)
      labels.dims{3}{i}=['X' num2str(i)];
  end
        else
  labels.dimslabel{3}={' '};
  labels.dims{3}{1}=[' '];
        end
        if size(stats,4)>1
  labels.dimslabel{4}={'STRATA'}; 
  for i=1:size(stats,4)
      labels.dims{4}{i}=['Strata ' num2str(i)];
  end
        end
        
     case {'ReM','MC'}
       labels.dimslabel{2}={'VARIABLES'};
       labels.dimslabel{3}={'COMPARISONS'};
       if nargin<3
           stats2=DES_ReM_std(size(stats,1)+1);
       elseif isstruct(stats2)
           stats2=stats2.DES;
       end
        for i=1:(size(stats2,2))
  vs=find(stats2(:,i)~=0);
  [no id]=sort(-stats2(vs,i));
  labels.dims{3}{i}=[ num2str(vs(id(no>0))') 'vs' num2str(vs(id(no<0))')  ];
  labels.dims{3}{i}=labels.dims{3}{i}(not(isspace( labels.dims{3}{i})));
        end
        for i=1:size(stats,2)
  labels.dims{2}{i}=['Y' num2str(i)];
        end
     case {'Strata'}
       labels.dimslabel{2}={'VARIABLES'};
       labels.dimslabel{3}={'STRATA'};
        for i=1:size(stats,1)
  labels.dims{3}{i}=['Strata' num2str(i)];
        end
        for i=1:size(stats,2)
  labels.dims{2}{i}=['Y' num2str(i)];
        end

     case {'StochOrd_strata','StochOrd_peack'}
       labels.dimslabel{2}={'VARIABLES'};
       labels.dimslabel{3}={'GLOBAL'};
       labels.dimslabel{4}={'STRATA'};
       labels.dims{3}={' '};
        for i=1:size(stats2,1)
  labels.dims{4}{i}=['Strata' num2str(i)];
        end
        for i=1:size(stats,2)
  labels.dims{2}{i}=['Y' num2str(i)];
        end

end
