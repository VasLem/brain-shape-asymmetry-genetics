function [expanded, orig,w,labels] = expand_categ(mat,var_type,namevar)
%[expanded, orig,w,labels] = expand_categ(MAT,VAR_TYPE,namevar)
%
%It espand a given dataset of categorical variables in dummy variables.
%Given C different numerical values in a single column, expand_categ expand the non
%ordinal variables in C dummy variables. The i-th dummy have value 1 wehn
%the i-th value is present. When the variable is ordinal, expand_categ
%create C-1 variable, with the i-th having value 1 wehn a value greater
%then i-th is present in the original variable.
%It is able to work with missing values.
%
%INPUT:
%MAT = nxm matrix, m number of variables. if mat is a struct array (in
%      NPClib format) also the output will be a struct array.
%VAR_TYPE = specify the type of categorical variable. 1: non ordinal, 2:
%   ordinal. If VAR_TYPE is a scalar, the same type is assumed for each
%   variable, whereas if is a vector of lentgh m, different type are 
%   specified for each variable.
%   LA DUMMYZAZIONE PRODUCE C VARIABILI DUMMY 1 E 0 PER LE VAR CATEG NON
%   ORDNINALI, MENTRE PRE LE ORDINALI FA 0 E 1 PER (I) VS (II, III ecc), (I,
%   II) VS (III,IV ecc)
%namevar = variables names
%const = if const=1 the function will ad a constant ones column 
%for constant original variables.  if const=0 it will omit constant columns. 0 by default.
%
%OUTPUT:
%expanded = matrix of dummy varaibles
%orig = vector having value 1 in the i-th position if the i-th dummy
%   variable come from the first variable, value 2 if it come from the
%   second, etc.
%w = is a vector with the same length of orig. each value is 1/#(dummy var
%   from the same variable).
%labels= is a mx1 cell array containing the labels of each varaibles
%
%Livio Finos
%e-mail: livio@stat.unipd.it


if ischar(mat) 
    namevar={mat};
    [mat]=v(mat);
elseif iscell(mat)
    namevar=mat;
    [mat]=v(mat);
elseif isstruct(mat)
    namevar=[mat.name];
    [mat]=[mat.vals];
end


[n m]=size(mat);
if nargin==1
   var_type=ones(1,m);
end

if length(var_type)==1
    var_type=var_type.*ones(1,m);
end

if not(nargin==3)
   for i=m:-1:1
       namevar{i}=['V' num2str(i)];
   end
end


expanded=zeros(n,0);
orig=zeros(1,0);
w=zeros(1,0);
nlab=0;
for i=1:m
    miss=find(not(isfinite(mat(:,i))));
    notmiss=find((isfinite(mat(:,i))));
    n_cat=unique(mat(notmiss,i));
   if var_type(i)==2
      strg='>=';
      strg=strg((length(n_cat)==1)+1);
      n_cat=n_cat(1:max(1,length(n_cat)-1));
      if not(isempty(n_cat))
          if isstruct(namevar)
              for c=1:length(n_cat)
                  labels{1,nlab+c}=[char(namevar(i).labels) '_' char(namevar(i).etich(c))];
              end 
          else 
              for c=1:length(n_cat)
                  labels{1,nlab+c}=[namevar{i} strg num2str(n_cat(c))];
              end
          end
          nlab=length(labels);
          expa=repmat(NaN,size(mat,1),length(n_cat));
          expa(notmiss,:)=[ repmat(mat(notmiss,i),1,size(n_cat,1))>repmat(n_cat',length(notmiss),1)];
          expanded=[expanded expa];
          orig=[orig i.*ones(1,size(n_cat,1))];
          w=[w ones(1,size(n_cat,1))./size(n_cat,1)];
      end
   else
       if not(isempty(n_cat))
           if isstruct(namevar)
               for c=1:length(n_cat)
                   labels{1,nlab+c}=[char(namevar(i).labels) '=' char(namevar(i).etich(c))];
               end 
           else 
               for c=1:length(n_cat)
                   labels{1,nlab+c}=[namevar{i} '_' num2str(n_cat(c))];
               end
           end
           nlab=length(labels);
           expa=repmat(NaN,size(mat,1),length(n_cat));
           expa(notmiss,:)=[       (repmat(mat(notmiss,i),1,size(n_cat,1))==repmat(n_cat',length(notmiss),1))];
       end
       expanded=[expanded expa];
       orig=[orig i.*ones(1,size(n_cat,1))];
       w=[w ones(1,size(n_cat,1))./size(n_cat,1)];
   end
end
