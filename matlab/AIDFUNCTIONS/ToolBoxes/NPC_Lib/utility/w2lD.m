function [longformat,labels]=w2lD(widedata, StartingColumns,ntimes);
%[longformat]=w2l(widedata,  StartingColumns,ntimes);
%reshape wide dataset to long  format
%
%StartingColumns=colonne in cui si trova il primo tempo di ogni variabile
%ntimes= numero di misure ripetute (uguale per ogni variabile)
%StartingColumns può essere anche una matrice in cui si riportano tutte le
%posizioni di ogni variabile ed ogni tempo:
%StartingColumns=[2 7 12],ntimes=5
%equivale a 
%     StartingColumns=[2     7    12
%                      3     8    13
%                      4     9    14
%                      5    10    15
%                      6    11    16];
%      ntimes=[];
%
%output:
%longformat=dataset in a long format. prima colonna numeri progressivi per
%ogni tempo, poi tutte le cononne non comprese nella trasformazione
%(ripetute per il numero di tempi) a seguire tutte le variabili considerate
%nel reshape.
%labels=cell of strings ={'Time' 'le labels della var' 'per ogni variabile i caratteri comuni
%ai tempi considerati'}
D=widedata;
widedata=[D.vals];
labels=[D.name];

[n ncols]=size(widedata);
if prod(size(StartingColumns)) >length(StartingColumns) %se matrice
    %[ ntimes nvars]=size(StartingColumns);
    temp=StartingColumns;
    clear StartingColumns
else
    all=StartingColumns(1)+(1:(length(StartingColumns)*ntimes))-1;
    temp(1,:)=StartingColumns;
    for j=2:ntimes
    for i=1:length(StartingColumns)
        t=setdiff(setdiff(all, temp(:)),1:StartingColumns(i));
        temp(j,i)=t(1);
    end
    end
end
    [ ntimes nvars]=size(temp);
%     %steps=StartingColumns-StartingColumns(1)+1;
%     all= (1:(length(StartingColumns)*ntimes))
%     %StartingColumns=StartingColumns-StartingColumns(1)+1;
%     StartingColumns
%     for i=1:length(StartingColumns)
%         all=setdiff(all,Sta
%         all(1:ntime)
%         
%         temp(i,:)=all(StartingColumns-StartingColumns(1)+1-(0:i:(length(StartingColumns)-1)));
%         all=setdiff(all,temp(i,:))
%     end
%     
    
%     if length(StartingColumns)==1
%         nvars=length(StartingColumns);
%         step=1;
%     else
%         nvars=length(StartingColumns);
%         step=StartingColumns(2)-StartingColumns(1);
%     end
%     if step==ntimes
%         temp=repmat(StartingColumns,ntimes,1)+repmat((1:ntimes)',1,nvars)-1;
%     elseif step==1
%         temp=reshape(StartingColumns(1)+ (0:(nvars*ntimes-1)),nvars,ntimes)';
%     end
% 
% end
  keepcode=[ setdiff(1:ncols,temp(:)) temp(1,:)];
  
labels2={};

for i=1:nvars
    t=strvcat([labels(temp(:,i))']);
    for ii=2:ntimes
        t=t(1:(end-1),find((t(end,:)-t(end-1,:))==0));
    end
    if isempty(t) | t(1)==' '
        t=['Var' num2str(i)];
    end
    labels2=[labels2 {t}];
end
    
    temp=temp(:);
  
labels2=[{'Time' labels{setdiff(1:ncols,temp)}} labels2 ];

% if exist('code')
% %    clear labels
%     labels3(1).labels='Time';
%     ttt=setdiff(1:ncols,temp);
%     for i=1:length(ttt)
%         labels3(i+1).labels=labels(ttt(i));
%     end
%     
%     
%     for ii=1:size(temp,2)
%         labels3(length(ttt)+1+ii).labels=labels2(length(ttt)+1+ii);
%         labels3(length(ttt)+1+ii).etich=code(temp(i,ii)).etich;
%         for i=2:size(temp,1)
%             [no, id]=setdiff([code(temp(i,ii)).etich{:,2}],[labels3(length(ttt)+1+ii).etich{:,2}]);
%             labels3(length(ttt)+1+ii).etich=[labels3(length(ttt)+1+ii).etich; code(temp(i,ii)).etich(id,:)];
%         end
%         labels3(length(ttt)+1+ii).etich(:,3)=labels3(length(ttt)+1+ii).labels;
%     end
%     labels2=labels3;
% end
labels=labels2;

hold=widedata(:,setdiff(1:ncols,temp));
temp=widedata(:,temp);
temp=reshape(temp,n*ntimes,nvars);
temp=[sort(repmat((1:ntimes)',n,1)) repmat(hold,ntimes,1) temp];


for i=size(temp,2):-1:2
    longformat(i).name=labels(i);
    longformat(i).vals=temp(:,i);
    longformat(i).code=D(keepcode(i-1)).code;
%    longformat(i).type=D(keepcode(i-1)).type;
end
    longformat(1).name=labels(1);
    longformat(1).vals=temp(:,1);
    for i=1:length(unique(longformat(1).vals))
        longformat(1).code(i,:)={['T' num2str(i)], i, 'Time'};
    end
%   longformat(1).type='Nomin';
