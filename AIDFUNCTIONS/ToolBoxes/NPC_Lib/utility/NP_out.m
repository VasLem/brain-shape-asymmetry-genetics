function []=NP_out(title,options)
%[]=NP_out(title,options)
%
%Livio Finos
%e-mail: livio@stat.unipd.it


maxlengthlabels=12;
maxlengthlabels=floor(maxlengthlabels/2)*2;

if options.OUT(1)>0
[no m C dim4]=size(options.p.raw);

  if(isfield(options,'labels')==0) options.labels=label_std(1:m,0); end
fprintf('\n________________________________________________________________\n')      
eval(['fprintf(''\n ' title ' \n'');'])
if isfield(options,'B')
fprintf(' Number of Conditional Montecarlo: %0.0f\n',options.B);
end

for i4=1:dim4
    if(dim4>1)
    fprintf('\n_____%-10s: %-10s_____\n',char(options.labels.dimslabel{4}),char(options.labels.dims{4}{i4}))
    end
tempp=repmat(' ',1,maxlengthlabels);
temp=tempp;
temp2=char(options.labels.dimslabel{2});
temp2=temp2(1:min(maxlengthlabels,length(temp2)));
temp(1:length(temp2))=temp2;
    fprintf([tempp '\t\t\t%-' num2str(maxlengthlabels) 's\n\t%-' num2str(maxlengthlabels) 's'],char(options.labels.dimslabel{3}),temp );

%     temp='          ';
% temp(1:length(char(options.labels.dimslabel{3})))=char(options.labels.dimslabel{3});
%     fprintf(['\t\t\t %-1s\n\t%' num2str(maxlengthlabels) 's\t'],char(options.labels.dimslabel{2}),temp );
for i=1:m
    if length(options.labels.dims{2}{i})>=maxlengthlabels
        fprintf(['\t%-' num2str(maxlengthlabels) 's'],[char(options.labels.dims{2}{i}(1:(maxlengthlabels/2-1))) '..' char(options.labels.dims{2}{i}(end-maxlengthlabels/2-1:end))])
    else
        fprintf(['\t%-' num2str(maxlengthlabels) 's'], char(options.labels.dims{2}(i)))
    end
end

if prod(size(options.p.raw))>0
    fprintf('\n');  
    if length(options.labels.dims)<3 
        for i=1:C
            options.labels.dims{3}(i)={''};
        end
    end
    for i=1:C      
        if length(options.labels.dims{3}{i})>=maxlengthlabels
            disp([sprintf(['\t%+' num2str(maxlengthlabels) 's'], [char(options.labels.dims{3}{i}(1:(maxlengthlabels/2-1))) '..' char(options.labels.dims{3}{i}(end-maxlengthlabels/2-1:end))])...
                sprintf(['\t%-' num2str(maxlengthlabels-2) '.4G  '], options.p.raw(1,:,i,i4))]);
        else
            disp([sprintf(['\t%+' num2str(maxlengthlabels) 's']', char(options.labels.dims{3}{i}))...
                sprintf(['\t%-' num2str(maxlengthlabels-2) '.4G  '], options.p.raw(1,:,i,i4))]);
        end
    end
end
    fprintf('\n');
end
end