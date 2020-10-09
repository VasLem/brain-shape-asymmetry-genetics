function []=NPC_out(title,options,comb_funct)

maxlengthlabels=12;
maxlengthlabels=floor(maxlengthlabels/2)*2;
options=title;
title=options.title;
comb_funct=options.CF;
 
if options.OUT(1)>0
  if(isfield(options,'labels')==0) options.labels=label_std(options.p.raw,0); end
  if(isfield(options,'Combdims')==0) options.Combdims=length(options.p.raw); end
  
  fprintf('\n________________________________________________________________\n')
fprintf(options.title);
if isfield(options,'B')
fprintf('\nNumber of Conditional Montecarlo: %0.0f',options.B);
end
if isfield(options,'threshold')
fprintf('\nThreshold: %3.4f',options.threshold);
end
    
if not(isempty(options.p.raw)), options.p.raw=split2dims2(options.p.raw,options.Combdims); end
if not(isempty(options.p.adj)), options.p.adj=split2dims2(options.p.adj,options.Combdims); end
if not(isempty(options.p.glb)), options.p.glb=split2dims2(options.p.glb,options.Combdims); end
options.labels.dims=[options.labels.dims(1) options.labels.dims(options.Combdims) options.labels.dims(setdiff(2:length(options.labels.dims),options.Combdims))];
options.labels.dimslabel=[options.labels.dimslabel(1) options.labels.dimslabel(options.Combdims) options.labels.dimslabel(setdiff(2:length(options.labels.dimslabel),options.Combdims))];

fprintf('\n')
for i4=1:size(options.p.raw,4)
    if size(options.p.raw,4)>1
        fprintf('\n\n _____ LEVEL= %-10s _____\n', [char(options.labels.dims{4}{i4})] )
    end
if (length(options.p.raw)+length(options.p.adj))>0
    
tempp=repmat(' ',1,maxlengthlabels);
temp=tempp;
temp2=char(options.labels.dimslabel{2});
temp2=temp2(1:min(maxlengthlabels,length(temp2)));
temp(1:length(temp2))=temp2;
    fprintf([tempp '\t\t\t%-' num2str(maxlengthlabels) 's\n\t%-' num2str(maxlengthlabels) 's'],char(options.labels.dimslabel{3}),temp );
    
    %fprintf('\n%-10s\t',char(options.labels.dimslabel{3}));
    
if isfield(options, 'w')
    fprintf('\tweights\t');
end
 if length(options.labels.dims)>2
    for i=1:size(options.p.raw,3)
        if length(options.labels.dims{3}{i})>=maxlengthlabels
            fprintf(['\t%-' num2str(maxlengthlabels) 's '], [char(options.labels.dims{3}{i}(1:(maxlengthlabels/2-1))) '..' char(options.labels.dims{3}{i}(end-maxlengthlabels/2-1:end))] )
        else
            fprintf(['\t%-' num2str(maxlengthlabels) 's '], char(options.labels.dims{3}(i)))
        end
    end
 fprintf('\n')
 end

for i=1:size(options.p.raw,2)
    
    %if size(options.p.raw,2)>1
        if length(options.labels.dims{2}{i})>=maxlengthlabels
            fprintf(['\t%-' num2str(maxlengthlabels) 's'] , [char(options.labels.dims{2}{i}(1:(maxlengthlabels/2-1))) '..' char(options.labels.dims{2}{i}(end-(maxlengthlabels/2-1):end))])
        else
            fprintf(['\t%-' num2str(maxlengthlabels) 's'], char(options.labels.dims{2}(i)))
        end
    %end

if length(options.p.raw)>0
    if length(options.p.adj)>0 
        fprintf(['\n\t%-' num2str(maxlengthlabels) 's'],'p-value');
    end
    
if isfield(options, 'w')
    fprintf(['\t%-' num2str(maxlengthlabels) '.4G '], options.w(i))
end
    fprintf(['\t%-' num2str(maxlengthlabels) '.5G '], squeeze(options.p.raw(1,i,:,i4)))
    
end
if length(options.p.adj)>0
    fprintf(['\n\t%-' num2str(maxlengthlabels) 's'],'   p-FWE');
    
if isfield(options, 'w')
    fprintf('\t           ', options.w(i))
end
    fprintf(['\t %-' num2str(maxlengthlabels) '.5G'], squeeze(options.p.adj(1,i,:,i4)));
end


fprintf('\n')
end

 
    fprintf('\n Comb Funct\t\t %-6s', (options.CF))
    fprintf('\n p-GLOBAL\t')
    
if isfield(options, 'w')
    fprintf(['\n\t' tempp]);
else
    fprintf(['\n\t' tempp]);
end
 if length(options.labels.dims)>2
    for i=1:size(options.p.raw,3)
        if length(options.labels.dims{3}{i})>=num2str(maxlengthlabels)
        fprintf(['\t%-' num2str(maxlengthlabels) 's'], [char(options.labels.dims{3}{i}(1:(maxlengthlabels/2-1))) '..' char(options.labels.dims{3}{i}(end-(maxlengthlabels/2-1):end))] )
        else
        fprintf(['\t%-' num2str(maxlengthlabels) 's'], char(options.labels.dims{3}(i)))
        end
    end
 end
fprintf(['\n\t' tempp])
    
if isfield(options, 'w')
    fprintf('\t          ')
end
    fprintf(['\t%-' num2str(maxlengthlabels) '.5G '], squeeze(options.p.glb(1,1,:,i4))')

end
end
end

   fprintf('\n')
