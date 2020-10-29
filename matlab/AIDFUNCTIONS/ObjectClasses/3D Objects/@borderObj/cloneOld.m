function cloneobj = cloneOld(obj)
    cloneobj = borderObj('Visible',false);
    fields = fieldnames(obj);
    for i=1:1:length(fields)
        if strcmp(fields{i},'ph')
               continue;           
        end
        setfield(cloneobj,fields{i},getfield(obj,fields{i}));
    end
    createpatch(cloneobj);
end

