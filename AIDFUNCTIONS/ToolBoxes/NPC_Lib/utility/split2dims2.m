function mat2=split2dims2(mat,fromdim)

if fromdim==2
    mat2=mat;
else
dims=size(mat);

 if length(size(mat))<fromdim
     dims=[dims ones(1, fromdim-length(size(mat)))];
 end    
 
predim=dims(1:fromdim-1);
postdim=dims(fromdim+1:end);
altre=prod(dims(setdiff(2:length(dims),fromdim)));

mat=reshape(mat,[prod(predim) dims(fromdim) prod(postdim)]);
mat2=zeros([dims(1) dims(fromdim) altre]);

for i=1:dims(fromdim)
    mat2(:,i,:)=reshape(mat(:,i,:),[dims(1) 1 altre]);
end

mat2=reshape(mat2,[dims(1) dims(fromdim) dims(setdiff(2:length(dims),fromdim))]);
end