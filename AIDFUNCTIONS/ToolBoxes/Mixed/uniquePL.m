function pl = uniquePL(pl)
   [pl.Location,index] = unique(pl.Location','rows');
   pl.Location = pl.Location';
   pl.Value = pl.Value(index);
   if isfield(pl,'Gradient'),pl.Gradient = pl.Gradient(:,index);end
   if isfield(pl,'Accuracy'),pl.Accuracy = pl.Accuracy(:,index);end
end