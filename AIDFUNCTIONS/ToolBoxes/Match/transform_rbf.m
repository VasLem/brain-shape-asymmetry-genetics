function [transformed_rbf] = transform_rbf(rbf,transform)

%eerst de translatie
trans = eye(4,4);
trans(1:3,4) = transform(1:3,4);
%rotatie
rot = transform(1:3,1:3);

%eerst transleren
base  = rbf.PolyBase;
base = [base';1];
newbase = transform*base;
newbase = newbase(1:3,:)';

Centres = rbf.Centres;
Centres = [Centres; ones(1,size(Centres,2))];
Centres = transform * Centres;
Centres = Centres(1:3,:);

polycoeff = rbf.Coeffs(1,end-2:end);
polycoeff = rot * polycoeff';


transformed_rbf = fastrbf_makerbf(Centres, rbf.Coeffs, 0,0,1,newbase);
transformed_rbf.Coeffs(1,end-2:end) = polycoeff';
