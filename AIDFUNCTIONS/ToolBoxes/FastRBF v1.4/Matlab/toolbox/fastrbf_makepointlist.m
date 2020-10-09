function pl = fastrbf_makepointlist(loc, val, acc, low, upp, grad)

% FastRBF_MakePointList: create a pointlist struct from existing data
%
%  PL = FASTRBF_MAKEPOINTLIST(LOC) creates a point list
%  struct PL with LOC in the Location field.
%
%  PL = FASTRBF_MAKEPOINTLIST(LOC, VAL, ACC,  LOW,  UPP, GRAD) creates a point list
%  struct PL with VAL, ACC, LOW, UPP, and GRAD in the 'Value', 'Accuracy', 'Lower' errorbar,
%  'Upper' errorbar and 'Gradient' fields respectivly if present and not empty.
%
%  See also: FastRBF, FastRBF_MakeRBF

error(nargchk(1,6, nargin));

s = size(loc);
if length(s) ~= 2
   error('Location argument must be 2 dimensional array');
end

if s(1) ~= 2 & s(1) ~= 3
   if s(2) == 2 | s(2) == 3
      loc = loc';
      s = size(loc);
   else
      error('Location must be dim-by-N where dim is 2 or 3');
   end
end
pl.Location = loc;

if nargin > 1 & ~isempty(val)
   val = val(:)';
   if length(val) ~= s(2)
      error('Value array must be same length as location');
   end
   pl.Value = val;
end

if nargin > 2 & ~isempty(acc)
   acc = acc(:)';
   if length(acc) ~= s(2)
      error('Accuracy array must be same length as location');
   end
   pl.Accuracy = acc;
end

if nargin > 3 & ~isempty(low)
   low = low(:)';
   if length(low) ~= s(2)
      error('Lower error bar array must be same length as location');
   end
   pl.Lower = low;
end

if nargin > 4 & ~isempty(upp)
   upp = upp(:)';
   if length(upp) ~= s(2)
      error('Upper error bar  array must be same length as location');
   end
   pl.Upper = upp;
end

if nargin > 5 & ~isempty(grad)
   if any(s ~= size(grad))
      if all( s == size(grad'))
         grad = grad';
      else
         error('Gradient array must be same size as location');
      end
   end
   pl.Gradient = grad;
end
