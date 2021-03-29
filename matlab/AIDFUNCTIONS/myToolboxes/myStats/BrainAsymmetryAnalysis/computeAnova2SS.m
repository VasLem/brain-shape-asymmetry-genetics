function SS = computeAnova2SS(X,reps)
% Compute Sum of Squares of a 3D matrix of landmarks. Last dimension is
% assumed to correspond to the landmark. This is a direct
% expansion/simplification of the builtin matlab anova2 function

if (nargin < 1)
   error(message('stats:anova2:TooFewInputs'));
end
if (any(isnan(X(:))))
   error(message('stats:anova2:DataNotBalanced'));
end
[r,c,nV] = size(X);
if nargin == 1,
  reps = 1;
  m=r;
  Y = X;
elseif reps == 1
  m=r;
  Y = X;
else
  m = r/reps;
  if (floor(m) ~= r/reps), 
      error(message('stats:anova2:BadSize'));
  end
  Y = zeros(m,c,nV);
  for i=1:m,
      j = (i-1)*reps;
      Y(i,:,:) = mean(X(j+1:j+reps,:,:));
  end
end
colmean = mean(Y,1);        % column means
rowmean = squeeze(mean(Y,2))';       % row means
gm = mean(colmean);         % grand mean
colmean = squeeze(colmean)';
gm = squeeze(gm);
CSS = m*reps*sum((colmean - gm).^2,2);              % Column Sum of Squares
RSS = c*reps*sum((rowmean - gm).^2,2);              % Row Sum of Squares
correction = (c*m*reps)*gm.^2;
TSS =  squeeze(sum(sum(X .^2))) - correction;                    % Total Sum of Squares
ISS = reps*squeeze(sum(sum(Y .^2))) - correction - CSS - RSS;   % Interaction Sum of Squares
if reps == 1,
  SSE = ISS;
else
  SSE = TSS - CSS - RSS - ISS;          % Error Sum of Squares
end
SS = [CSS RSS ISS SSE]';