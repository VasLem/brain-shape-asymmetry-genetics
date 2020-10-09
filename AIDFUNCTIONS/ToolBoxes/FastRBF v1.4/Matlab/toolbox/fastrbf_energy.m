function E = fastrbf_energy( S, Acc )
% FASTRBF_ENERGY: Compute the energy in a FastRBF solution
%    E = FASTRBF_ENERGY(S) computes the energy for the RBF S.
%
%    E = FASTRBF_ENERGY(S,ACC) uses an evaluation accuracy of ACC
%    rather than the default accuracy from the RBF.
%
%    Note that if two RBFs have different basic functions then the
%    energies given by this function are incomparable.
%
% The Energy in an RBF is given by
%
%    E = scale \sum_{i=1}^{N} \lambda_i s(x_i)
%
% where
%    scale      is a scale factor that depends on the harmonicity,
%    N          is the number of centres,
%    \lambda_i  is the i-th coefficient,
%    x_i        is the i-th centre.
%
% The scale is given by (-1)^m times the factor that relates the basic 
% function to the fundamental solution, E(x), of the appropriate iterated 
% Laplaces equation. These fundamental solutions are:
%
%    thin-plate spline m=2, d=2, E(x) =  1/(8pi)   |x|^2 \log |x|,
%    2-d triharmonic   m=3, d=2, E(x) =  1/(128pi) |x|^4 \log |x|,
%    3-d biharmonic    m=2, d=3, E(x) = -1/(8pi)   |x|,
%    3-d triharmonic   m=3, d=3, E(x) = -1/(96pi)  |x|^3.

% Jon Cherrie
% Applied Research Associates Nz Ltd
% 18 July, 2001

%
% #defines
%
FALSE = 0;
TRUE  = 1;

%
% Parse input arguments
%
if nargin < 1,
   error(['fastrbf_energy: insufficent input arguments.']);
end

if ~isfield( S, 'Centres'     ) | ...
   ~isfield( S, 'BasicFunc'   ) | ...
   ~isfield( S, 'AchievedAcc' ),
   error(['fastrbf_energy: S must be a FastRBF Solution.']);
end

if nargin < 2,
   Acc = max( S.AchievedAcc, 1e-6 );
elseif Acc <= 0,
   error(['fastrbf_energy: Acc must be greater than zero.']);
end

%
% Actual stuff starts here...
%

% find harmonicity
switch (size(S.Centres,1))
   case 2,  % 2D
	switch S.BasicFunc
	case 0, % biharmonic
   	scale = 1/(8*pi);
	case 1,  % triharmonic
   	scale = 1/(128*pi);
	otherwise
   	error(['fastrbf_energy: only biharmonic and triharmonic cases are handled.']);
	end
   case 3, % 3D
	switch S.BasicFunc
	case 0, % biharmonic
   	scale = -1/(8*pi);
	case 1,  % triharmonic
   	scale = 1/(96*pi);
	otherwise
   	error(['fastrbf_energy: only biharmonic and triharmonic cases are handled.']);
	end
end
   

% Eval RBF at all centres
tmp.Location = S.Centres;
tmp = fastrbf_pointeval( S, tmp, Acc, 'quiet' );

% compute energy
N = size( S.Centres, 2 );
E = scale * S.Coeffs(1:N) * tmp.Value';

% EOF
