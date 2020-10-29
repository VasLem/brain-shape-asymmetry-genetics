function T = fastrbf_transformrbf( varargin )

% FASTRBF_TRANSFORMRBF: Applies an affine transformation to an RBF
%    T = FASTRBF_TRANSFORMRBF(S, M) applies the affine transformation  
%    matrix M to the RBF S and returns the transformed RBF T.
%
%    The matrix M must represent an affine (rigid body) transformation. 
%    It should be a square matrix of size 1 larger than the dimension 
%    of the RBF.  Any smaller matrix will be expanded to this size by 
%    adding in the appropriate rows and columns of the identity matrix.  
%
%    The RBF S is transformed so that if S is an interpolatant to data 
%    (x_i, f_i), then the new RBF T is an interpolant to the transformed 
%    data (M x_i, f_i) where x_i is a homogeneous column vector. 
%     
%    See also: FASTRBF, FASTRBF_MAKERBF

T = FastRBF_MEX( 'TransformRBF', varargin{:} );
