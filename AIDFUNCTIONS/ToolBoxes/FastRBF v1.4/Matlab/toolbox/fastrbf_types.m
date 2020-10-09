% FASTRBF_TYPES Fast 2D and 3D RBF toolbox data types.
%   The FastRBF toolbox is a collection of RBF (radial basis
%   function) tools.  The toolbox uses the following structures.
%
%   All the data structures used are 1-by-1 struct arrays. Square
%   brackets are used to indicate an optional field.
%   Dim is either 2 or 3, but must be constant within a given struct.
%
%     Name              Size
%  --------------------------------------------
%  Attributes:
%     [Value]           1,Np
%     [Accuracy]        1,Np
%     [Lower]           1,Np
%     [Upper]           1,Np
%     [Gradient]        Dim,Np
%     (at least one must be present)
%
%  Point list:
%     Location          Dim,Np
%     [Value]           1,Np
%     [Accuracy]        1,Np
%     [Lower]           1,Np
%     [Upper]           1,Np
%     [Gradient]        Dim,Np
%     (any number of Attribute structs)
%
%  Mesh:
%     (all 3D point list fields)
%     [Tri]             3,Nt
%     [Quad]            4,Nq
%     (one of Tri or Quad must be present)
%
%  Scan:
%     (all 3D point list fields)
%     Size              1,Ns
%     Origin            3,Ns
%
%  Grid:
%     Min               1,Dim
%     Max               1,Dim
%     Spacing           1,Dim
%     Value             X,Y[,Z]
%     [Gradient]        Dim,X,Y[,Z]
%
%  RBF Solution:
%     AchievedAcc       scalar
%     DefaultEvalAcc    scalar
%     Centres           Dim,Np
%     Coeffs            1,Nc
%     PolyBase          1,Dim
%     DataMin           1,Dim
%     DataMax           1,Dim
%     BasicFunc         scalar
%     [BasicFuncParam]  scalar
%     [BasicFuncParam2] scalar
%     PolyDegree        scalar
%     Rho               scalar
%     FitType           scalar
%     Version           string
%     [Indices]         1,Np
%       
%  See also: FASTRBF

% Copyright 2001 Applied Research Associates NZ Ltd
% Author: Tim Evans

helpwin fastrbf_types;
