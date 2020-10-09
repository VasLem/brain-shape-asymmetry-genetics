function P = fastrbf_pointeval( varargin )
% FASTRBF_POINTEVAL Evaluate an RBF on a list of points
%    P = FASTRBF_POINTEVAL(S,P) evaluates the solution S at the
%    point locations in P and returns a copy of P containing those
%    evaluations in the Value field.
%
%    P = FASTRBF_POINTEVAL(S,P,ACC) sets the evaluation accuracy
%    to ACC.  The default accuracy is S.AchievedAcc/100.
%
%    FASTRBF_POINTEVAL(...,'gradient') also evaluates the RBF
%    gradient at each point.  The gradients are placed in the
%    Gradient field of P.
%
%    P = FASTRBF_POINTEVAL(...,'attr','Name') returns output in the
%    attribute field P.Name.Value instead of the default P.Value
%    (and P.Name.Gradient instead of P.Gradient).
%
%    FASTRBF_POINTEVAL(...,'smooth', WIDTH) filters the rbf values 
%    with a low pass filter of width WIDTH. 3D Only.
%
%    See also: FASTRBF, FASTRBF_GRIDEVAL
%

[altname, val, grad] = FastRBF_MEX( 'PointEval', varargin{:} );
P = varargin{2};
if isempty(altname)
  P.Value = val;
else
  P = setfield(P, altname, 'Value', val);
end
if ~isempty(grad)
  if isempty(altname)
    P.Gradient = grad;
  else
    P = setfield(P, altname, 'Gradient', grad);
  end
end
