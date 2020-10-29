function u=phi(r)
%
% Function phi of r
%
% Syntax [u] = phi(r)
%

% Remember if using something like the thin-plate spline in 2D
% u = r^2 * log(r)
% you will need to test for r nonzero before taking the log.

u = r;
