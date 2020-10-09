function ProgressTitle( Value )
%PROGRESSTITLE update fastrbf progress bar title string
%   PROGRESSUPDATE(V) updates the current progress bar title to the
%   string V.  If there is no open progress dialog, does nothing.
%
%   See also PROGRESSCREATE, PROGRESSCLOSE

% Copyright 2001 Applied Research Associates NZ Ltd

h = findobj(allchild(0), 'flat', 'tag', 'FastRBFProgress');
if isempty(h)
  return
end

h = h(1);
t = findobj(h, 'tag', 'text');

set(t, 'string', Value);
