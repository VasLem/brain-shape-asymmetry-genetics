function ProgressUpdate( Value, Quiet )
%PROGRESSUPDATE update fastrbf progress bar
%   PROGRESSUPDATE(V) updates the current progress bar to the value V.
%   V should be between 0 and 1.
%
%   See also PROGRESSCREATE, PROGRESSCLOSE

% Copyright 2001 Applied Research Associates NZ Ltd

if nargin == 1
  Quiet = 0;
end

h = findobj(allchild(0), 'flat', 'tag', 'FastRBFProgress');
if isempty(h)
  if Quiet
    return
  else
    error('no progress dialog');
  end
else
  h = h(1);
end

if ischar(Value)
  t = findobj(h, 'tag', 'bartext');
  set(t, 'string', Value, 'visible', 'on');
else
  Value = max(0,min(100*Value,100));

  p = findobj(h, 'tag', 'bar');
  xpatch = [0 Value Value 0];
  set(p, 'xdata', xpatch);

  % if it moved down, we need to call refresh to update properly
  oldvalue = get(p, 'UserData');
  if oldvalue > Value
    refresh(h);
  end
  % store new value
  set(p, 'UserData', Value);
end

% check for abort
if ~Quiet & get(h, 'UserData')
  error('progress dialog aborted')
end
