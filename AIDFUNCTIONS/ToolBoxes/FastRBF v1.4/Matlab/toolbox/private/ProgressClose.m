function ProgressClose
%PROGRESSCLOSE close fastrbf progress bar
%   PROGRESSCLOSE closes the current progress bar.  If no progress bar
%   is found then nothing is done.
%
%   See also PROGRESSCREATE, PROGRESSUPDATE

% Copyright 2001 Applied Research Associates NZ Ltd

h = findobj(allchild(0), 'flat', 'tag', 'FastRBFProgress');
if ~isempty(h)
  h = h(1);
  close(h);
end
