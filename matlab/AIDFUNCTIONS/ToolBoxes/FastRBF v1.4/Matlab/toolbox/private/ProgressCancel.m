function ProgressCancel(Handle)
%PROGRESSCANCEL set fastrbf progress bar to cancelled
%   H = PROGRESSCANCEL(H) sets the cancel flag on the progress bar H.
%   It is normally used only by the cancel button on the dialog.
%
%   See also PROGRESSCREATE, PROGRESSUPDATE, PROGRESSCLOSE

% Copyright 2001 Applied Research Associates NZ Ltd

if ishandle(Handle)
  set(Handle, 'UserData', logical(1));
  c = findobj(Handle, 'tag', 'abort');
  if ~isempty(c)
    set(c, 'style', 'text', 'string', 'Aborting...')
  end
end
