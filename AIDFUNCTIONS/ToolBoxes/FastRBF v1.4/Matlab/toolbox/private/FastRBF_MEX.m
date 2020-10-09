function varargout = FastRBF(varargin)

out = cell(1,nargout);

if ~isempty(findstr('(R11.1)', version))
  if nargout == 0
    FastRBF_MEXR11p1(varargin{:});
  else
    [out{:}] = FastRBF_MEXR11p1(varargin{:});
  end
elseif ~isempty(findstr('(R12.1)', version))
  if nargout == 0
    FastRBF_MEXR12p1(varargin{:});
  else
    [out{:}] = FastRBF_MEXR12p1(varargin{:});
  end
elseif ~isempty(findstr('(R13)', version))
  if nargout == 0
    FastRBF_MEXR13(varargin{:});
  else
    [out{:}] = FastRBF_MEXR13(varargin{:});
  end
elseif ~isempty(findstr('(R14)', version))
  if nargout == 0
    FastRBF_MEXR14(varargin{:});
  else
    [out{:}] = FastRBF_MEXR14(varargin{:});
  end
elseif ~isempty(findstr('(R2006a)', version))
  if nargout == 0
    FastRBF_MEXR2006a(varargin{:});
  else
    [out{:}] = FastRBF_MEXR2006a(varargin{:});
  end
elseif ~isempty(findstr('(R2006b)', version))
  if nargout == 0
    FastRBF_MEXR2006b(varargin{:});
  else
    [out{:}] = FastRBF_MEXR2006b(varargin{:});
  end
elseif ~isempty(findstr('(R2007a)', version))
  if nargout == 0
    FastRBF_MEXR2007a(varargin{:});
  else
    [out{:}] = FastRBF_MEXR2007a(varargin{:});
  end
else
 % error(sprintf('Unsupported version of Matlab.\n    The FastRBF toolbox
 % currently supports Matlab release 11.1, 12.1, 13, 14 or 2006a/b.\n    Attemping to use 2006b binary.'))  
  if nargout == 0
    FastRBF_MEXR2007a(varargin{:});
  else
    [out{:}] = FastRBF_MEXR2007a(varargin{:});
  end
end

for n=1:nargout
  varargout{n} = out{n};
end
