function smoothf = smoothFunctionv2(obj,f,naver,type,varargin)
% smoothFuntion - smooth a function defined on a mesh by averaging
% Smooth a function f on a width of naver vertices
%
% f = smoothFunction(obj,f,naver,type)
%
% naver = getoptions(options, 'niter_averaging', 1);

    if size(f,1)<size(f,2)
        f = f';
    end
    smoothf = f;
    % compute normalized averaging matrix
    W = computeSmoothWeightsv2(obj,type,f);  
    % do averaging to smooth the fieldthick
    for i=1:1:size(f,2)       
        for k=1:naver
            f(:,i) = W*f(:,i);
        end
    end
    Input = find(strcmp(varargin, 'Index'));
    if isempty(Input)
       smoothf = f;
    else
       smoothf(varargin{Input+1},:) = f(varargin{Input+1},:);
    end
end