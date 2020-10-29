function [rdata, thresh] = categorize(varargin)

% CATEGORIZE turns continuous data into categorized one.
% -------------------------------------------------------------------
% [rdata boundaries] = categorize(data, no_ranks, p_name, p_val, ...)
% -------------------------------------------------------------------
% Description: turns continuous data into categorized one.
% Input:       {data} continuous data (matrix).
%              {no_ranks} number of categories.
%              <{p_name, p_val}> pairs to control the function. Currently
%                   available:
%                   'mode',mode - supports 'normal' mode (default), and
%                       'safe' mode, when caution is taken for data with
%                       many identical values. Under 'normal' mode, a
%                       warning is issued if there is a need in 'safe'
%                       mode.
%                   'dimension',dim - dimension along which to compute
%                       ranks. By default, {dim} = 1.
%                   'algorithm' - determines how to categorize the data.
%                       This is a structure with a mandatory field 'name'.
%                       The other fields change as a function of the
%                       algorithm. Currently available:
%                           alg.name = 'flat'. The data is categorized
%                               such that each category has has
%                               approximately the same number of samples
%                               (default).
%                           alg.name = 'flat_max'. The categories are set
%                               such that, if maximum is taken across
%                               samples, the maximal values are distributed
%                               evenly.
%                           alg.name = 'power'. The categories have
%                               constantly decreasing number of samples,
%                               obeying a power law. Here, alg.pow is the
%                               power.
% Output:      {rdata} categorized data.
%              {boundaries} between categories. First category includes all
%                   samples that are <=boundaries(1), second category
%                   includes all samples that are >boundaries(1) and <=
%                   boundaries(2), etc.

% © Liran Carmel
% Classification: Data manipulations
% Last revision date: 13-May-2008

% parse input line
[data no_ranks dim alg safe_mode] = parseInput(varargin{:});

% flip data
if dim == 2
    data = data';
end

% find thresholds
thresh = compute_threshold(data,alg,no_ranks,dim);

% rank data
rdata = rank_data(data,thresh);

% test if indeed {no_ranks} is the number of distinct rankings
[isSafeModeRequired missing_rank] = checknoranks(rdata,no_ranks);

% modify threshold if safemode is 'on'
if any(isSafeModeRequired)
    % there are missing ranks in the data
    if safe_mode
        warning('caregorize:SafeMode','Safe mode is on.');
        % correct all columns that require correction
        for col = find(isSafeModeRequired)
            coldata = data(:,col);  % data of problematic column
            colth = thresh(:,col);  % thresholds associated with this column
            % We descriminate between two situations:
            % (1) REPEAT. there is one value (or more) that repeats many times.
            %     In that case, we allocate to it a category of its own, remove
            %     it, and proceed recursively.
            % (2) MALTHRESH. there is a certain threshold value, t(i), such that
            %     no data exists that is >=t(i) and <t(i+1). Every such
            %     occurance, divides the data into two blocks: {data<t(i)} and
            %     {data>=t(i+1)}. If the missing rank in MALTHRESH is the first,
            %     this is a special case, similar to REPEAT.
            coldf = diff(colth);
            if any(~coldf)  % REPEAT
                repeat_val = colth(find(~coldf,1)); % the value that repeats
                left_idx = find(coldata < repeat_val);
                repeat_idx = find(coldata == repeat_val);
                right_idx = find(coldata > repeat_val);
                idx_blocks = {left_idx repeat_idx right_idx};
                data_blocks = {coldata(idx_blocks{1}) ...
                    coldata(idx_blocks{2}) coldata(idx_blocks{3})};
                th_blocks = [0.5*(repeat_val + max(data_blocks{1})) ...
                    0.5*(repeat_val + min(data_blocks{3}))];
                rerankable = [true false true]; % designate by FALSE data blocks
                    % for which the number of ranks should not be re-computed
                    % (because the rank is known to be one there)
                % remove empty blocks
                if isempty(left_idx)
                    data_blocks(1) = [];
                    idx_blocks(1) = [];
                    rerankable(1) = [];
                elseif isempty(right_idx)
                    data_blocks(3) = [];
                    idx_blocks(3) = [];
                    rerankable(3) = [];
                end
            else    % MALTHRESH
                mrank = missing_rank{col};  % ranks that are missing in {rdata}
                if mrank(1) == 1
                    repeat_val = colth(mrank(1)); % the value that repeats
                    repeat_idx = find(coldata == repeat_val);
                    right_idx = find(coldata > repeat_val);
                    idx_blocks = {repeat_idx right_idx};
                    data_blocks = {coldata(idx_blocks{1}) ...
                        coldata(idx_blocks{2})};
                    th_blocks = 0.5*(repeat_val + min(data_blocks{2}));
                    rerankable = [false true];
                else
                    no_mranks = length(mrank);
                    data_blocks = cell(1,no_mranks+1);  % data in each blocks
                    idx_blocks = cell(1,no_mranks+1);   % indices of data in each block
                    th_blocks = colth(mrank);           % thresholds separating blocks
                    idx_blocks{1} = find(coldata<colth(mrank(1)));
                    data_blocks{1} = coldata(idx_blocks{1});
                    idx_blocks{end} = find(coldata>=colth(mrank(end)));
                    data_blocks{end} = coldata(idx_blocks{end});
                    for rr=2:no_mranks
                        idx_blocks{rr} = find(coldata>=colth(mrank(rr-1)) & ...
                            coldata<colth(mrank(rr)));
                        data_blocks{rr} = coldata(idx_blocks{rr});
                    end
                    rerankable = true(1,no_mranks+1);
                    % remove zero-size blocks
                    to_remove = [];
                    for rr = 1:length(data_blocks)
                        if isempty(data_blocks{rr})
                            to_remove = [to_remove rr]; %#ok<AGROW>
                        end
                    end
                    if ~isempty(to_remove)
                        data_blocks(to_remove) = [];
                        idx_blocks(to_remove) = [];
                        rerankable(to_remove) = [];
                        if to_remove(1) == 1
                            to_remove(1) = 2;
                        end
                        th_blocks(to_remove-1) = [];
                    end
                end
            end
            % number of blocks
            no_blocks = length(data_blocks);
            % determine how many ranks should be allocated to each data block
            rank_blocks = allocateranks(data_blocks,no_ranks,rerankable);
            % separately rank each block and combine results
            rdata_blocks = cell(1,no_blocks);
            rthresh_blocks = cell(1,no_blocks);
            for bb = 1:no_blocks
                [rdata_blocks{bb} rthresh_blocks{bb}] = ...
                    categorize(data_blocks{bb},rank_blocks(bb),...
                    'mode','safe','dimension',dim,'algorithm',alg);
            end
            % combine all thresholds and data
            rdata(idx_blocks{1},col) = rdata_blocks{1};
            to_add = rank_blocks(1);
            th_tmp = rthresh_blocks{1};
            for bb = 2:length(data_blocks)
                rdata(idx_blocks{bb},col) = rdata_blocks{bb} + to_add;
                th_tmp = [th_tmp; th_blocks(bb-1); rthresh_blocks{bb}]; %#ok<AGROW>
                to_add = to_add + rank_blocks(bb);
            end
            thresh(:,col) = th_tmp;

%             coldf = diff(colth);
%             colrank = no_ranks;
%             idx = find(~coldf,1);
%             fixedth = [];
%             while ~isempty(idx)
%                 fixedth = [fixedth colth(idx)]; %#ok<AGROW>
%                 coldata = coldata(coldata~=colth(idx));
%                 colrank = colrank - 1;
%                 colth = compute_threshold(coldata,alg,colrank,dim);
%                 coldf = diff(colth);
%                 idx = find(~coldf,1);
%             end
%             % distribute ranking as evenly as possible
%             fixedth = sort(fixedth);
%             catsize = length(coldata) / colrank;
%             fixedth = [min(data(:,col))-1 fixedth max(data(:,col))+1]; %#ok<AGROW>
%             intersize = zeros(1,length(fixedth)-1);
%             for ii = 2:length(fixedth)
%                 intersize(ii-1) = sum(coldata>fixedth(ii-1) & ...
%                     coldata<fixedth(ii));
%             end
%             intersize = intersize / catsize;
%             intersize(intersize>0 & intersize<0.5) = 0.5;
%             intersize = round(intersize);
%             [mx mxidx] = max(intersize);
%             intersize(mxidx(1)) = intersize(mxidx(1)) + ...
%                 colrank - sum(intersize);
%             colth = [];
%             fixedth([1 end]) = fixedth([1 end]) + [1 -1];
%             for ii = 2:length(fixedth)
%                 if fixedth(ii) == fixedth(ii-1)
%                     continue;
%                 end
%                 if intersize(ii-1) == 0
%                     colth = [colth; mean(fixedth(ii-1:ii))]; %#ok<AGROW>
%                     continue;
%                 end
%                 partdata = coldata(coldata>=fixedth(ii-1) & ...
%                     coldata<fixedth(ii));
%                 colth = [colth; compute_threshold(partdata,alg,...
%                     intersize(ii-1),dim)]; %#ok<AGROW>
%                 if ii > 2
%                     colth = [colth; 0.5*(fixedth(ii-1) + min(partdata))]; %#ok<AGROW>
%                 end
%                 if ii < length(fixedth)
%                     colth = [colth; 0.5*(fixedth(ii) + max(partdata))]; %#ok<AGROW>
%                 end
%             end
%             thresh(:,col) = sort(colth);
        end
    else
        warning('categorize:UnRepresentedRanks',...
            'some ranks are not represented. Try using ''mode'',''safe''.');
    end
end

% return to original dimensions
if dim == 2
    rdata = rdata';
end

% #########################################################################
function [data, no_ranks, dim, alg, safe_mode] = parseInput(varargin)

% PARSEINPUT parses input line.
% --------------------------------------------------------
% [data no_ranks dim alg safe_mode] = parseInput(varargin)
% -------------------------------------------------------
% Description: parses the input line.
% Input:       {varargin} original input line.
% Output:      {data} data to categorize.
%              {no_ranks} number of categories.
%              {dim} dimension along which to categorize.
%              {alg} algorithm to use.
%              {safe_mode} mode of operation.

% first argument is always {data}.
data = varargin{1};

% second argument is always {no_ranks}.
no_ranks = varargin{2};

% defaults
dim = 1;
alg = struct('name','flat');
safe_mode = false;

% loop on instructions
for ii = 3:2:nargin
    switch str2keyword(varargin{ii},6)
        case 'mode  '   % instruction: mode
            switch str2keyword(varargin{ii+1},4)
                case 'safe'
                    safe_mode = true;
                case 'norm'
                    safe_mode = false;
                otherwise
                    error('%s: unfamiliar instruction for MODE',...
                        varargin{ii+1});
            end
        case 'dimens'   % instruction: dimension
            dim = varargin{ii+1};
        case 'algori'   % instruction: algorithm
            alg = varargin{ii+1};
            switch str2keyword(alg.name,6)
                case 'flat  '   % flat
                    alg.name = 'flat';
                case 'flat_m'   % flat_max
                    alg.name = 'flat_max';
                case 'power '   % power
                    alg.name = 'power';
            end
        otherwise
            error('%s: undefined instruction',varargin{ii});
    end
end

% modify algorithm 'flat_max'
if strcmp(alg.name,'flat_max')
    alg = struct('name','power','pow',1/size(data,3-dim));
end

% #########################################################################
function thresh = compute_threshold(data,alg,no_ranks,dim)

% COMPUTE_THRESHOLD computes threshold values.
% ----------------------------------------------------
% thresh = compute_threshold(data, alg, no_ranks, dim)
% ----------------------------------------------------
% Description: computes threshold values.
% Input:       {data} continuous data (matrix).
%              {alg} algorithm structure.
%              {no_ranks} number of categories.
%              {dim} dimension along which to compute.
% Output:      {thresh} threshold values.

% take care of a single category
if no_ranks < 2
    thresh = [];
    return;
end

% compute percentiles
switch alg.name
    case 'flat'
        ptiles = linspace(100/no_ranks,100-100/no_ranks,no_ranks-1);
    case 'power'
        ptiles = 100 * ((1:(no_ranks-1))/no_ranks) .^ alg.pow;
end

% compute threshold values
thresh = prctile(data,ptiles,dim);

% correct for binary data
if no_ranks == 2
    for col = 1:size(data,2)
        if any(thresh(col) == minmax(data(:,col)))
            thresh(col) = mean(data(:,col));
        end
    end
end

% flip dimensions
if dim == 2
    thresh = thresh';
end

% #########################################################################
function rdata = rank_data(data,thresh)

% RANK_DATA computes ranked data.
% -------------------------------
% rdata = rank_data(data, thresh)
% -------------------------------
% Description: columnwise computation of ranked data.
% Input:       {data} continuous data (matrix).
%              {thresh} boundary values between categories (see above).
% Output:      {rdata} ranked data.

% fill-in ranks
rdata = ones(size(data));
rdata(isnan(data)) = nan;
for col = 1:size(data,2)
    for rr = 1:size(thresh,1)
        idx = find(data(:,col) >= thresh(rr,col));
        rdata(idx,col) = rdata(idx,col) + 1;
    end
end

% #########################################################################
function [isSafeModeRequired, missing_rank] = checknoranks(rdata,no_ranks)

% CHECKNORANKS checks if the correct number of ranks was produced.
% ----------------------------------------------------------------
% [isSafeModeRequired missing_rank] = checknoranks(rdata,no_ranks)
% ----------------------------------------------------------------
% Description: checks if the correct number of ranks was produced.
% Input:       {rdata} ranked data.
%              {no_ranks} number of desired ranks.
% Output:      {missing_rank} what rank number is missing.

% initialize
no_cols = size(rdata,2);
isSafeModeRequired = false(1,no_cols);
missing_rank = cell(1,no_cols);

% check all data columns
for col = 1:no_cols
    uranks = unique(rdata(:,col));
    missing_rank{col} = allbut(uranks,no_ranks);
    if ~isempty(missing_rank{col})
        isSafeModeRequired(col) = true;
    end
end

% #########################################################################
function rank_blocks = allocateranks(data_blocks,no_ranks,rerankable)

% ALLOCATERANKS tells how many ranks should be given to each data block.
% ------------------------------------------------------------
% rank_blocks = allocateranks(data_blocks,no_ranks,rerankable)
% ------------------------------------------------------------
% Description: tells how many ranks should be given to each data block.
% Input:       {data_blocks} blocked data.
%              {no_ranks} number of desired ranks.
%              {rerankable} ranks are reallocated only to blocks where this
%                   vector is TRUE. If it is FALSE, the number of ranks in this
%                   block is taken as 1.
% Output:      {rank_blocks} how many ranks are allocated to each block.

% initialize
no_blocks = length(data_blocks);
rank_blocks = ones(1,no_blocks);

% remove blocks for which ranks should not be computed
mod_blocks = find(rerankable);
eff_no_blocks = length(mod_blocks);
eff_no_ranks = no_ranks - sum(~rerankable);

% find number of data points in each block
size_blocks = zeros(1,eff_no_blocks);
for bb = 1:eff_no_blocks
    size_blocks(bb) = length(data_blocks{mod_blocks(bb)});
end

% optimal category size
ocs = sum(size_blocks) / eff_no_ranks;

% compute the number of ranks allocated for each block
eff_rank_blocks = size_blocks / ocs;
eff_rank_blocks(eff_rank_blocks>0 & eff_rank_blocks<0.5) = 0.5;
eff_rank_blocks = round(eff_rank_blocks);
[mx mxidx] = max(eff_rank_blocks);
eff_rank_blocks(mxidx(1)) = eff_rank_blocks(mxidx(1)) + ...
    eff_no_ranks - sum(eff_rank_blocks);

% update {rank_blocks}
rank_blocks(mod_blocks) = eff_rank_blocks;