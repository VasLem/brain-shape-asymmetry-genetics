function install(to_whom)

% INSTALL sets the Matlab path to recognize my toolboxes.
% ----------------
% install(to_whom)
% ----------------
% Description: sets the Matlab path to recognize my toolboxes. Should be
%              run upon first installation of the toolbox.
% Input:       <{to_whom}> either 'user' (def) or 'admin'.

% © Liran Carmel
% Last revision date: 18-Feb-2010

% set to_whom
if nargin == 0
    to_whom = 'user';
end

% get current directory
base = pwd;

% set "administrator" path
if strcmp(to_whom,'admin')
    addpath(base);
end

% add general tools
dirname = [base '\General Tools'];
if exist(dirname,'dir')
    addpath(dirname);
end

% add comparative genomic toolbox
dirname = [base '\Comparative Genomics'];
if exist(dirname,'dir')
    addpath(dirname);
end

% add multivariate analysis toolbox
dirname = [base '\Multivariate Analysis'];
if exist(dirname,'dir')
    addpath(dirname);
    addpath([dirname '\Data Handling']);
    addpath([dirname '\Dimensionality Reduction']);
    addpath([dirname '\Graph theory']);
    addpath([dirname '\Information theory']);
    addpath([dirname '\Pairwise Relationships']);
    addpath([dirname '\Statistics']);
end

% add EREM utils toolbox
dirname = [base '\EREM utils'];
if exist(dirname,'dir')
    addpath(dirname);
    addpath([dirname '\Not in Website']);
end

% save everything
savepath;