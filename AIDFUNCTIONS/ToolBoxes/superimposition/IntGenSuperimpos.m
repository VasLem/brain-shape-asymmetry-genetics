function [Bs,belief_inlier,weights,sigmas] = IntGenSuperimpos(Bs,kappa)

% intelligentere manier om gezichten op elkaar te plaatsen: per punt wordt
% een verdeling gemaakt, maar er wordt nog geen onderscheid gemaakt tussen
% inliers en outliers.

%% Initialize
change = 1;
n = 0;
nbScans = size(Bs,1);

%% Load initial reference scan & normalise
% disp('initial reference scan...')
ref_scan = randi(nbScans,1);
A = squeeze(Bs(ref_scan,:,:));
nbpts = size(A,2);
k = size(A,1);

B_initial = Bs([1:ref_scan-1,ref_scan+1:nbScans],:,:);

%% initial least squares
disp('initial Least Squares');
parfor i = 1:nbScans-1;
    B = squeeze(B_initial(i,:,:));
    B = [B;ones(1,nbpts)];      % homogeneous coordinates
    % superimpose
    if k==2
        T = LS_2D(A, B);
    else
        T = superimpose_LS_SVD(A, B);
    end
    B = T*B;
    B_initial(i,:,:) = B(1:k,:);
end

Bs(1:ref_scan-1,:,:) = B_initial(1:ref_scan-1,:,:);
Bs(ref_scan,:,:) = A(1:k,:);
Bs(ref_scan+1:end,:,:) = B_initial(ref_scan:end,:,:);

% calculate new mean
A = squeeze(sum(Bs,1)./nbScans);

%% inlier distribution per point
belief_inlier = ones(nbpts,nbScans);
weights = ones(nbpts,nbScans);
% inlier-verdeling
sigmas = zeros(nbpts,1);
% Calculate difference vector E for each face
E = zeros(nbpts, nbScans);
parfor i = 1:nbScans
    B = squeeze(Bs(i,:,:));
    delta = A(1:k,:)-B;
    E2 = sum(delta.^2,1);
    E(:,i) = sqrt(E2);
end
parfor j = 1:nbpts
    [w,sigma] = updateWeights(E(j,:),weights(j,:),kappa);
    belief_inlier(j,:) = w;
    weights(j,:) = w./(2*sigma);
    sigmas(j) = sigma;
end

ObjectiveFct = mean(sum(weights.*E.^2,1));
while change==1
    n = n+1;
    disp(num2str(n));
    ObjectiveFct_old = ObjectiveFct;
    A = [A;ones(1,nbpts)];          % homogeneous coordinates

    % weighted LS
    parfor i = 1:nbScans
        B = squeeze(Bs(i,:,:));
        B = [B;ones(1,nbpts)];      % homogeneous coordinates

        % superimpose
        if k==2
            T = wLS_2D(A, B, weights(:,i)');
        else
            T = superimpose_wLS_SVD(A, B, weights(:,i)');
        end
        B = T*B;
        Bs(i,:,:) = B(1:k,:);
    end
    
    % calculate new mean
    A = squeeze(sum(Bs,1)./nbScans);
    
    % inlier-verdeling
    sigmas = zeros(nbpts,1);
    % Calculate difference vector E for each face
    E = zeros(nbpts, nbScans);
    parfor i = 1:nbScans
        B = squeeze(Bs(i,:,:));
        delta = A(1:k,:)-B;
        E2 = sum(delta.^2,1);
        E(:,i) = sqrt(E2);
    end
    parfor j = 1:nbpts
        [w,sigma] = updateWeights(E(j,:),weights(j,:),kappa);
        belief_inlier(j,:) = w;
        weights(j,:) = w./(2*sigma^2);
        sigmas(j) = sigma;
    end
    
    % change?
    ObjectiveFct = max(sum(weights.*E.^2,1));
    if ObjectiveFct_old-ObjectiveFct <= 10^(-16)
        change = 0;
    end
end
end


function [w,sigma] = updateWeights(E,w,kappa)

% Calculate the new weights for the scans
% input:    E = vector with the differences between the landmarks
%           w = vector with the pervious weights of the landmarks

nbScans = length(E);

% update inlier-destribution parameters
sigma = sqrt(sum(w.*(E.^2),2)/sum(w,2));

% update outlier-destribution parameters
lambda = 1/(sqrt(2*pi)*sigma)*exp(-1/2*kappa^2);

% update weights
w = normpdf(E,0,sigma)./(normpdf(E,0,sigma)+repmat(lambda,1,nbScans));
end