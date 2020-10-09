function [W] = Rv2Weights(RV,show)
        if nargin < 2,show = true;end
        % PCA of RV
        [P,phi]=eigen(RV);
        G=P*diag(phi.^(1/2));
        tau_RV=round(100 *(phi./sum(phi)));
        % 3.4  Plot the PCA of the Assessors (RV) 
        if show
            figure;
            titre=['PCA of R_V Matrix ',...
                '\lambda_1=',num2str(phi(1)),...
                ' \tau_1=',int2str(tau_RV(1))....
                '\lambda_2=',num2str(phi(2)),...
                ' \tau_2=',int2str(tau_RV(2))];
            %plotxyabs(G,1,2,titre,nom_juge')
            plotxyabs(G,1,2,titre);
        end
        W=P(:,1)/sum(abs(P(:,1)))  ;
end