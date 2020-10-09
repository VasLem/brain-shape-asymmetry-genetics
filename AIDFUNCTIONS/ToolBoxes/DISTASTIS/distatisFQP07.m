%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Script: distatisFQP07.m
% requires Matlab 7 and up 
% for the plot to work correctly 
% -> intruction convhull
% Apart from that will run on matlab 5 and up
% see below for list of functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Herve Abdi (2007)
% This script is an additional material
%  for
%  Abdi, H., Valentin, D., 
%  Chollet, S., & Chrea C. (2007)
%  Analyzing assessors and products in sorting tasks:
%  DISTATIS, theory and practice.
%  Food Quality and Preference, 18, ***-***.
%  Obtained from Herve Abdi's home page
%  www.utd.edu/~herve e-mail: herve@utdallas.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% Script to read and analyze the
% distatis example for the beers
% HA October 4, 2005.
% Revision HA December 8, 2006.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% need matlab m-files:
%   eigen.m plotxyabs.m corr2mat.m plotxyrange.m
%   plotline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Figure 1 gives Figure 6  of AVCC07
%   Figure 2 gives Figure 7a of AVCC07 (ie., dim1*dim2)
%   Figure 4 gives Figure 8  of AVCC07
%   Figure 3 gives an alternate version 
%     of Figure 8 (not in the paper)
%   The program print eps versions of the figures
%   if the parameter prinfigure=1 (see l3 of progran)
%      default is printfigure=0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The program creates the following matrices
%  L1 to L10     Indicator matrices        (cf. Eq. 16)
%  R1 to R10     Co-occurence matrices     (cf. Eq. 17)
%  D1 to D10     Distance matrices         (cf. Eq. 18)
%  cent          Centering matrix          (cf  Eq.  4)
%  r_S1 to r_S10 Cross-product matrices    (cf. Eq. 19)
%  S1 to S2      Normalized cp matrices    (cf. Eq. 20)
%  mat_RV        RV coefficient matrix (C) (cf. Eq. 21)
%  P             Eigenvectors of C         (cf. Eq.  7)
%  phi           Eigenvalues (vector) of C (cf. Eq.  7)
%  G             Factor scores for C (cf. Figure 6)
%  weigths       Alpha weights             (cf. Eq. 22)
%  compromise    Compromise (S+)           (cf. Eq. 23)
%  pc            Eigenvectors of S+ (V)    (cf. Eq. 12)
%  lc            Eigenvalues (vector) 
%                    of S+ (diag(Lambda)   (cf. Eq. 12)
%  F             Factor scores for the
%                   beers (compromise)     (cf. Eq. 13)
%  Rsup4I        Biplot projection of
%                  supplementary variable (cf Figure 7)
%  RProj         Projection Operator       (cf. Eq. 14)
%  RF1 to RF10   Projections of the assessors
%                 as supplementary scores  (cf. Eq. 14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Herve Abdi. December 08, 2006.
% herve@utdallas.edu
% www.utd.edu/~herve
%



% clear everything
clear
clc
printfigure=1; % print figure
printfigure=0; % do not print figures
%
gen_name='beers_norm'; 
% prefix for the graphs/files


% 1. Get the beers/sorting results
%    Get the descriptions of assessors and beers
%    (cf. Table 1 of Abdi et al. 2007).
%%%%%%%%%% Name of the beers %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nom_de_row{1}='Affligen';
Nom_de_row{2}='Budweiser';
Nom_de_row{3}='Buckler Blonde';
Nom_de_row{4}='Killian';
Nom_de_row{5}='St Landelin';
Nom_de_row{6}='Buckler Highland';
Nom_de_row{7}='Fruit Defendu';
Nom_de_row{8}='EKU28';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% These are the data (cf. Table 1 of Abdi et al. 2007).
le_sort=[...
   1,4,3,4,1,1,2,2,1,3
   4,5,2,5,2,3,1,1,4,3
   3,1,2,3,2,4,3,1,1,2
   4,2,3,3,1,1,1,2,1,4
   1,5,3,5,2,1,1,2,1,3
   2,3,1,1,3,5,4,4,3,1
   1,4,3,4,1,1,2,2,2,4
   5,2,4,2,4,2,5,3,4,5
   ];
%  Here we describe the assessors
%  Col1 Gender 1=M 2=F Col2 age  
%  (cf. Table 1 of Abdi et al. 2007).
des_juges=[...
  2 22;1 23;2 23;2 23;1 26;... 
  1 29;1 35;1 25;2 39;1 28];
%%% Here we describe the beers
%%% Col.1 Alcohol*10 Col.2 Average Hedonic*10 
% (cf. Table 1 of Abdi et al. 2007).
des_bieres=[...
    65 55; 50 25; 12 18; 65 43;...
    65 47; 12 10; 88 62; 110 18];         
% Done data entry 
[nbeers,njuges]=size(le_sort);
%
% Now go for the analysis per se
% 2. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1. Create, for each assessor, the 
%  disjonctive coding matrices L1 to L10
%  distance matrices D1 to D10
%  co-occurrence matrices R1 to R10
 for k=1:njuges;
    la_col=le_sort(:,k);
    ncol=max(la_col);
    Lwork=zeros(nbeers,ncol);
    work=zeros(nbeers,nbeers);
    for j=1:nbeers;
        Lwork(j,la_col(j))=1;
    end
    eval(['L',int2str(k),'=Lwork;']);
    for j1=1:nbeers-1;
        for j2=j1+1:nbeers;
            work(j1,j2)=(la_col(j1)~=la_col(j2));
        end
    end
    work=work+work';
    eval(['D',int2str(k),'=work;']);
    eval(['R',int2str(k),'=1-work;']);
 end
 % 2.2 Create the name of the assessors
 %    First letter is F/M for female/male 
for k=1:njuges;
    if des_juges(k,1)==1;sj='M';else,sj='F';end
    nom_juge{k}=['Judge(',sj,'-',...
        int2str(des_juges(k,2)),')',int2str(k)]; 
end
% 3. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go for distatis as said
%
% 3.1. Transform the D matrices
%      in SCP matrices (MDS transform)
%      Create the total distance matrix
 m=ones(nbeers,1)./nbeers;
 cent=eye(nbeers,nbeers)-ones(nbeers,1)*m';
 Dsum=zeros(nbeers,nbeers);
 for k=1:njuges;
     eval(['D=D',int2str(k),';']);
     eval(['Dsum=Dsum+D',int2str(k),';']);
     % transform D -> SCP 
     r_S=-(1/2)*cent*D*cent;
     [PS,lS]=eigen(r_S);
     % normaliz by first eig
     S=r_S.*(lS(1).^(-1));% normalize
     % S=r_S; % don't normalize
     eval(['S',int2str(k),'=S;']);
     eval(['r_S',int2str(k),'=r_S;']);
 end
 
 
 % 3.2. Compute the RV coefficients
 ncells=nbeers*nbeers;
 Grosse_mat=zeros(ncells,njuges);
 for k=1:njuges;
     eval(['la_c=S',int2str(k),'(:);']);
     Grosse_mat(:,k)=la_c;
 end
 num_RV=Grosse_mat'*Grosse_mat;
la_diag=diag(num_RV);
nd=length(la_diag);
den_RV=(repmat(la_diag',nd,1).*repmat(la_diag,1,nd)).^(1/2);
mat_RV=num_RV./den_RV; 
% 3.3.   PCA of the RV matrix
[P,phi]=eigen(mat_RV);
G=P*diag(phi.^(1/2));
tau_RV=round(100 *(phi./sum(phi)));

% 3.4  Plot the PCA of the Assessors (RV) 
nfig=0;
nfig=nfig+1;
figure(nfig);clf
%subplot(131)
titre=['PCA of R_V Matrix ',...
    '\lambda_1=',num2str(phi(1)),...
    ' \tau_1=',int2str(tau_RV(1))....
    '\lambda_2=',num2str(phi(2)),...
    ' \tau_2=',int2str(tau_RV(2))];
plotxyabs(G,1,2,titre,nom_juge')

boys=find(des_juges(:,1)==1);
mean_boys=mean(G(boys,:));
girls=find(des_juges(:,1)==2);
mean_girls=mean(G(girls,:));
text(mean_boys(1),mean_boys(2),'B');
text(mean_girls(1),mean_girls(2),'G');

if printfigure==1
 print('-depsc',[gen_name,'BetStudies.eps']);
end 

disp('Compute: weights for the compromise');
% 3.3.1 Compute the alpha weights
weights=P(:,1)/sum(abs(P(:,1)))  ;

% 3.4. Commpute the compromise
compromise=zeros(nbeers,nbeers);
for k=1:njuges;
  eval(['compromise=compromise+weights(k)*S',...
      int2str(k),';'])
end

% 3.5 PCA of the compromise
[pc,lc]=eigen(compromise);

disp('Compute: coordinates of the obs on the compromise')
F=pc*diag(sqrt(lc))

disp('eigenvalues of the compromise');
disp('    Eigen    Sum_Eig       %     sum % ')
[ lc cumsum(lc) lc./sum(lc) cumsum( lc./sum(lc)  ) ]

   minF=min(F);
   maxF=max(F);
   
   
% 3.6. Project supplementary variables
% Biplots: Project the Alcohol and Hedonic
%          as supplementary variables
%
% Rsup4I are the correlations (loadings)
% between the supplementary variables for the I set
% and the Factor Scores for the I set
% They could be used for creating 
% the circle of correlation
Rsup4I=corr2mat(des_bieres,F);
% ProjSup4I are scaled/normalized by
% half of the singular value. 
% This mean that a correlation of 1
% would correspond to the maximum
% projection of a variable
% This normalization may need more work
%    (or even some thinking ... Gasp!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ProjSup4I=Rsup4I*diag(lc.^(1/2));
ProjSup4I=Rsup4I*diag(lc.^(1/2)/2);

% 3.6.1. Plot all that
nfig=nfig+1;
figure(nfig);clf
ax1=1;ax2=2;
% ax1=3;ax2=4;
% 
plotxyrange(F,ax1,ax2,...
   'Compromise as barycenter of Studies',...
   char(Nom_de_row),...
   [minF(ax1),maxF(ax1),minF(ax2),maxF(ax2)]);
if printfigure==1
   print('-depsc ',...
       [gen_name,'CompAsBary',int2str(ax1),int2str(ax2),'.eps'])
end
%% Project the sup var as vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nJsup=length(ProjSup4I(:,1));
for is=1:nJsup;
plot([0,ProjSup4I(is,ax1)],[0,ProjSup4I(is,ax2)],'->r')
text(ProjSup4I(is,ax1),ProjSup4I(is,ax2),int2str(is),'Color','r')
end   
   
% 
% 3.7. Now Go for the assessort   
%      Project the individual assessors/solutions onto the 
%      compromise
disp('Compute: Projection Operator ProjOp=P*Lambda^(-1/2)')
 RProj=pc*diag(lc.^(-1/2)) ;
 for k=1:njuges;
  eval(['RF',int2str(k),'=S',int2str(k),'* RProj;']); 
  eval(['minF=min([minF;RF',int2str(k),']);']);
  eval(['maxF=max([maxF;RF',int2str(k),']);']); 
 end

% 3.7.1 Make common plot
code_couleur=['r','m','y','b','g','c','k']; 
nfig=nfig+1;
figure(nfig);clf;
plotxyrange(F,ax1,ax2,...
   'Compromise as barycenter of Studies',char(Nom_de_row),...
   [minF(ax1),maxF(ax1),minF(ax2),maxF(ax2)])
mark=code_couleur';
for t=1:njuges;
   eval(['plotline(F,RF',int2str(t),...
      ' ,ax1,ax2,mark(mod(t-1,7)+1,:))'])
end
if printfigure==1
  print('-depsc ', [gen_name,'CompProjAlgo.eps'])
end
% 3.7.2. Make the convex hulls for the beers 
% for Dimensions 1 and 2
axisbeers=zeros(nbeers,njuges,2);
for i=1:nbeers; 
    for t=1:njuges;
        eval(['axisbeers(i,t,:)=RF',int2str(t),...
            '(i,1:2);']);
    end
end

nfig=nfig+1;
figure(nfig);clf;
plotxyrange(F,ax1,ax2,...
   'Compromise & convex hulls',char(Nom_de_row),...
   [minF(ax1),maxF(ax1),minF(ax2),maxF(ax2)])
mark=code_couleur';

for i=1:nbeers
  coul=mark(mod(i-1,7)+1,:); 
  x=axisbeers(i,:,1);
  y=axisbeers(i,:,2);
  KK=convhull(x,y);
  plot(x(KK),y(KK),[coul,'-'],x,y,[coul,'+'],'LineWidth',2)
  lenom=char(Nom_de_row(i));
  text(F(i,ax1),F(i,ax2),lenom,'Color',coul)
  K2=length(KK);
  for k2=1:K2;
       plot([x(KK(k2)),F(i,ax1)],[y(KK(k2)),F(i,ax2)],...
         [coul,':'],'LineWidth',1.5)
  end
end








