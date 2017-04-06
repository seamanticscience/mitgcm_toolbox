function [bolusU,bolusV,bolusW]=mit_calc_bolus(GM_PsiX,GM_PsiY,grd);
%object:    compute bolus velocty field (bolusU,bolusV,bolusW)
%           from gm streamfunction (GM_PsiX,GM_PsiY).
%input:     GM_PsiX,GM_PsiY
%output:    bolusU,bolusV,bolusW

%gcmfaces_global;
if strcmpi(mit_getparm('data.gmredi','GM_AdvForm'),'false')
   GM_PsiX=GM_PsiX./2;
   GM_PsiY=GM_PsiY./2;
end

nr=length(grd.zc);

%replace NaNs with 0s:
GM_PsiX(isnan(GM_PsiX))=0;
GM_PsiY(isnan(GM_PsiY))=0;

%compute bolus velocity:
GM_PsiX(:,:,nr+1)=0*GM_PsiX(:,:,end);
GM_PsiY(:,:,nr+1)=0*GM_PsiY(:,:,end);
bolusU=0*grd.hfacw;
bolusV=0*grd.hfacs;
for k=1:nr;
    bolusU(:,:,k)=(GM_PsiX(:,:,k+1)-GM_PsiX(:,:,k))/grd.dz(k);
    bolusV(:,:,k)=(GM_PsiY(:,:,k+1)-GM_PsiY(:,:,k))/grd.dz(k);
end;
bolusU=bolusU.*(~isnan(grd.hfacw));
bolusV=bolusV.*(~isnan(grd.hfacs));

%and its vertical part
%   (seems correct, leading to 0 divergence)
tmp_x=GM_PsiX(:,:,1:nr).*repmat(grd.dyg,[1 1 nr]);
tmp_y=GM_PsiY(:,:,1:nr).*repmat(grd.dxg,[1 1 nr]);
%[tmp_x,tmp_y]=exch_UV(tmp_x,tmp_y,grd);
tmp_a=repmat(grd.rac,[1 1 nr]);

tmp_wx=0*tmp_x;
    tmp_wx(1:end-1,:,:)=tmp_x(2:end,:,:)-tmp_x(1:end-1,:,:);
    tmp_wx(end    ,:,:)=tmp_x(1,:,:)-tmp_x(end,:,:);

tmp_wy=0*tmp_y;
    tmp_wy(:,1:end-1,:)=tmp_y(:,2:end,:)-tmp_y(:,1:end-1,:);
    tmp_wy(:,end    ,:)=tmp_y(:,end,:).*0;

tmp_w=(tmp_wx+tmp_wy)./tmp_a;
% tmp_w=( tmp_x(2:end,:,:)-tmp_x(1:end-1,:,:) )+...
%         ( tmp_y(:,2:end,:)-tmp_y(:,1:end-1,:) );
%     tmp_w=tmp_w./tmp_a;
 bolusW=tmp_w.*(~isnan(grd.hfacc));
end

function [FLDU,FLDV]=exch_UV(fldU,fldV,grd);

FLDUtmp=exch_T_N(fldU,grd);
FLDVtmp=exch_T_N(fldV,grd);

FLDU=FLDUtmp;
FLDV=FLDVtmp;

FLDU=FLDUtmp(2:end,2:end-1,:);  
FLDV=FLDVtmp(2:end-1,2:end,:);

end

function [FLD]=exch_T_N(fld,grd,varargin);

%if nargin==2; N=varargin{1}; else; N=1; end;
N=1;
%if ~isfield(grd,'domainPeriodicity'); 
%fprintf('\nexch_T_N_ll.m init: different 1 face configurations may \n');
%fprintf('  differ with respect to domain periodicity. By default gcmfaces \n');  
%fprintf('  assumes 0 periodicity, except that if the first dimension has n*90 \n');
%fprintf('  points it is assumed to be periodic. If this is inadequate, \n');
%fprintf('  you can set domainPeriodicity yourself as explained below.\n\n');
%%domainPeriodicity=[0 0];%no periodidicity
%%domainPeriodicity=[1 0];%1st dimension only is periodic
%%domainPeriodicity=[0 1];%2nd dimension only is periodic
%%domainPeriodicity=[1 0];%both dimensions are periodic
%if mod(size(fld{1},1),90)==0; 
  domainPeriodicity=[1 0];
%else; 
%  grd.domainPeriodicity=[0 0];
%end;
%end;

FLD=fld;
s=size(FLD); s(1:2)=s(1:2)+2*N; FLD=NaN*zeros(s);

n3=max(size(fld,3),1); n4=max(size(fld,4),1);
for k3=1:n3; % depth loop
    for k4=1:n4; % time lopp
        
        f1=fld(:,:,k3,k4);
        nan1=NaN*ones(N,size(fld,2));
        nan2=NaN*ones(size(FLD,1),N);
        %face 1:
        if domainPeriodicity(1); F1=[f1(end-N+1:end,:);f1;f1(1:N,:)];
        else; F1=[nan1;f1;nan1]; end;
        if domainPeriodicity(2); F1=[F1(:,end-N+1:end) F1 F1(:,1:N)];
        else; F1=[nan2 F1 nan2]; end;
        
        %store:
        FLD(:,:,k3,k4)=F1;
        
    end;
end;

end