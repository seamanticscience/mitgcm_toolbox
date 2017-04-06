clear

setgrid

%h=f77read('bathymetry.r4',[nxc nyc],'real*4','b');
h=rdda('bathymetry.bin',[144 100],1);

%tnx=2*3*3*3
%tny=2*5*5
tnx=36
tny=20

ntx=nxc/tnx;
nty=nyc/tny;

tiles=zeros(ntx,nty);
ntiles=0;
for bj=1:nty,
 jj=(1:tny)+(bj-1)*tny;
 for bi=1:ntx,
  ii=(1:tnx)+(bi-1)*tnx;
  if prod(size( find(h(ii,jj)) ))==0
   h(ii,jj)=NaN;
  else
   ntiles=ntiles+1;
   tiles(bi,bj)=1;
   ht(:,:,ntiles)=h(ii,jj)';
  end
 end
end

acttiles=sum(sum(tiles))

colormap('default')
colormap([[1 1 1]' [0 0 0]' colormap']')

pcol(xc,yc,sq(h,0,-7800)');colorbar

xf(end+1)=xf(end)+dlon(end);
yf(end+1)=yf(end)+dlat(end);
hold on
for bj=1:nty,
 plot([xf(1) xf(end)],[1 1]*yf(1+(bj-1)*tny),'k--')
end
for bi=1:ntx,
 plot([1 1]*xf(1+(bi-1)*tnx),[yf(1) yf(end)],'k--')
end
hold off

xlabel('Longitude')
ylabel('Latitude')
title(sprintf(['%i x %i tiles of size %i x %i'...
      '   Active = %i (%i%%)  Inactive = %i (%i%%)'],...
  ntx,nty,tnx,tny,acttiles,round(100*acttiles/ntx/nty),...
   ntx*nty-acttiles, round(100*(ntx*nty-acttiles)/ntx/nty) ))
