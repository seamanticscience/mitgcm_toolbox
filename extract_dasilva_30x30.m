function [T,x,y]=extract_dasilva_30x30(suff,months,lon_lo,lon_hi,lat_lo,lat_hi)
%
% Extract daSilva data from archive file '/d01/adcroft/dasilva/30x30/raw/*.r4'
%   1 degree by 1 degree dataset
%
% e.g.
%   [Tx,x,y]=extract_dasilva_30x30('Taux',1,-180,180,-90,90);
%   [Qnet,x,y]=extract_dasilva_30x30('Qnet',1:12,0,360,-90,90);

archivedirectory

fid=fopen([ARCHIVEDIR '/dasilva/30x30/raw/' suff '.r4'],'r','b');
if fid <0 
 error('I failed to open the file')
end

% Loop over all the month arguments
jjrec=0;
for month=months,
jjrec=jjrec+1;

status=fseek(fid,(month-1)*(360*2*180*2*4),'bof');
if status ~= 0
 sprintf(' fseek status=%i',status);
 error(ferror(status))
end

i0=round( lon_lo*2 )+1;
i1=round( lon_hi*2 );
j0=round( (lat_lo+90)*2 )+1;
j1=round( (lat_hi+90)*2 );

if i1>i0
 ii=i0:i1;
else
 ii=[i0:360*2 1:i1];
end
ii=mod(ii+360*2-1,360*2)+1;
jj=j0:j1;

nx=prod(size(ii));
ny=prod(size(jj));

LD=fread(fid,360*2*180*2,'real*4');
nodata=LD(1);
LD( find(LD==nodata) )=NaN;
LD=reshape( LD, [360*2 180*2]);

T(:,:,jjrec)=LD(ii,jj,:);

end
fclose(fid);

xl=0.25:0.5:359.75;
if lon_lo<0 | lon_hi<0
 xl=mod(xl+180,360)-180;
end
yl=-89.75:0.5:89.75;

x=xl(ii);
y=yl(jj);
