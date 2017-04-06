% WRDA( file, arr, irec [options] )
%
% This routine writes the array 'arr' as a record in a direct-access
% binary file (float or double precision) with filename 'file'.
% The data will be written as the irec'th record. The file is created
% if it does not exists.
%
% Required arguments:
%
%   file  - string  - name of file to create/append
%   arr   - array   - matlab matrix or multi-dimensional array
%   irec  - integer - record number in file in which to write data
%
% Optional arguments (must appear after the required arguments):
%   prec  - string  - precision of storage in file. Default = 'real*8'.
%   ieee  - string  - IEEE bit-wise representation in file. Default = 'b'.
%
% 'prec' may take the values:
%       'real*4' - floating point, 32 bits.
%       'real*8' - floating point, 64 bits - the default.
%
% 'ieee' may take values:
%    'ieee-be'     or 'b' - IEEE floating point with big-endian
%                           byte ordering - the default
%    'ieee-le'     or 'l' - IEEE floating point with little-endian
%                           byte ordering
%    'native'      or 'n' - local machine format
%    'cray'        or 'c' - Cray floating point with big-endian
%                           byte ordering
%    'ieee-le.l64' or 'a' - IEEE floating point with little-endian
%                           byte ordering and 64 bit long data type
%    'ieee-be.l64' or 's' - IEEE floating point with big-endian byte
%                           ordering and 64 bit long data type.
%
% eg.   wrda('t.data',T,4);
%       wrda('t.data',T,4,'real*4');
%       wrda('t.data',T,4,'real*4','b');
function [] = wrda(file,arr,k,varargin)

% Defaults
WORDLENGTH=8;
rtype='real*8';
ieee='b';

% Check optional arguments
args=char(varargin);
while (size(args,1) > 0)
 if deblank(args(1,:)) == 'real*4'
  WORDLENGTH=4;
  rtype='real*4';
 elseif deblank(args(1,:)) == 'real*8'
  WORDLENGTH=8;
  rtype='real*8';
 elseif deblank(args(1,:)) == 'n' | deblank(args(1,:)) == 'native'
  ieee='n';
 elseif deblank(args(1,:)) == 'l' | deblank(args(1,:)) == 'ieee-le'
  ieee='l';
 elseif deblank(args(1,:)) == 'b' | deblank(args(1,:)) == 'ieee-be'
  ieee='b';
 elseif deblank(args(1,:)) == 'c' | deblank(args(1,:)) == 'cray'
  ieee='c';
 elseif deblank(args(1,:)) == 'a' | deblank(args(1,:)) == 'ieee-le.l64'
  ieee='a';
 elseif deblank(args(1,:)) == 's' | deblank(args(1,:)) == 'ieee-be.l64'
  ieee='s';
 else
  sprintf(['Optional argument ' args(1,:) ' is unknown'])
  return
 end
 args=args(2:end,:);
end

N=size(arr);
nnn=prod(N);

[fid mess]=fopen(file,'r+',ieee);
if fid == -1
 [fid mess]=fopen(file,'w',ieee);
 if fid == -1
  sprintf('Error while opening file:\n%s',mess)
  return
 end
end
st=fseek(fid,nnn*(k-1)*WORDLENGTH,'bof');
if st ~= 0
 mess=ferror(fid);
 sprintf('There was an error while positioning the file pointer:\n%s',mess)
 return
end
count=fwrite(fid,arr,rtype);
if count ~= nnn
 sprintf('Not enough data was available to be read: off EOF?')
 return
end
st=fclose(fid);
