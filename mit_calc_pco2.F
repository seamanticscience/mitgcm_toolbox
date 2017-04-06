#include "fintrf.h"
C======================================================================
C     CADVECT calls advection schemes to calculate scalar transports
C     [pco2oc,co2aq,hco3,co3,pH,iter]=mit_calc_pco2(theta,salt,pressc,mask,
C								             dic,alk,po4,sit);
C
C======================================================================
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer*4 nlhs, nrhs

C     Function declarations:
      mwPointer mxCreateNumericArray
      integer*4 mxClassIDFromClassName
      integer mxIsDouble, mxIsNumeric
      mwPointer mxGetM, mxGetN
      mwSize mxGetNumberOfDimensions
      mwPointer mxGetDimensions
      mwPointer mxGetPr

C     Arguments for mxCreateNumericArray
      integer*4 classid, complexflag
      mwSize ndim, dims(4)
    
C     Pointers to input/output mxArrays:
      mwPointer pco2oc, co2aq, hco3, co3, ph, iter
      mwPointer theta, salt, pressc, mask
      mwPointer dic, alk, po4, sit

C     Array information:
      mwPointer mrows, ncols
      mwSize size

CC    For writing error to the screen
C     mexErrMsgIdAndTxt causes crash
      character*250 line
      integer*4 k
      integer*4 mexPrintf
      integer*4 errCode

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 8) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nInput     ',
C     &                           'Eight inputs required. ')
         write(line,*)'ERROR: Fourteen inputs required.'
         k=mexPrintf(line//achar(13)) 
         return
      elseif(nlhs .gt. 6) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nOutput    ',
C     &                           'Too many output arguments.')
         write(line,*)'ERROR: Too many output arguments.'
         k=mexPrintf(line//achar(13)) 
         return
      endif

C     Validate inputs
C     Check that the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:NonNumeric',
C     &                           'Input must be a number.  ')
         write(line,*)'ERROR: Input must be numeric.'
         k=mexPrintf(line//achar(13)) 
         return     
      endif
      
C     Check data type of input argument.
      if (mxIsDouble(prhs(1)) .eq. 0) then
C         call mexErrMsgTxt('Input argument must be of type double.')
         write(line,*)'ERROR: Input must be double precision.'
         k=mexPrintf(line//achar(13)) 
         return     
      endif
      
C     Get the size of the input arrays using temperature as template.
C     Use mxGetDimensions to get dimensions into dims
      ndim = mxGetNumberOfDimensions(prhs(1))
      call mxCopyPtrToInteger4(mxGetDimensions(prhs(1)),  
     &           dims, mxGetNumberOfDimensions(prhs(1)))

C      write(line,*)'ndims: ',ndim
C      k=mexPrintf(line//achar(13)) 
C      write(line,*)'dimensions: ',dims(1),'x',dims(2),
C     &                        'x',dims(3),'x',dims(4)
C      k=mexPrintf(line//achar(13)) 

C     Create Fortran arrays from the input arguments.
      theta = mxGetPr(prhs(1))
      salt  = mxGetPr(prhs(2))
      pressc= mxGetPr(prhs(3))
      mask  = mxGetPr(prhs(4))
      dic   = mxGetPr(prhs(5))
      alk   = mxGetPr(prhs(6))
      po4   = mxGetPr(prhs(7))
      sit   = mxGetPr(prhs(8))
      
C     Create matrix for the return argument.
      classid = mxClassIDFromClassName('double')
      complexflag = 0

      plhs(1) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)
     
      plhs(2) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)
     
      plhs(3) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)
     
      plhs(4) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)
     
      plhs(5) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)
     
      plhs(6) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)

      pco2oc= mxGetPr(plhs(1))
      co2aq=mxGetPr(plhs(2))
      hco3= mxGetPr(plhs(3))
      co3 = mxGetPr(plhs(4))
      pH  = mxGetPr(plhs(5))
      iter= mxGetPr(plhs(6))

CC     Call the computational subroutines.
C calc_pco2 calculates the pco2 and ph, as well as carbon species
         call calc_pco2(dims(1),dims(2),dims(3),
     I           %val(mask),%val(pressc),
     I           %val(theta),%val(salt),%val(dic),
     I           %val(po4),%val(sit),%val(alk),
     O           %val(ph),%val(pco2oc),%val(co2aq),
     O           %val(hco3),%val(co3),%val(iter))          
      return
      end