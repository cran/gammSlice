ccccccccc Fortran subroutine: gammslice.f cccccccccccccccccc

c Routine to call other main subroutines 
c It is named gammSlice
c All paramethers are the same as parameters of 
c other main subroutines.
c The routindex indicates the subroutine needs
c to be called.


c Last changed:  02 MAY 2012


      subroutine gammSlice(rind,y,Cmat,CTy,nobs,ncC,ncX,ngps,nreb,
     +     rbei,ngb,ubkid,lbkid,lennu,lenssu,sumdist,sigsqnu,nucurr,
     +     Cnurest,Su,nu,sigsqu)

      integer rind,nobs,ncC,ncX,ngps,nreb,rbei(nreb),ngb,
     +     ubkid(ngps),lbkids(ngps),lennu,lenssu

      double precision y(nobs),Cmat(nobs,ncC),CTy(ncC),sumdist(ncC),
     +     sigsqnu(ncC),nucurr(ncC),Cnurest(nobs),Su(nreb),
     +     nu(lennu),sigsqu(lenssu)


      if(rind.eq.0) then 
        call slbsp(y,Cmat,CTy,nobs,ncC,ncX,ngps,nreb,rbei,ngb,
     +           ubkid,lbkid,lennu,lenssu,sumdist,sigsqnu,
     +           nucurr,Cnurest,Su,nu,sigsqu)
      endif

      if(rind.eq.1) then 
        call slbsprd(y,Cmat,CTy,nobs,ncC,ncX,ngps,nreb,rbei,ngb,
     +           ubkid,lbkid,lennu,lenssu,sumdist,sigsqnu,
     +           nucurr,Cnurest,Su,nu,sigsqu)
      endif


      if(rind.eq.2) then 
        call slpsp(y,Cmat,CTy,nobs,ncC,ncX,ngps,nreb,rbei,ngb,
     +           ubkid,lbkid,lennu,lenssu,sumdist,sigsqnu,
     +           nucurr,Cnurest,Su,nu,sigsqu)
      endif

      if(rind.eq.3) then 
        call slpsprd(y,Cmat,CTy,nobs,ncC,ncX,ngps,nreb,rbei,ngb,
     +           ubkid,lbkid,lennu,lenssu,sumdist,sigsqnu,
     +           nucurr,Cnurest,Su,nu,sigsqu)
      endif

      return
      end

ccccccccccccccccccccccc-gammslice.f-ccccccccccccccccccc
