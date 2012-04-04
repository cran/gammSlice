ccccccccc Fortran subroutine: logfbi cccccccccccccccccc

c For computing the logarithm of
c the unnormalized density function
c for slice sampling for GLMM

c Last changed: 24 MAY 2010

      subroutine logfbi(nu,j,CTy,Cmat,Cnurest,sigsqnu,ncC,n,logf)

      double precision nu,Cnurest(n),Cmat(n,ncC),sigsqnu(ncC),
     +                 logf,CTy(ncC)
      integer i,j,n,ncC

      logf = CTy(j)*nu - nu*nu/(2.0*sigsqnu(j))
  
      do 10 i = 1,n
         logf = logf - log(1+exp(Cmat(i,j)*nu + Cnurest(i)))     
10    continue  

      return
      end

cccccccccc End of subroutine logfbi cccccccccc
