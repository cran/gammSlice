ccccccccc Fortran subroutine: logfpsrd cccccccccccccccccc

c For computing the logarithm of
c the unnormalized density function
c for slice sampling for GLMM

c Last changed: 09 SEP 2010

      subroutine logfpsrd(nu,j,CTy,Cmat,Cnurest,sigsqnu,ncC,n,logf,
     +                      ubkid,lbkid,ngps,ncX)

      double precision nu,Cnurest(n),Cmat(n,ncC),sigsqnu(ncC),
     +     logf, CTy(ncC)
      integer i,j,n,ncC,ncX,ngps,idingp,ubkid(ngps),lbkid(ngps)
  
      idingp = j - ncX
      
      if (sigsqnu(j).ne.0) then
	     logf = CTy(j)*nu - nu*nu/(2.0*sigsqnu(j))
      endif

      if (idingp.lt.1.or.idingp.gt.ngps) then 
         do 10 i = 1,n 
            logf = logf - exp(Cmat(i,j)*nu + Cnurest(i))     
 10      continue  
     
      else if (idingp.eq.1) then
         do 20 i = 1,ubkid(1)
            logf = logf - exp(nu + Cnurest(i))
 20      continue
         do 30 i = lbkid(2),n
            logf = logf - exp(Cnurest(i))
 30      continue
         
      else if (idingp.gt.1.and.idingp.lt.ngps) then
         do 40 i = 1,ubkid(idingp-1)
            logf = logf - exp(Cnurest(i))            
 40      continue
         
         do 50 i = lbkid(idingp),ubkid(idingp)
            logf = logf - exp(nu+Cnurest(i))
 50      continue
         
         do 60 i = lbkid(idingp+1),n
            logf = logf - exp(Cnurest(i))
 60      continue
    
      else 
         do 70 i = 1,ubkid(idingp-1)
            logf = logf - exp(Cnurest(i))
 70      continue
         
         do 80 i = lbkid(ngps),n
            logf = logf - exp(nu+Cnurest(i))
 80      continue
      endif
      
      return
      end

cccccccccc End of subroutine logfpsrd cccccccccc
