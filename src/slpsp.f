ccccccccc Fortran subroutine: slpsp.f cccccccccccccccccc

c For implementing a Gibbs sampling loop for 
c slice sampling for general X and Z matrices.
c Note that Cmat=[X Z] is inputted.


c Last changed: 21 JAN 2012


      subroutine slpsp(y,Cmat,CTy,nobs,ncC,ncX,ngps,nreb,rbei,ngb,
     +     ubkid,lbkid,lennu,lenssu,sumdist,sigsqnu,nucurr,Cnurest,
     +     Su,nu,sigsqu)

      integer nobs,ncC,ncX,ngps,nreb,rbei(nreb),ngb,
     +     ubkid(ngps),lbkid(ngps),lennu,lenssu,
     +     ig,i,ind,irb,klow,kupp,k,idingp,glow,gupp,kgap,gpid

      double precision y(nobs),Cmat(nobs,ncC),CTy(ncC),sumdist(ncC),
     +     sigsqnu(ncC),nucurr(ncC),Cnurest(nobs),Su(nreb),
     +     nu(lennu),sigsqu(lenssu),gammaparm(2),
     +     nucurrj,nujnew,logfcr,uvar,a,w,Lnuj,Unuj,logfval,
     +     Lnujt,Unujt,unorm,nujold,gvar,sum,scalgvar

      do 10 ig = 2,ngb
         do 20 j = 1,ncC
            ind = (ig-2)*ncC + j
            nucurr(j) = nu(ind)
 20      continue
        
c     Set vector denoted by "Cnurest"
        
         do 30 idingp = 1,nobs
            Cnurest(idingp) = 0.0
            do 40 j = 2,ncC
               ind = (ig-2)*ncC + j
               Cnurest(idingp) = Cnurest(idingp) + 
     +              Cmat(idingp,j)*nu(ind)
 40         continue
 30      continue
            
c     Update the Cnurest
         
         do 50 j = 1,ncC
            nucurrj = nucurr(j)
            
            if(j.gt.1)then
               do 60 idingp = 1,nobs
                  Cnurest(idingp) = Cnurest(idingp) +
     +                 Cmat(idingp,j-1)*nujnew - Cmat(idingp,j)*nucurrj
 60            continue
            endif
            
            call logfps(nucurrj,j,CTy,Cmat,Cnurest,sigsqnu,ncC,
     +           nobs,logfcr)
            call urand(1,uvar)
            a =  logfcr + log(uvar)
            
            if(ig.eq.2)then
               w = 1.0
            else
               w = sumdist(j)/(ig-2)
            endif
            
            call urand(1,uvar)
            
            Lnuj = nucurrj - w*uvar
            Unuj = Lnuj + w
            
            
 70         call logfps(Lnuj,j,CTy,Cmat,Cnurest,sigsqnu,ncC,
     +           nobs,logfval)
            if(a.lt.logfval)then
               Lnuj = Lnuj - w
               goto 70
            endif
 80        call logfps(Unuj,j,CTy,Cmat,Cnurest,sigsqnu,ncC,
     +           nobs,logfval)
            if(a.lt.logfval)then
               Unuj = Unuj + w
               goto 80
            endif
            

            Unujt = Unuj
            Lnujt = Lnuj
            
 90         call urand(1,uvar)
            nujnew = Lnujt + uvar*(Unujt-Lnujt)
            
            call logfps(nujnew,j,CTy,Cmat,Cnurest,sigsqnu,ncC,
     +           nobs,logfval)
            if (a.lt.logfval)then
               ind = (ig-1)*ncC + j
               nu(ind) = nujnew
               nucurr(j) = nujnew
               goto 100
            endif
            
            if(nujnew.lt.nucurrj)then
               Lnujt = nujnew
            else
               Unujt = nujnew
            endif
            goto 90

 100        ind = (ig-2)*ncC + j
            nujold = nu(ind)
            sumdist(j) = sumdist(j) + abs(nujnew-nujold)
        
 50      continue
        
c     Update variance components
         if (nreb.gt.0) then
            do 110 irb = 1, nreb
               unorm = 0.0
               
               if (irb.eq.1) then
                  klow = ncX + 1
               else
                  klow = ncX + rbei(irb-1) + 1
               endif
               kupp = ncX + rbei(irb)
               
               do 120 k = klow,kupp
                  ind = (ig-1)*ncC + k
                  unorm = unorm + nu(ind)*nu(ind)
 120           continue
               
               ind = (ig-2)*nreb + irb
               scalgvar = 1/sigsqu(ind) + 1/(Su(irb)**2)

               gammaparm(1) = 1
               gammaparm(2) = 1
               call gammarand(1,gammaparm,gvar)
               gvar = gvar/scalgvar
 
               kgap = kupp - klow
               scalgvar = gvar + 0.5 * unorm

               gammaparm(1) = (1 + kgap)/2
               gammaparm(2) = 1
               call gammarand(1,gammaparm,gvar)
               
               ind = (ig-1)*nreb + irb
               sigsqu(ind) = scalgvar/gvar
               
               do 130 k = klow,kupp
                  sigsqnu(k) = sigsqu(ind)
 130           continue
            
 110        continue
         endif
         
 10   continue

      return
      end

ccccccccccccccccccccccc-slpsp.f-ccccccccccccccccccc
