ccccccccc Fortran subroutine: slbsp.f cccccccccccccccccc

c For implementing a Gibbs sampling loop for 
c slice sampling for general X and Z matrices.
c Note that Cmat=[X Z] is inputted.


c Last changed: 21 JAN 2012


      subroutine slbsp(y,Cmat,CTy,nobs,ncC,ncX,ngps,nreb,rbei,ngb,
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
       
         do 30 gpid = 1,ngps
            glow = lbkid(gpid)
            gupp = ubkid(gpid)
           
            do 40 idingp = glow,gupp
               Cnurest(idingp) = 0.0
               do 50 j = 2,ncC
                  ind = (ig-2)*ncC + j
                  Cnurest(idingp) = Cnurest(idingp) + 
     +                 Cmat(idingp,j)*nu(ind)
 50            continue
 40         continue
 
 30      continue   
c     Update the 
     
         do 70 j = 1,ncC
            nucurrj = nucurr(j)
            
            if(j.gt.1)then
               do 80 idingp = 1,nobs
                  Cnurest(idingp) = Cnurest(idingp) +
     +                 Cmat(idingp,j-1)*nujnew - Cmat(idingp,j)*nucurrj
 80            continue
            endif
            
            call logfbi(nucurrj,j,CTy,Cmat,Cnurest,sigsqnu,ncC,
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
            
            
 90         call logfbi(Lnuj,j,CTy,Cmat,Cnurest,sigsqnu,ncC,
     +           nobs,logfval)
            if(a.lt.logfval)then
               Lnuj = Lnuj - w
               goto 90
            endif
 100        call logfbi(Unuj,j,CTy,Cmat,Cnurest,sigsqnu,ncC,
     +           nobs,logfval)
            if(a.lt.logfval)then
               Unuj = Unuj + w
               goto 100
            endif
            

            Unujt = Unuj
            Lnujt = Lnuj
            
 110        call urand(1,uvar)
            nujnew = Lnujt + uvar*(Unujt-Lnujt)
            
            call logfbi(nujnew,j,CTy,Cmat,Cnurest,sigsqnu,ncC,
     +           nobs,logfval)
            if (a.lt.logfval)then
               ind = (ig-1)*ncC + j
               nu(ind) = nujnew
               nucurr(j) = nujnew
               goto 120
            endif
            
            if(nujnew.lt.nucurrj)then
               Lnujt = nujnew
            else
               Unujt = nujnew
            endif
            goto 110

 120        ind = (ig-2)*ncC + j
            nujold = nu(ind)
            sumdist(j) = sumdist(j) + abs(nujnew-nujold)
          
 70      continue
        
c     Update variance components
         do 130 irb = 1, nreb
            unorm = 0.0
            
            if (irb.eq.1) then
               klow = ncX + 1
            else
               klow = ncX + rbei(irb-1) + 1
            endif
            kupp = ncX + rbei(irb)
            
            do 140 k = klow,kupp
               ind = (ig-1)*ncC + k
               unorm = unorm + nu(ind)*nu(ind)
 140        continue

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

            do 150 k = klow,kupp
               sigsqnu(k) = sigsqu(ind)
 150        continue
            
 130     continue

 10      continue

      return
      end

ccccccccccccccccccccccc-slbsp.f-ccccccccccccccccccc
