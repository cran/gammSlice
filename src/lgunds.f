cccccccccc FORTRAN subroutine lgunds cccccccccc

c For evaluating the logarithm of the unnormalised density
c function of the `nu' vector that arises in slice sampling
c for generalised additive mixed models.

c Last changed: 27 APR 2018

      subroutine lgunds(j,xnucrj,xncnoj,lennu,lnumn1,Cspl,y,CsplTy,
     +                  ncCspl,numObs,ssqnu,idPres,nVec,idnum,numGps,
     +                  idStt,idEnd,ifam,buSpl,betauO,idxSpl,uSbj,uSbjO,
     +                  idxSbj,ans)

      double precision xnucrj,xncnoj(lnumn1),y(numObs),CsplTy(numObs),
     +                 Cspl(numObs,ncCspl),ssqnu(lennu),buSpl(ncCspl),
     +                 betauO(lennu),uSbj(numGps),uSbjO(numGps),
     +                 suminn,term,ans,fstTrm,scdTrm,thdTrm

      integer i,j,jd,k,lennu,lnumn1,numObs,ncCspl,numGps,ifam,
     +        idStt(numGps),idEnd(numGps),idPres,idxSpl(ncCspl),
     +        idxSbj(numGps),nVec(numGps),idnum(numObs)

c     Obtain the first term:

      if (j.le.ncCspl) then
         ans = CsplTy(j)*xnucrj
      endif
      if (idPres.eq.1) then
         if (j.gt.ncCspl) then
            jd = j - ncCspl
            term = 0.0
            do 10 k = idStt(jd),idEnd(jd)
               term = term + y(k)   
10          continue
            ans = term*xnucrj
         endif
      endif

      fstTrm = ans

c     Obtain the second term, and update the answer starting with determination 
c     of the current `buSpl' vector and, if required the current 'uSbj' vector:

      if (j.le.ncCspl) then

         betauO(1) = xnucrj
         do 20 k = 1,lnumn1
            betauO(k+1) = xncnoj(k)
20       continue
  
         if (j.eq.1) then
            do 30 k = 1,ncCspl
               idxSpl(k) = k 
30          continue
         endif
         if (j.gt.1.and.j.lt.ncCspl) then
            do 40 k = 1,(j-1)
               idxSpl(k) = k + 1
40          continue
            idxSpl(j) = 1
            do 50 k = (j+1),ncCspl
               idxSpl(k) = k 
50          continue
         endif
         if (j.eq.ncCspl) then
            do 60 k = 1,(ncCspl-1)
               idxSpl(k) = k + 1 
60          continue
            idxSpl(ncCspl) = 1
         endif

         do 70 k = 1,ncCspl
            buSpl(k) = betauO(idxSpl(k))
70       continue  

         if (idPres.eq.1) then
             do 80 k = 1,numGps
                uSbj(k) = betauO(ncCspl+k)
80           continue
         endif
      endif

      if (idPres.eq.1) then
         if (j.gt.ncCspl) then
            do 90 k = 1,ncCspl
               buSpl(k) = xncnoj(k)
90          continue

            uSbjO(1) = xnucrj
            do 100 k = 2,numGps         
               uSbjO(k) = xncnoj(ncCspl+k-1)
100         continue

            jd = j - ncCspl 
            if (jd.eq.1) then
               do 110 k = 1,numGps
                  idxSbj(k) = k               
110            continue
            endif
            if (jd.gt.1.and.jd.lt.numGps) then
               do 120 k = 1,(jd-1)
                  idxSbj(k) = k + 1
120            continue
               idxSbj(jd) = 1
               do 130 k = (jd+1),numGps
                  idxSbj(k) = k 
130            continue
            endif 
            if (jd.eq.numGps) then
               do 140 k = 1,(numGps-1)
                  idxSbj(k) = k + 1 
140            continue
               idxSbj(numGps) = 1
            endif

            do 150 k = 1,numGps
               uSbj(k) = uSbjO(idxSbj(k))
150         continue  
         endif
      endif

      scdTrm = 0.0

      if (ifam.eq.1) then
         do 160 i = 1,numObs
            term = 0.0
            do 170 k = 1,ncCspl
               term = term + Cspl(i,k)*buSpl(k)
170         continue
          if (idPres.eq.1) then
             term = term + uSbj(idnum(i))   
          endif
          if (term.lt.100) then
             term = log(1.0 + exp(term))
          endif
          scdTrm = scdTrm - term
160       continue
      endif

      if (ifam.eq.2) then
         do 180 i = 1,numObs
            term = 0.0
            do 190 k = 1,ncCspl
               term = term + Cspl(i,k)*buSpl(k)
190         continue
            if (idPres.eq.1) then
               term = term + uSbj(idnum(i))
            endif
         scdTrm = scdTrm - exp(term)
180      continue
      endif

c     Obtain the third term:

      thdTrm = -0.5*xnucrj*xnucrj/ssqnu(j)

c     Combine all three terms for the answer:   

      ans = fstTrm + scdTrm + thdTrm
 
      return
      end

cccccccccc End of lgunds cccccccccc
