cccccccccc FORTRAN subroutine gSlcMC cccccccccc

c For doing slice sampling for generalised additive
c mixed model structures.

c Last changed: 19 MAY 2018

      subroutine gslcmc(y,Cspl,CsplTy,Ahyp,ibkstt,ibkend,nMCMC,numObs,
     +                  ncCspl,idPres,idnum,numGps,idStt,idEnd,nVec,
     +                  xnu,sigsq,xnucur,lennu,lenbki,lenssq,lnumn1,
     +                  ifam,gparm,xncnoj,buSpl,betauO,idxSpl,uSbj,
     +                  uSbjO,idxSbj,aux,ssqnu,msgCod)

      double precision y(numObs),Cspl(numObs,ncCspl),CsplTy(numObs),  
     +                 xnu(lennu,nMCMC),sigsq(lenssq,nMCMC),xnucrj,
     +                 xlgdcr,xnucur(lennu),xncnoj(lnumn1),
     +                 buSpl(ncCspl),betauO(lennu),uSbj(numGps),
     +                 uSbjO(numGps),sumDts,currPrc,a,w,Asq,ans,
     +                 unfvar,gamvar,gparm(2),xLnuj,xUnuj,curTsh,
     +                 xLnujt,xUnujt,xnuNwj,aux(lenssq,nMCMC),shpcur,
     +                 ratcur,unorm,ssqnu(lennu),Ahyp(lenssq)

      integer ibkstt(lenbki),ibkend(lenbki),numObs,ncCspl,nMCMC,lennu,
     +        lenbki,lenssq,lnumn1,idPres,idnum(numGps),numGps,nVec,
     +        ipcCnt,idxSpl(ncCspl),idxSbj(numGps),ig,j,k,ipos,newptf,
     +        msgCod,incmnt

      sumDts = 0.0
      ipcCnt = 0

c     Make possible adjustments to the message code variable:

      if (nMCMC.lt.100.and.msgCod.eq.1) then
         msgCod = 0
      endif
      if (nMCMC.lt.100.and.msgCod.eq.2) then
         msgCod = 0
      endif
      if (nMCMC.lt.10.and.msgCod.eq.3) then
         msgCod = 0
      endif

c     Perform Markov chain Monte Carlo iterations:

      if (msgCod.gt.0) then
          call intpr("",0,msgCod,0)
      endif

      do 10 ig = 2,nMCMC

         currPrc = 100.0*ig/nMCMC

         if (currPrc.ge.ipcCnt) then
            if (msgCod.gt.0) then
               if (ipcCnt.gt.0) then
                  if (ipcCnt.lt.10) then
                      call intpr("Percentage completed:",21,ipcCnt,1)              
                  endif
                  if (ipcCnt.ge.10.and.ipcCnt.le.99) then
                      call intpr("Percentage completed:",21,ipcCnt,1)              
                  endif
                  if (ipcCnt.eq.100) then
                      call intpr("Percentage completed:",21,ipcCnt,1)
                      call intpr("",0,msgCod,0)  
                  endif
                  if (ipcCnt.eq.10.and.msgCod.eq.1) then
                      call intpr("======================",22,msgCod,0)  
                      call intpr("From now on will only ",22,msgCod,0)
                      call intpr("flag multiples of 10%.",22,msgCod,0)
                      call intpr("======================",22,msgCod,0)  
                  endif
               endif
               if (msgCod.eq.1) then
                  if (ipcCnt.lt.10) then
                     incmnt = 1
                  endif
                  if (ipcCnt.ge.10) then
                     incmnt = 10
                  endif
               endif
               if (msgCod.eq.2) then
                   incmnt = 1
               endif
               if (msgCod.eq.3) then
                   incmnt = 10
               endif
               ipcCnt = ipcCnt + incmnt
            endif
         endif
 
c Update the nu chains:

         do 20 j = 1,lennu
            xnucur(j) = xnu(j,ig-1)
20       continue

          do 30 j = 1,lennu

            xnucrj = xnucur(j)
            ipos = 1
            do 40 k = 1,lennu
               if (k.ne.j) then
                  xncnoj(ipos) = xnucur(k)
                  ipos = ipos + 1 
               endif
40          continue

            ans = 0.0
            call lgunds(j,xnucrj,xncnoj,lennu,lnumn1,Cspl,y,CsplTy,
     +                  ncCspl,numObs,ssqnu,idPres,nVec,idnum,numGps,
     +                  idStt,idEnd,ifam,buSpl,betauO,idxSpl,uSbj,
     +                  uSbjO,idxSbj,ans)  
            xlgdcr = ans
            call urand(1,unfvar)
            a = xlgdcr + log(unfvar)

            if (ig.eq.2) then
               w = 1.0
            endif

            if (ig.gt.2) then
               w = sumDts/(ig-2)
            endif

            call urand(1,unfvar)

            xLnuj = xnucrj - w*unfvar
            xUnuj = xLnuj + w

            ans = 0.0
            call lgunds(j,xLnuj,xncnoj,lennu,lnumn1,Cspl,y,CsplTy,
     +                  ncCspl,numObs,ssqnu,idPres,nVec,idnum,numGps,
     +                  idStt,idEnd,ifam,buSpl,betauO,idxSpl,uSbj,
     +                  uSbjO,idxSbj,ans)  
            curTsh = ans

            if (a.lt.curTsh) then
50             xLnuj = xLnuj - w
               ans = 0.0
               call lgunds(j,xLnuj,xncnoj,lennu,lnumn1,Cspl,y,CsplTy,
     +                     ncCspl,numObs,ssqnu,idPres,nVec,idnum,
     +                     numGps,idStt,idEnd,ifam,buSpl,betauO,idxSpl,
     +                     uSbj,uSbjO,idxSbj,ans)  
               curTsh = ans
               if (a.lt.curTsh) then
                  goto 50
               endif                      
            endif

            ans = 0.0
            call lgunds(j,xUnuj,xncnoj,lennu,lnumn1,Cspl,y,CsplTy,
     +                  ncCspl,numObs,ssqnu,idPres,nVec,idnum,numGps,
     +                  idStt,idEnd,ifam,buSpl,betauO,idxSpl,uSbj,
     +                  uSbjO,idxSbj,ans)

            curTsh = ans
 
            if (a.lt.curTsh) then
60             xUnuj = xUnuj + w
               ans = 0.0
               call lgunds(j,xUnuj,xncnoj,lennu,lnumn1,Cspl,y,CsplTy,
     +                     ncCspl,numObs,ssqnu,idPres,nVec,idnum,
     +                     numGps,idStt,idEnd,ifam,buSpl,betauO,idxSpl,
     +                     uSbj,uSbjO,idxSbj,ans)  
               curTsh = ans
               if (a.lt.curTsh) then
                  goto 60
               endif                      
            endif

            xLnujt = xLnuj
            xUnujt = xUnuj

            newptf = 0

70             call urand(1,unfvar)

               xnuNwj = xLnujt + unfvar*(xUnujt - xLnujt)

               ans = 0.0
               call lgunds(j,xnuNwj,xncnoj,lennu,lnumn1,Cspl,y,CsplTy,
     +                     ncCspl,numObs,ssqnu,idPres,nVec,idnum,
     +                     numGps,idStt,idEnd,ifam,buSpl,betauO,idxSpl,
     +                     uSbj,uSbjO,idxSbj,ans)  
               curTsh = ans
               if (a.lt.curTsh) then
                  xnu(j,ig) = xnuNwj
                  xnucur(j) = xnuNwj
                  newptf = 1
               endif

               if (xnuNwj.lt.xnucrj) then
                  xLnujt = xnuNwj
               endif

               if (xnuNwj.ge.xnucrj) then
                  xUnujt = xnuNwj
               endif

               if (newptf.eq.0) then
                  goto 70
               endif

            sumDts = sumDts + abs(xnuNwj - xnucrj)
         
30       continue

c Update the sigsq chains:

         do 80 j = 1,lenssq

            Asq = Ahyp(j)*Ahyp(j)
            
            gparm(1) = 1.0
            gparm(2) = 1.0
            call gammarand(1,gparm,gamvar)

            aux(j,ig) = ((1/sigsq(j,ig-1)) + (1/Asq))/gamvar
            
            unorm = 0.0
            do 90 k = ibkstt(j+1),ibkend(j+1)
                unorm = unorm + xnu(k,ig)*xnu(k,ig)
90          continue

            shpcur =  0.5*(ibkend(j+1) - ibkstt(j+1)) + 1.0
            ratcur = 0.5*unorm + (1/aux(j,ig))

            gparm(1) = shpcur
            gparm(2) = 1.0
            call gammarand(1,gparm,gamvar)

            sigsq(j,ig) = ((1/aux(j,ig)) + 0.5*unorm)/gamvar

            do 100 k = ibkstt(j+1),ibkend(j+1)
                ssqnu(k) = sigsq(j,ig)
100         continue

80       continue

10    continue

      return
      end

cccccccccc End of gSlcMC cccccccccc
