cccccccccc Fortran subroutine: frgamma cccccccccc

c Emulates R function rgamma() for gamma random
c number generation inside Fortran 77. Current
c algorithm has problems for the shape parameter
c less than about 0.01.

c Last changed: 19 DEC 2003

      subroutine frgamma(n,alpha,svec)  
      integer n,acc
      double precision alpha,a,b,c,d,p,r,u,v,w,x,y,uf,vr,
     +                 vsmall,svec(n)

c     Set value for smallest possible alpha

      vsmall = 0.01

c     Initialise variables

      if (alpha.lt.1.0) then

         a = 1.0 - alpha
         p = a/(a + alpha*exp(-a))
         c = 1.0/alpha
         uf = p*(vsmall/a)**alpha
         vr = 1.0 - vsmall
         d = a*log(a)

      elseif (alpha.gt.1) then
 
         b = alpha - 1
         c = (12*alpha-3)/4

      endif

      do 10 i = 1,n

         if (alpha.eq.1) then

            call urand(1,u)
            x = -log(u)

         elseif (alpha.lt.1) then

20          call urand(1,r)
            if (r.ge.vr) then
               goto 20
            elseif (r.gt.p) then
               x = a - log((1.0 - r)/(1.0 - p))
               w = a*log(x)-d
            elseif (r.gt.uf) then
               x = a*(r/p)**c
               w = x
            else
               x = 0.0
               goto 40
            endif

            call urand(1,r)
            if (1.0-r.le.w.and.r.gt.0) then
               if (r*(w + 1.0).ge.1.0) then
                  goto 20 
               endif
               if (-log(r) <= w) then
                  goto 20 
               endif
            endif
             
         elseif (alpha.gt.1) then

            acc = 0
 30            call urand(1,u)
               call urand(1,v)
             
               w = u*(1-u)
               y = sqrt(c/w)*(u-0.5)
               x = b + y

               if (x.gt.0.0) then
                  z = 64*v*v*w*w*w
 
                  if (z.le.(1-2*y*y/x)) then
                     acc = 1
                  endif
 
                  if (acc.eq.0) then
                     if (2*(b*log(x/b)-y).ge.log(z)) then
                        acc = 1
                     endif
                  endif   
               endif  
         
            if (acc.eq.0) then
               goto 30
            endif
         endif 

 40      svec(i) = x
 10   continue

      return
      end

cccccccccc End of frgamma.f cccccccccc
