c        ********************** absorplh**********************
c        *                      -----                       *
c        *  this subroutine calculates the imaginary part   *
c        *  of refractive index using formula from:   	    *
c        *  P.Bonoli, Linear theory of lower hybrid heating *
c        *  IEEE Transaction on Plasma Science,Vol. PS-12,  *
c        *  No.2,(1984), p 95-107                           *
c        *  See also Paul T. Bonoli and Ronald C. Englade,  *
c        *  Phys. Fluids 29, 2937 (1986).                   *
c        ****************************************************
c         input parameters: u(6)  -z,r,phi,n_z,n_r,m  
c                           cnpar -N_par			  
c                           cnper -N_pre (real part)   
c                           tempe -electron temperature  (keV)  
c                           dense -electron density   (10*13/cm**3
c                           tempiar -ions temperature  (keV)      *
c                           b_z,b_r,b_phi,bmod-magnetic field	  *
c                           nbulk - number of plasma components	  *
c                           frqncy- wave friquency [GHz]	  *
c                           zeff  - plasma charge      		  *
c         output parameter: cnprim_e -imaginary part of N_perp    *
c                             (uncollisional electron absorption) *
c                           cnprim_i -imaginary part of N_perp    *
c                                (uncollisional ion absorption)   *
c                           cnprim_cl-imagionary part of N_perp	  *
c                                collisional(ions and electron)	  *
c-----------------------------------------------------------------*
      subroutine absorplh(u,cnpar,cnper,tempe,dense,tempiar
     1 ,ifsd,Lambda0,DLambda,b_z,b_r,b_phi,nbulk,bmod,frqncy,zeff,
     1 cnprim_e,cnprim_i,cnprim_cl,cnprim_s)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      dimension u(6),deruu(6)
      dimension tempiar(nbulka)
      real*8 ifsd(nbulka),Lambda0(nbulka),DLambda(nbulka)
      dimension cnprim_s(nbulka),di_is(nbulka)
      real*8 BLL,Xresn,Yresn,Xres,Yres
********************************************************************
c electron absorption (Landau damping, formula (21a) in the Bonoli article)
c Di_e=2*sqrt(pi)*(omegape/omega)**2*(N_par*N_perp)**2*abs(x_oe**3)
c      *exp(-x_oe**2)
c--------------------------------------------------------------------
      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)
c--------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      spi=dsqrt(pi)
      cnpar2=cnpar*cnpar
      cnpar4=cnpar2*cnpar2
      cnper2=cnper*cnper
      cnper4=cnper2*cnper2
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c----------------------------------------
c     ve=(sqrt(Te/me)),    ve in (cm/sec),Te(keV)
      ve=1.87d9*dsqrt(0.5d0*tempe)
      sqrt2=dsqrt(2.d0)
c----------------------------------------
c     x_oe=omega/((sqrt(2))*k_par*ve)=cvac/(sqrt(2)*N_par*ve)
      x_oe=dabs(cvac/(sqrt2*cnpar*ve))
      x_oe2=x_oe*x_oe
      x_oe3=x_oe2*x_oe
      ye=y(z,r,phi,1)
      xe=x(z,r,phi,1)
c--------------------------------------
c     di_e is imaginary part of the dispersion function(uncollisional)
      di_e=2.d0*spi*xe*cnper2*cnpar2*x_oe3*dexp(-x_oe2)
c*********************************************************************
c     di_i is imaginary part of the dispersion function(uncollisional),
c     formula (21b) in the Bonoli article.
c----------------------------------------------------------
      di_i=0.d0
      do i=1,nbulka
            di_is(i) = 0.d0
      enddo
      di_is(1) = di_e

c      write(*,*)' absorplh.f nbulk',nbulk
      vc = 0.d0
      do i=2,nbulk
            if (ifsd(i) < 0.5d0) then
                  tempi=tempiar(i)
c        write(*,*)'i,tempi,dmas(i)',i,tempi,dmas(i)
c-------vi=(sqrt(Ti/mi)),    vi in (cm/sec),Ti(keV)
                  vi=1.87d9*dsqrt(0.5d0*tempi/dmas(i))
                  
c-------x_oi=omega/((sqrt(2))*k_par*vi)=cvac/(sqrt(2)*N_par*vi)
                  x_oi=dabs(cvac/(sqrt2*cnper*vi))
                  x_oi2=x_oi*x_oi
                  x_oi3=x_oi2*x_oi
                  x_oi4=x_oi2*x_oi2
                  yi=y(z,r,phi,i) ! yi = omega_cyclotron_i/omega
                  xi=x(z,r,phi,i) ! xi = (omega_pl_i/omega)**2

                  sum=0.d0
                  
c        write(*,*)'absorplh.f i,spi,xi,cnper4,vi',i,spi,xi,cnper4,vi
c        write(*,*)'absorplh.f x_oi,x_oi2,x_oi3,x_oi4',
c     &x_oi,x_oi2,x_oi3,x_oi4
                  


                  pdi_i=2.d0*spi*xi*cnper4*x_oi3*dexp(-x_oi2)
c        write(*,*)'pdi_i',pdi_i
                  psum=yi*dabs(cvac/(sqrt2*cnpar*vi))/spi
                  if (pdi_i.lt.1.d-10) goto 20
                  n0=1/yi
                  pn=10.d0*sqrt2*vi*dabs(cnpar)/cvac
                  x_ni0=(1.d0-n0*yi)*cvac/(sqrt2*cnpar*vi)
                  n1=(1-pn)/yi-1
                  x_ni1=(1.d0-n1*yi)*cvac/(sqrt2*cnpar*vi)
                  n2=(1+pn)/yi+1
                  x_ni2=(1.d0-n2*yi)*cvac/(sqrt2*cnpar*vi)
                  
c        write(*,*)'absorplh.f n1,n2',n1,n2 
                  do n=n1,n2
c--------------------------------------
c         x_ni=(omega-n*omegas_ci)/sqrt(2)*k_par*vi)
                  x_ni=(1.d0-n*yi)*cvac/(sqrt2*cnpar*vi)
                  x_ni2=x_ni*x_ni
                  if (x_ni2.gt.100)then
                        goto 10
                  endif
                  sum=sum+dexp(-x_ni2)
10                continue
                  enddo !n
c        write(*,*)'absorplh.f,pdi_i,sum',pdi_i,sum
20                sum=sum*psum
                  di_is(i) = pdi_i*sum
                  di_i=di_i + di_is(i)
                  vc = vc + xi
            endif
      enddo !i

      vc = 0.75d0*spi*vc/xe
      vc3 = (ve*dsqrt(2.d0))**3*vc
      vc = vc3**(1.d0/3.d0)
      do i=2,nbulk
            if (ifsd(i) >= 0.5d0) then
                  ! ! Unmagnetized ImD
                  ! ! Calculate w/k_perp/v_0 and v_c/v_0
                  ! tempi = tempiar(i)
                  ! vi = 1.87d9*dsqrt(tempi/dmas(i))
                  ! x_oi = dabs(cvac/(cnper*vi))
                  ! if (x_oi > 1.d0) goto 30
                  ! x_uc = vc/vi
                  ! ! Check unisotropic distribution
                  ! if (abs(gImD(i) - 1.d0) > 1.d-5) then
                  !       x_oi = x_oi*gzeta(i)
                  !       x_uc = x_uc*guc(i)
                  ! endif
                  ! ! Calculate ImD
                  ! x_oi3 = x_oi*x_oi*x_oi
                  ! x_uc3 = x_uc*x_uc*x_uc
                  ! xi=x(z,r,phi,i)
                  ! fac = 1.5d0*pi/dlog(1+1/x_uc3)
                  ! psum = x_oi3*(1/(x_oi3+x_uc3)-1/(1+x_uc3))
                  ! psum = psum*fac*cnper4*xi*gzeta(i)
                  ! ! Accumulate ImD of ions
                  ! if (abs(gImD(i) - 1.d0) > 1.d-5) then
                  !       psum = psum*gImD(i)
                  ! endif
                  ! di_is(i) = psum
                  ! di_i = di_i + psum




                  ! Calculate pinch angle normalization parameter
                  if (dabs(DLambda(i)).lt.1.d-5) then
                        BLL = 2
                  else
                        call Pinch_BLL(Lambda0(i),DLambda(i),BLL)
                  endif
                  ! Calculate w/k_perp/v_0 and v_c/v_0
                  yi=y(z,r,phi,i)                           ! yi = omega_cyclotron_i/omega
                  Yharm = 1/yi
                  xi=x(z,r,phi,i)                           ! xi = (omega_pl_i/omega)**2
                  tempi = tempiar(i)
                  vi = 1.87d9*dsqrt(tempi/dmas(i))          ! v_0
                  x_oi = dabs(cvac/(cnper*vi))
                  x_uc = vc/vi                              ! u_c
                  x_uc3 = x_uc*x_uc*x_uc
                  bpar = dabs(cnpar) * Yharm * vi/cvac      ! k_para*v0/wc
                  bper = dabs(cnper) * Yharm * vi/cvac      ! k_perp*v0/wc
                  ! Calculate general parameters
                  dharmn = Yharm - bpar
                  iharmn = idint(dharmn)              ! ceil(w-kpar*vpar)/wc
                  if (dble(iharmn).lt.dharmn) iharmn = iharmn + 1
                  dharmp = Yharm + bpar
                  iharmp = idint(dharmp)              ! floor(w+kpar*vpar)/wc
                  if (dble(iharmp).gt.dharmp) iharmp = iharmp - 1
                  Xres = 0
                  Yres = 0
                  do j=iharmn,iharmp
                        unpar = (Yharm-j)/dmax1(bpar,1.d-5)
                        call Pinch_Slowing_Down_Integral(unpar,
     & x_uc3,j,Yharm,Lambda0(i),DLambda(i),bpar,bper,Xresn,Yresn)
                        Xres = Xres + Xresn
                        Yres = Yres + Yresn
                  enddo
                  psum = 3*Yharm*Xres + 4*Yres
                  psum = psum/dmax1(bpar,1.d-5)*x_oi**2
                  psum = psum*3*pi/dlog(1+1/x_uc3)/BLL
                  psum = psum*(cnpar2+cnper2)*cnper2*xi
                  di_is(i) = psum
                  di_i = di_i + psum
            endif
30          continue
      enddo

c*********************************************************************
c     di_ic is the imaginary part of the dispersion function(collisional)
c     Collisional, ions and electron damping coefficients: (cgs units)
c     formula (26) from the Bonoli article
c     1.47e-9 = (4/3)*sqrt(2*pi)*16*q**4/(sqrt(me)*(1.6e-9 erg/keV)**1.5)
      rnuei = 1.47d-9*1.d+13*dense*zeff/tempe**1.5 ! [1/sec]
      omega=2.d0*pi*frqncy*1.d9	  ! [1/sec], frqncy in GHz
      di_ic=rnuei*(xe*cnper2/(ye*ye)+xe*cnpar2)*cnper2/omega
c--------------------------------------------------------------
c     calculation of (dD_real/dN_par)
c     D is from the P.Bonoli article  , formula (7),(9) and (17) 
      ! xi=x(z,r,phi,2)
      xi = 0
      xiT = 0
      do i=2,nbulk
            xi0 = x(z,r,phi,i)
            xi = xi + xi0
            if (ifsd(i) >= 0.5d0) then
                  tempi = tempiar(2)
                  vi2 = 1.87d9*1.87d9*0.5d0*tempi/dmas(i)
                  xiT = xiT + xi0*vi2/cvac/cvac
            endif
      enddo
      xiT = 1.5d0*xiT
      epsper=1.d0+xe/(ye*ye)-xi
      epspar=1.d0-xe-xi
      epsxy=xe/ye
c------------------------------------------
c p6 calculations
      ye2=ye*ye
      ye4=ye2*ye2
      pe=ve/cvac
      pe2=pe*pe
c     vi=(sqrt(Ti/mi)),    vi in (cm/sec),Ti(keV)
      ! tempi=tempiar(2)
      ! vi=1.87d9*dsqrt(0.5d0*tempi/dmas(2))
      ! xi=x(z,r,phi,2)
      ! pi=vi/cvac
      ! pi2=pi*pi
      ! p6=-(3.d0*xi*pi2+0.75d0*xe/ye4*pe2)
      p6=-(3.d0*xiT+0.75d0*xe/ye4*pe2)
c------------------------------------------
      p4=epsper
      p2=(epsper+epspar)*(cnpar2-epsper)+epsxy*epsxy
      p0=epspar*((cnpar2-epsper)*(cnpar2-epsper)-epsxy*epsxy)
c------------------------------------------

      dddnper=(6.d0*p6*cnper4+4.d0*p4*cnper2+2.d0*p2)*cnper
      cnprim_e=-(di_e/dddnper)
      cnprim_i=-(di_i/dddnper)
      cnprim_cl=-(di_ic/dddnper) !   *coll_mult !BH: added 191017, Needs passing?
      !YuP: coll_mult is applied after call of this subroutine [in prep3d]
      do i=1,nbulk
            cnprim_s(i) = -(di_is(i)/dddnper)
      enddo
      cnprim_e=dabs(cnprim_e)
      cnprim_i=dabs(cnprim_i)
      cnprim_cl=dabs(cnprim_cl)
      do i=1,nbulk
            cnprim_s(i) = dabs(cnprim_s(i))
      enddo

c      write(*,*)'lh cnprim_e,cnprim_i,cnprim_cl',
c     1 cnprim_e,cnprim_i,cnprim_cl

      return
      end


cSZYang0392 2026.3.21
! Calculate the pinch-angle normalization factor
      subroutine Pinch_BLL(L0,DL,BLL)
      implicit none
      real*8 L0,DL,pi
      real*8 h,x,y,BLL
      integer i,imax

      pi=4.d0*datan(1.d0)
      h = 0.01
      imax = dint(pi/h)
      if (h*dble(imax).gt.pi) imax = imax - 1

      BLL = 0
      do i=0,imax-1
            x = i*h
            y = dexp(-((dsin(x)**2-L0)/DL)**2)*dsin(x)
            BLL = BLL + y
      enddo
      BLL = BLL*h
      
      return
      end



cSZYang0392 2026.3.21
! Calculate the integral of pinch angle slowing-down damping
! with Simpson method
      subroutine Pinch_Slowing_Down_Integral(unpar,
     & uc3,nharm,Y,L0,DL,bpar,bper,Xres,Yres)
      implicit none
      include 'param.i'
      include 'write.i'
      integer nharm,i,mn,mnpar
      real*8 uper(nuperint),Xint(nuperint),Yint(nuperint)
      real*8 unpar, uc3, Y, L0, DL, bpar, bper
      real*8 a,b,Xa,Xb,Ya,Yb,h,H0x,H0y
      real*8 Snx,Sny,S2nx,S2ny,Hnx,Hny,H2nx,H2ny,Rnx,Rny
      real*8 Reltol
      real*8 Xres,Yres
      logical ifOKx,ifOKy

      ! Set the relative error
      Reltol = 5.e-3
      ifOKx = .false.
      ifOKy = .false.

      ! Integrand value at the start and end points
      a = 0.d0
      b = dsqrt(dmax1(1.d0-unpar**2,0.d0))
      uper(1) = a
      call Pinch_Slowing_Down_Integrand(uper(1),1,unpar,uc3,nharm,
     & Y,L0,DL,bpar,bper,Xint(1),Yint(1))
      Xa = Xint(1)
      Ya = Yint(1)
      uper(1) = b
      call Pinch_Slowing_Down_Integrand(uper(1),1,unpar,uc3,nharm,
     & Y,L0,DL,bpar,bper,Xint(1),Yint(1))
      Xb = Xint(1)
      Yb = Yint(1)

      ! Initialize the integral
      mn = 1000
      mnpar = idint(b*1.0d3)
      if (mnpar.lt.1) mnpar = 1
      if (mn>mnpar) mn = mnpar
      h = (b-a)/2/mn
      Snx = Xa + Xb
      Sny = Ya + Yb
      ! Part 1 : fx_2i
      do i=1,mn-1
            uper(i) = a + 2*i*h
      enddo
      call Pinch_Slowing_Down_Integrand(uper(1),mn-1,unpar,
     & uc3,nharm,Y,L0,DL,bpar,bper,Xint(1),Yint(1))
      do i=1,mn-1
            Snx = Snx + 2*Xint(i)
            Sny = Sny + 2*Yint(i)
      enddo
      ! Part 2 : fx_{2i+1}
      do i = 0,mn-1
            uper(mn+i) = a + (2*i+1)*h
      enddo
      call Pinch_Slowing_Down_Integrand(uper(mn),mn,unpar,
     & uc3,nharm,Y,L0,DL,bpar,bper,Xint(mn),Yint(mn))
      H0x = 0
      H0y = 0
      do i=0,mn-1
            H0x = H0x + Xint(mn+i)
            H0y = H0y + Yint(mn+i)
      enddo
      Snx = Snx + 4.d0*H0x
      Sny = Sny + 4.d0*H0y
      ! Multiply h
      Snx = Snx * (h/3)
      Sny = Sny * (h/3)
      S2nx = Snx
      S2ny = Sny
      
      ! Loop to find the accurate value
      do while ((4*mn-1) .lt. nuperint)
            ! Calculate Hnx,Hny
            Hnx = H0x * (h/3)
            Hny = H0y * (h/3)
            ! Set new mesh
            h = h/2
            mn = 2*mn
            do i=0,mn-1
                  uper(mn+i) = a + (2*i+1)*h
            enddo
            ! Calculate H2nx, H2ny
            call Pinch_Slowing_Down_Integrand(uper(mn),
     & mn,unpar,uc3,nharm,Y,L0,DL,bpar,bper,
     & Xint(mn),Yint(mn))
            H0x = 0
            H0y = 0
            do i=0,mn-1
                  H0x = H0x + Xint(mn+i)
                  H0y = H0y + Yint(mn+i)
            enddo
            H2nx = H0x * (h/3)
            H2ny = H0y * (h/3)
            ! Calculate S2nx
            S2nx = Snx/2 + 4*H2nx - Hnx
            S2ny = Sny/2 + 4*H2ny - Hny
            ! Calculate the error and update data
            Rnx = (S2nx - Snx)/15
            Rny = (S2ny - Sny)/15
            ifOKx = dabs(Rnx).lt.Reltol*dmax1(dabs(S2nx),1.d-30)
            ifOKy = dabs(Rny).lt.Reltol*dmax1(dabs(S2ny),1.d-30)
            if (ifOKx.and.ifOKy) exit
            Snx = S2nx
            Sny = S2ny
      end do

      if (.not.(ifOKx.and.ifOKy)) then
            iintgr_err_one_ray = 1
      endif

      Xres = S2nx
      Yres = S2ny

      return
      end



cSZYang0392 2026.3.21
! Calculate the integrand of pinch angle slowing-down damping
      subroutine Pinch_Slowing_Down_Integrand(uper,nuper,unpar,
     & uc3,nharm,Y,L0,DL,bpar,bper,Xint,Yint)
      implicit none
      integer nuper, nharm, i
      real*8 uper(*), Xint(*), Yint(*)
      real*8 unpar, uc3, Y, L0, DL, bpar, bper
      real*8 uper2,u2,u,fsd,xpinch1,xpinch2,fpinch,Jnb,Jnbp,fcommon,krho

      do i=1,nuper
         uper2 = uper(i)**2
         u2 = uper2 + unpar**2
         u = dsqrt(u2)
         fsd = 1/max(u**3 + uc3, 1.d-5)
         xpinch1 = uper2/max(u2,1.d-5) - L0
         krho = bper*uper(i)
         call bes_calc(krho, nharm, Jnb, Jnbp)
         if (dabs(DL).gt.1.d-5) then
            xpinch2 = xpinch1/max(dabs(DL),1.d-5)
            fpinch = exp(-xpinch2**2)
            fcommon = uper(i)*fsd*fpinch*Jnb**2
            Yint(i) = fcommon*xpinch2/max(dabs(DL),1.d-5)*unpar
            Yint(i) = Yint(i)/max(u2**2,1.d-5)
            Yint(i) = Yint(i)*(Y*unpar - bpar*u2)
         else
            fcommon = uper(i)*fsd*Jnb**2
            Yint(i) = 0.d0
         endif
         Xint(i) = fcommon*u*fsd
      enddo

      return
      end

      subroutine absorp_collisional(tempe,dense,frqncy,zeff,
     & v_gr_perp,coll_mult,
     & cnprim_cl)
c      calculates the imaginary part of refractive index
c      due to collisional electron-ion damping 
c      using the formula 
c      Im(N)= nu_ei/(omega*Vgr_perp/clight)
c
c      input parameters: 
c                        tempe -electron temperature  (keV)  
c                        dense -electron density   (10*13/cm**3	  
c                        frqncy- wave friquency [GHz]	  
c                        zeff  - plasma charge
c        		 v_gr_perp -perpendicular to the magnetic field
c                           group velocity normalized by clight
c      output parameter: 
c                        cnprim_cl-imaginary part of N_perp	  
c                                collisional(ions and electron)	  
      implicit none

c-----input
      real*8 tempe,dense,frqncy,zeff,v_gr_perp,coll_mult

c-----output
      real*8 cnprim_cl

c-----local
      real*8 arg,cln,r_nu_ei,pi,omega
      
c-----Coulomb Logarithm = cln      
      arg=(1.d3/tempe)*dsqrt(10.d0*dense) !tempe KeV, 
                                           !dense 10**13/sm**3   
      cln= 24.d0 -dlog(arg)

c-----electron-ion collision rate =r_nu_ei [1/sec] 
      r_nu_ei=dense/(tempe*dsqrt(tempe))*zeff*cln*0.919d3 ![1/sec]

c      write(*,*)'dense,tempe,zeff,cln,r_nu_ei',
c     &dense,tempe,zeff,cln,r_nu_ei

c-----wave angle frequency 
      pi=4.d0*datan(1.d0)
      omega=2.d0*pi*frqncy*1.d9	  ! [1/sec], frqncy in GHz

      cnprim_cl=dabs(r_nu_ei/v_gr_perp/omega)*coll_mult
c      write(*,*)'cln,r_nu_ei,v_gr_perp,omega,cnprim_cl',
c     & cln,r_nu_ei,v_gr_perp,omega,cnprim_cl
      
      return
      end
