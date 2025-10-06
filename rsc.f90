!*******************************************************************************
! Reid93 soft core phenomenological potential, updated version,
! including one-pion-exchange with neutral-pion and charged-pion
! masses; coupling f^2=0.075. Tensor potential is regularized
! to equal zero at r=0
! Reference: V.G.J. Stoks et al., Phys. Rev. C 49, 2950 (1994).
!-------------------------------------------------------------------------------
! Fortran 90 version 2.1: August 2001
!
! E-mail: info@nn-online.org
! WWW: http://nn-online.org
!-------------------------------------------------------------------------------
! Integer variable 'wp' defines the precision of 'real' variables. It is by
! default put to double precision. This code is made with double-precision
! of 8 bytes in mind.  We do not know about the behavior of the code for
! other precisions, but do not expect anything anomalous.
!-------------------------------------------------------------------------------
! IN:  real(wp)          :: r         = distance in fermi
!      character(2)      :: type      = reaction, 'PP', 'NP', 'PN', or 'NN'
!      character(3)      :: pname     = name of the partial wave (see below),
!                                       maximum total angular momentum j = 9
! OUT: real(wp)          :: vpot(2,2) = potential matrix in MeV on LSJ-basis
!-------------------------------------------------------------------------------
! The variable 'pname' contains the name of the partial wave in spectral
! notation:  singlets:                 1S0  1P1  1D2  1F3  1G4 ...
!            triplets uncoupled:       3P0  3P1  3D2  3F3  3G4 ...
!            triplets coupled:              3C1  3C2  3C3  3C4 ...
! where 3C1 denotes   3S1 - epsilon_1 - 3D1 channel
!       3C2 denotes   3P2 - epsilon_2 - 3F2 channel
!       ...
!*******************************************************************************

      subroutine rreid93(r,pname,type,vpot)

      implicit none

      integer, parameter          :: wp = kind(1.0d0)

      real(wp), intent(in)        :: r
      character(2), intent(in)    :: type
      character(3), intent(in)    :: pname
      real(wp), intent(out)       :: vpot(2,2)
      real(wp)                :: ri,rj,vc,vt,vls,vtpi,f,x,mpi,vspi,vspis
      integer                     :: l,spin,j,nchan,iso
      character(3)                :: name
      character(*), parameter     :: large = 'CSPDFGHIKLMN'
      character(*), parameter     :: small = 'cspdfghiklmn'
      real(wp), parameter         :: r1 = 1.0_wp
      real(wp), parameter         :: r2 = 2.0_wp
      real(wp), parameter         :: r3 = 3.0_wp
      real(wp), parameter         :: r4 = 4.0_wp
      real(wp), parameter         :: r5 = 5.0_wp
      real(wp), parameter         :: r6 = 6.0_wp
      real(wp), parameter         :: r8 = 8.0_wp
      real(wp), parameter         :: f0pi = 0.075_wp
      real(wp), parameter         :: fcpi = 0.075_wp
      real(wp), parameter         :: hbc = 197.327053_wp
      real(wp), parameter         :: mpi0 = 134.9739_wp
      real(wp), parameter         :: mpic = 139.5675_wp
      real(wp), parameter         :: mpis = mpic
      real(wp), parameter    :: a(5,5) = reshape( (/                    &
     &  0.1756084e0_wp,-0.1414234e2_wp, 0.1518489e3_wp,-0.6868230e3_wp,  &
     &  0.1104157e4_wp,-0.4224976e2_wp, 0.2072246e3_wp,-0.3354364e3_wp,  &
     &  -0.1989250e1_wp,-0.6178469e2_wp, 0.2912845e2_wp, 0.1511690e3_wp,  &
     &  0.8151964e1_wp, 0.5832103e2_wp,-0.2074743e2_wp,-0.5840566e0_wp,  &
     &  -0.1029310e2_wp, 0.2263391e2_wp, 0.2316915e2_wp,-0.1959172e1_wp,  &
     &  -0.2608488e1_wp, 0.1090858e2_wp,-0.4374212e0_wp,-0.2148862e2_wp,  &
     &  -0.6584788e0_wp /) , (/5,5/) , order=(/2,1/) )
      real(wp), parameter         :: b(5,5) = reshape( (/                    &
     &  -0.2234989e2_wp, 0.2551761e3_wp,-0.1063549e4_wp, 0.1609196e4_wp,  &
     &  -0.3505968e1_wp,-0.4248612e1_wp,-0.5352001e1_wp, 0.1827642e3_wp,  &
     &  -0.3927086e3_wp, 0.5812273e2_wp,-0.2904577e1_wp, 0.3802497e2_wp,  &
     &   0.3395927e0_wp, 0.8318097e0_wp, 0.1923895e1_wp, 0.0913746e0_wp,  &
     &  -0.1274773e2_wp, 0.1458600e3_wp,-0.6432461e3_wp, 0.1022217e4_wp,  &
     &  -0.0461640e0_wp, 0.7950192e1_wp,-0.1925573e1_wp, 0.5066234e2_wp,  &
     &  0.8359896e1_wp /) , (/5,5/) , order=(/2,1/) )

!     Determine the quantumnumbers
      name = pname
      do l=1,12
          if (name(2:2) == small(l:l)) name(2:2) = large(l:l)
      enddo
      nchan = 1
      if (name(2:2) == 'C') nchan = 2
      if (name(1:1) == '1') spin = 0
      if (name(1:1) == '3') spin = 1
      read(name,'(2x,i1)') j
      l = j
      if (name == '3P0') l = 1
      if (nchan == 2) l = j - 1
      iso = mod(spin+l+1,2)

      ri = real(iso,kind=wp)
      rj = real(j,kind=wp)

!     OPE potential
      x     = mpi0/hbc * r
      f     = f0pi*(mpi0/mpis)**2*mpi0/r3
      vspis = f*yp(r1,r8,x)
      vspi  = f*y(r1,r8,x)
      vtpi  = f*z(r1,r8,x)
      select case(type)
      case('NP','PN','np','pn')
          x     = mpic/hbc * r
          f     = (r4*ri-r2)*fcpi*(mpic/mpis)**2*mpic/r3
          vspis = f*yp(r1,r8,x) - vspis
          vspi  = f*y(r1,r8,x)  - vspi
          vtpi  = f*z(r1,r8,x)  - vtpi
      endselect

!     Start clean
      vpot(:,:) = 0.0_wp

!     Non-OPE contribution to the potential
      mpi = (mpi0 + r2*mpic)/r3
      x = mpi/hbc * r

      select case(name)
      case('1S0')
          select case(type)
          case('PP','NN','pp','nn')
              vpot(1,1) =   a(1,1)*y(r2,r8,x)                              &
     &                     + a(1,2)*y(r3,r8,x)                           &
     &                     + a(1,3)*y(r4,r8,x)                           &
     &                     + a(1,4)*y(r5,r8,x)                           &
     &                     + a(1,5)*y(r6,r8,x)
          case('NP','PN','pn','np')
              vpot(1,1) =   b(1,1)*y(r3,r8,x)                           &
     &                     + b(1,2)*y(r4,r8,x)                           &
     &                     + b(1,3)*y(r5,r8,x)                           &
     &                     + b(1,4)*y(r6,r8,x)
          endselect
      case('1D2')
          vpot(1,1) =   a(2,1)*y(r4,r8,x)                                 &
     &                 + a(2,2)*y(r5,r8,x)                                 &
     &                 + a(2,3)*y(r6,r8,x)
      case('1G4')
          vpot(1,1) =   a(2,4)*y(r3,r8,x)
      case('3P0')
          vpot(1,1) =   a(3,1)*y(r3,r8,x)                               &
     &                 + a(3,2)*y(r5,r8,x)                               &
     &                 + a(2,5)*z(r3,r8,x)/r3
      case('3P1')
          vpot(1,1) =   a(3,3)*y(r3,r8,x)                               &
     &                 + a(3,4)*y(r5,r8,x)                               &
     &                 + a(3,5)*z(r3,r8,x)/r3
      case('3F3')
          vpot(1,1) =   a(4,5)*y(r3,r8,x)
      case('3C2','3C4')
          vc =  a(4,1)*y(r3,r8,x)                                       &
     &         + a(4,2)*y(r4,r8,x)                                        &
     &         + a(4,3)*y(r5,r8,x)                                       &
     &         + a(4,4)*y(r6,r8,x)
          vt = (   a(5,1)*z(r4,r8,x)                                    &
     &            + a(5,2)*z(r6,r8,x) )/r3
          select case(name)
          case('3C2')
              vls = a(5,3)*w(r3,r8,x) + a(5,4)*w(r5,r8,x)
          case('3C4')
              vls = a(5,5)*w(r3,r8,x)
          endselect
          vpot(1,1) = vc + (rj-r1)*vls - r2*(rj-r1)/(r2*rj+r1)*vt
          vpot(2,2) = vc - (rj+r2)*vls - r2*(rj+r2)/(r2*rj+r1)*vt
          vpot(1,2) = r6*sqrt(rj*(rj+r1))/(r2*rj+r1)*vt
          vpot(2,1) = vpot(1,2)
      case('1P1')
          vpot(1,1) =   b(2,1)*y(r3,r8,x)                               &
     &                 + b(2,2)*y(r4,r8,x)                                &
     &                 + b(2,3)*y(r5,r8,x)                               &
     &                 + b(2,4)*y(r6,r8,x)
      case('1F3')
          vpot(1,1) =   b(1,5)*y(r3,r8,x)                               &
     &                 + b(2,5)*y(r5,r8,x)
      case('3D2')
          vpot(1,1) =   b(3,1)*y(r3,r8,x)                               &
     &                 + b(3,2)*y(r5,r8,x)                               &
     &                 + b(3,3)*z(r3,r8,x)/r3
      case('3G4')
          vpot(1,1) =   b(3,4)*y(r3,r8,x)
      case('3C1','3C3')
          vc =  b(4,1)*y(r2,r8,x)                                       &
     &         + b(4,2)*y(r3,r8,x)                                       &
     &         + b(4,3)*y(r4,r8,x)                                       &
     &         + b(4,4)*y(r5,r8,x)                                       &
     &         + b(4,5)*y(r6,r8,x)
          vt = ( b(3,5)*z(r4,r8,x) + b(5,5)*z(r6,r8,x) )/r3
          select case(name)
          case('3C1')
              vls = b(5,1)*w(r3,r8,x) + b(5,2)*w(r5,r8,x)
          case('3C3')
              vls = b(5,3)*w(r3,r8,x) + b(5,4)*w(r5,r8,x)
          endselect
          vpot(1,1) = vc + (rj-r1)*vls - r2*(rj-r1)/(r2*rj+r1)*vt
          vpot(2,2) = vc - (rj+r2)*vls - r2*(rj+r2)/(r2*rj+r1)*vt
          vpot(1,2) = r6*sqrt(rj*(rj+r1))/(r2*rj+r1)*vt
          vpot(2,1) = vpot(1,2)
      case default
!         j > 5
          select case(spin)
          case(0)
              select case(iso)
              case(1)
                  vpot(1,1) =  a(1,1)*y(r2,r8,x)                        &
     &                        + a(1,2)*y(r3,r8,x)                        &
     &                        + a(1,3)*y(r4,r8,x)                        &
     &                        + a(1,4)*y(r5,r8,x)                        &
     &                        + a(1,5)*y(r6,r8,x)
              case(0)
                  vpot(1,1) =  b(2,1)*y(r3,r8,x)                        &
     &                        + b(2,2)*y(r4,r8,x)                         &
     &                        + b(2,3)*y(r5,r8,x)                         &
     &                        + b(2,4)*y(r6,r8,x)
              endselect
          case(1)
              select case(iso)
              case(1)
                  vc =  a(4,1)*y(r3,r8,x)                               &
     &                 + a(4,2)*y(r4,r8,x)                               &
     &                 + a(4,3)*y(r5,r8,x)                               &
     &                 + a(4,4)*y(r6,r8,x)
                  vt = ( a(5,1)*z(r4,r8,x) + a(5,2)*z(r6,r8,x) )/r3
              case(0)
                  vc =  b(4,1)*y(r2,r8,x)                               &
     &                 + b(4,2)*y(r3,r8,x)                               &
     &                 + b(4,3)*y(r4,r8,x)                               &
     &                 + b(4,4)*y(r5,r8,x)                               &
     &                 + b(4,5)*y(r6,r8,x)
                  vt =  ( b(3,5)*z(r4,r8,x) + b(5,5)*z(r6,r8,x) )/r3
              endselect
              select case(nchan)
              case(1)
               if (l == (j-1)) vpot(1,1) = vc - r2*(rj-r1)/(r2*rj+r1)*vt
               if (l == j)     vpot(1,1) = vc + r2*vt
               if (l == (j+1)) vpot(1,1) = vc - r2*(rj+r2)/(r2*rj+r1)*vt
              case(2)
                  vpot(1,1) = vc - r2*(rj-r1)/(r2*rj+r1)*vt
                  vpot(2,2) = vc - r2*(rj+r2)/(r2*rj+r1)*vt
                  vpot(1,2) = r6*sqrt(rj*(rj+r1))/(r2*rj+r1)*vt
                  vpot(2,1) = vpot(1,2)
              endselect
          endselect
      endselect
      vpot(:,:) = mpi * vpot(:,:)

!     All together
      select case(nchan)
      case(1)
          if (spin == 0) then
              if (l == 0) then
                  vpot(1,1) = vpot(1,1) - r3*vspis
              else
                  vpot(1,1) = vpot(1,1) - r3*vspi
              endif
          elseif (l == j) then
              vpot(1,1) = vpot(1,1) + vspi + r2*vtpi
          elseif (name == '3P0') then
              vpot(1,1) = vpot(1,1) + vspi - r2*(rj+r2)/(r2*rj+r1)*vtpi
          endif
      case(2)
          if (l == 0) then
              vpot(1,1) = vpot(1,1) + vspis - r2*(rj-r1)/(r2*rj+r1)*vtpi
          else
              vpot(1,1) = vpot(1,1) + vspi - r2*(rj-r1)/(r2*rj+r1)*vtpi
          endif
          vpot(2,2) = vpot(2,2) + vspi - r2*(rj+r2)/(r2*rj+r1)*vtpi
          vpot(1,2) = vpot(1,2) + r6*sqrt(rj*(rj+r1))/(r2*rj+r1)*vtpi
          vpot(2,1) = vpot(1,2)
      endselect

      contains

!*******************************************************************************

      function y(n,m,x)

      implicit none
      real(wp), intent(in) :: n,m
      real(wp), intent(in) :: x
      real(wp)             :: y
      real(wp), parameter  :: r1 = 1.0_wp
      real(wp), parameter  :: r2 = 2.0_wp


      if (x < 1.0e-4_wp) then
          y = - n + m/r2 + n**2/(r2*m)
      else
          y = exp(-n*x)/x - exp(-m*x)/x*(r1+(m**2-n**2)*x/(r2*m))
      endif

      end function y

!*******************************************************************************

      function yp(n,m,x)

      implicit none
      real(wp), intent(in) :: n,m
      real(wp), intent(in) :: x
      real(wp)             :: yp
      real(wp), parameter  :: r1 = 1.0_wp
      real(wp), parameter  :: r2 = 2.0_wp
      real(wp)             :: n2,m2,d


      n2 = n**2
      m2 = m**2

      if (x < 1.0e-4_wp) then
          d  = m*(m2-n2)/(r2*n2)
          yp = - n + m - d + x*(n2/r2 + m2/r2 + m*d)
      else
          yp = exp(-n*x)/x - exp(-m*x)/x*(r1+(m2-n2)*m*x/(r2*n2))
      endif

      end function yp

!*******************************************************************************

      function w(n,m,x)

      implicit none
      real(wp), intent(in) :: n,m
      real(wp), intent(in) :: x
      real(wp)             :: w
      real(wp)             :: n2,m2,x2
      real(wp), parameter  :: r1 = 1.0_wp
      real(wp), parameter  :: r2 = 2.0_wp
      real(wp), parameter  :: r3 = 3.0_wp
      real(wp), parameter  :: r6 = 6.0_wp
      real(wp), parameter  :: r8 = 8.0_wp


      n2 = n**2
      m2 = m**2

      if (x < 1.0e-4_wp) then
          w = (r2*n - r3*m + m2*m/n2)/r6 + x*(r2*m2 - n2 - m2*m2/n2)/r8
      else
          x2 = x*x
          w =  exp(-n*x)/x*(r1/(n*x)+r1/(n2*x2))                         &
     &        - exp(-m*x)/x*(r1/(m*x)+r1/(m2*x2))*m2/n2                  &
     &        - exp(-m*x)/x*(m2-n2)/(r2*n2)
      endif

      end function w

!*******************************************************************************

      function z(n,m,x)

      implicit none
      real(wp), intent(in) :: n,m
      real(wp), intent(in) :: x
      real(wp)             :: z
      real(wp)             :: n2,m2,x2
      real(wp), parameter  :: r1 = 1.0_wp
      real(wp), parameter  :: r2 = 2.0_wp
      real(wp), parameter  :: r3 = 3.0_wp
      real(wp), parameter  :: r4 = 4.0_wp
      real(wp), parameter  :: r8 = 8.0_wp


      n2 = n**2
      m2 = m**2

      if (x < 1.0e-4_wp) then
          z  = x*(n2 + r3*m2*m2/n2 - r4*m2)/r8
      else
          x2 = x*x
          z  =  exp(-n*x)/x*(r1+r3/(n*x)+r3/(n2*x2))                    &
     &         - exp(-m*x)/x*(r1+r3/(m*x)+r3/(m2*x2))*m2/n2               &
     &         - exp(-m*x)*(r1+r1/(m*x))*m*(m2-n2)/(r2*n2)
      endif

      end function z

!*******************************************************************************

      end subroutine rreid93
