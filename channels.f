      module cdccchannels
      implicit none

       integer :: indexthoxmax,exmax
       integer,dimension(:,:),allocatable :: indexthox



!Jacobi set for 3b channels
       Type nch_cdcc
!  | (l s_2b)j_2b, (lam j? ) j_spect; J >
        integer :: nchmax
        integer,dimension(:),allocatable :: l,lam
        real*8,dimension(:),allocatable :: s2b,j2b
        real*8,dimension(:),allocatable  ::j,j_spect
        integer,dimension(:),allocatable :: n
       End type



       type(nch_cdcc) :: incdcc

       contains


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine alpha_cdcc_in()
c     index of  | \alpha >_{in}
c     index of {(alpha_2b (lam jt)j_spect ; J}
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: l,lambda,l2bmax,lmax_cdcc
      real*8 :: jmin,jmax
      real*8 :: s2b,j2b,j_spect,J,s
      integer :: nJ,nJ2b,nJ_spect
      integer :: nch !index of {(l(jxjb)s2b)J2b (lam jt)J_spect ; J}
      real*8 :: jtmin, jtmax
      integer,dimension(10) :: l2b
      integer :: nl
      namelist /cdccwf/ l2b,lmax_cdcc,jtmin, jtmax
      l2b=-99
      open (unit=551,file='iav.in')
      read(551,nml=cdccwf)

      s2b=0.0d0 ! by assuming all spin zero particle, need to improved
      s=0.0d0 ! spin of the target


      incdcc%nchmax=0
C     do l=0,l2bmax
      do nl=1,10
        l=l2b(nl)
        if(l<0) cycle
        do nJ2b=nint(2.*(l+s2b)),nint(2.*abs(l-s2b)),2
           j2b=nJ2b/2.0d0
         do lambda=0,lmax_cdcc
           do nJ_spect=nint(2*abs(lambda-s)),nint(2*(lambda+s)),2
            J_spect=nJ_spect/2.0d0
            jmin=max(jtmin,abs(J2b-J_spect))
            jmax=min(jtmax,J2b+J_spect)
      	    do nj=nint(2.*jmin),nint(2.*jmax),2
      	      incdcc%nchmax=incdcc%nchmax+1
      	      write(0,*)'incdcc%nchmax=',incdcc%nchmax
      	    end do
           end do
         end do
        end do
      end do
      write(8,10)incdcc%nchmax
10    format('there are',I3,1X,'cdcc channels')

       if (allocated(incdcc%l)) deallocate(incdcc%l)
       if (allocated(incdcc%lam)) deallocate(incdcc%lam)
       if (allocated(incdcc%J)) deallocate(incdcc%J)
       if (allocated(incdcc%s2b)) deallocate(incdcc%s2b)
       if (allocated(incdcc%J2b)) deallocate(incdcc%J2b)
       if (allocated(incdcc%j_spect)) deallocate(incdcc%j_spect)
       if (allocated(incdcc%n)) deallocate(incdcc%n)



       allocate(incdcc%l(1:incdcc%nchmax),incdcc%lam(1:incdcc%nchmax))
       allocate(incdcc%J(1:incdcc%nchmax),incdcc%s2b(1:incdcc%nchmax))
       allocate(incdcc%J2b(1:incdcc%nchmax))
       allocate(incdcc%j_spect(1:incdcc%nchmax))
       allocate(incdcc%n(1:incdcc%nchmax))
       incdcc%n=0

      nch=1
      write(8,20)
      write(8,30)
C     do l=0,l2bmax
      do nl=1,10
        l=l2b(nl)
        if(l<0) cycle
        do nJ2b=nint(2.*(l+s2b)),nint(2.*abs(l-s2b)),2
           j2b=nJ2b/2.0d0
         do lambda=0,lmax_cdcc
           do nJ_spect=nint(2*abs(lambda-s)),nint(2*(lambda+s)),2
            J_spect=nJ_spect/2.0d0
            jmin=max(jtmin,abs(J2b-J_spect))
            jmax=min(jtmax,J2b+J_spect)
      	    do nj=nint(2.*jmin),nint(2.*jmax),2
               j=nj/2.0d0
      	       incdcc%l(nch)=l
               incdcc%lam(nch)=lambda
               incdcc%J(nch)=J
               incdcc%s2b(nch)=s2b
               incdcc%J2b(nch)=J2b
               incdcc%j_spect(nch)=j_spect

               write(8,40)nch,0,l,s,s,s2b,J2b,
     +                       lambda,s,j_spect,j
      	           nch=nch+1
      	    end do
         end do
      end do
      end do
      end do

20    format('---For CDCC channels the coupling coefficients are')
30    format(' a3b','|', 'a2b','|','(',' l ','(',' jb ',' jx ',')',
     +      ' s2b ',')',' J2b ', '(',' lam ',' jt ',')',' J3 ',',',
     +       ' Jtot ')
40    format(I4,I4,1x,I3,2x,f3.1,2x,f3.1,2x,f3.1,2x,
     +       f4.1,3x,I3,1x,f3.1,1x,f4.1,3x,f4.1)
      end subroutine
c-----------------------------------------------------------------------
      subroutine recoupling()
        use channels, only: jpsets,jpiset
        implicit none
        integer :: ne
        integer :: jset,ie,acdcc

        exmax=0
        do jset=1,jpsets
          ne=jpiset(jset)%nex
          exmax=max(exmax,ne)
        end do

         allocate(indexthox(1:exmax,1:incdcc%nchmax))
         indexthoxmax=0
        do ie=1, exmax
          do acdcc=1, incdcc%nchmax
             indexthoxmax=indexthoxmax+1
             indexthox(ie,acdcc)=indexthoxmax
          end do
        end do












      end subroutine





      end module
