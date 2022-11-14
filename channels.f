      module cdccchannels
      implicit none

       integer :: indexthoxmax,exmax
       integer,dimension(:,:),allocatable :: indexthox
       
       
       integer,dimension(:,:),allocatable :: abar 



!Jacobi set for 3b channels
       Type nch_cdcc
!  | (l s_2b)j_2b, (lam j? ) j_spect; J >
        integer :: nchmax
        integer,dimension(:),allocatable :: l,lam
        real*8,dimension(:),allocatable :: s2b,j2b
        real*8,dimension(:),allocatable  ::j,j_spect
        integer,dimension(:),allocatable :: n
       End type
       
       Type n2ch_index 
          integer,dimension(:),allocatable :: a3b, ain 
          integer :: max
       End type


       type(n2ch_index) :: n2ch
       type(nch_cdcc) :: incdcc,alpha_in

       contains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine abar_cdcc()
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none 
      integer :: i2ch ! coupling index for alpha_0 and alpha_cdcc 
      integer :: a_in, a_cdcc 
      integer :: l3b, lam3b, lin, lamin, nj_3b, nj_in
      
      n2ch%max = 0 
      if(allocated(abar)) deallocate(abar)
      allocate(abar(1:alpha_in%nchmax,1:incdcc%nchmax))
      abar=0 
      
      do a_cdcc=1, incdcc%nchmax
         nj_3b = nint(2. * incdcc%J(a_cdcc))
         l3b = incdcc%l(a_cdcc)
         lam3b = incdcc%lam(a_cdcc)
         
         do a_in=1, alpha_in%nchmax 
            nj_in = nint(2. * alpha_in%J(a_in))
            if(nj_3b /= nj_in) cycle 
            lin =  alpha_in%l(a_in)
            lamin = alpha_in%lam(a_in)
            if( ((-1)**(l3b+lam3b)) /= ((-1)**(lin+lamin)) ) cycle 
            n2ch%max = n2ch%max + 1 
            abar(a_in,a_cdcc) = n2ch%max
         
      
         
         end do
         
          
      end do 
      
      if(allocated(n2ch%a3b)) deallocate(n2ch%a3b)
      if(allocated(n2ch%ain)) deallocate(n2ch%ain)
      
      write(8,*) "there are n2ch%max=", n2ch%max, "coupling channels"
      
      allocate(n2ch%a3b(n2ch%max),n2ch%ain(n2ch%max) )
      n2ch%a3b=0
      n2ch%ain=0
      i2ch=0
      do a_cdcc=1, incdcc%nchmax
         nj_3b = nint(2. * incdcc%J(a_cdcc))
         l3b = incdcc%l(a_cdcc)
         lam3b = incdcc%lam(a_cdcc)
         
         do a_in=1, alpha_in%nchmax 
            nj_in = nint(2. * alpha_in%J(a_in))
            if(nj_3b /= nj_in) cycle 
            lin =  alpha_in%l(a_in)
            lamin = alpha_in%lam(a_in)
            if( ((-1)**(l3b+lam3b)) /= ((-1)**(lin+lamin)) ) cycle 
            i2ch = i2ch + 1 
            
            n2ch%a3b(i2ch) = a_cdcc
            n2ch%ain(i2ch) = a_in
         end do
         
          
      end do 




      end subroutine 
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine alpha3b_in()
c     index of  | \alpha >_{in}
c     index of {(alpha_2b (lam jt)j_spect ; J}
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: l,lambda,lmax_cdcc,lbound
      real*8 :: jmin,jmax
      real*8 :: s2b,j2b,j_spect,J,s
      integer :: nJ,nJ2b,nJ_spect
      integer :: nch !index of {(l(jxjb)s2b)J2b (lam jt)J_spect ; J}
      real*8 :: jtmin, jtmax
      integer,dimension(10) :: l2b
      namelist /cdccwf/ l2b,lmax_cdcc,jtmin, jtmax,lbound
      l2b=-99
      lbound=0
      open (unit=551,file='iav.in')
      rewind(551)
      read(551,nml=cdccwf)

      s2b=0.0d0 ! by assuming all spin zero particle, need to improved
      s=0.0d0 ! spin of the target


      alpha_in%nchmax=0
        l=lbound
        do nJ2b=nint(2.*(l+s2b)),nint(2.*abs(l-s2b)),2
           j2b=nJ2b/2.0d0
         do lambda=0,lmax_cdcc
           do nJ_spect=nint(2*abs(lambda-s)),nint(2*(lambda+s)),2
            J_spect=nJ_spect/2.0d0
            jmin=max(jtmin,abs(J2b-J_spect))
            jmax=min(jtmax,J2b+J_spect)
      	    do nj=nint(2.*jmin),nint(2.*jmax),2
      	      alpha_in%nchmax=alpha_in%nchmax+1
      	    end do
           end do
         end do
        end do

      write(8,10)alpha_in%nchmax
10    format('there are',I3,1X,'cdcc channels')

       if (allocated(alpha_in%l)) deallocate(alpha_in%l)
       if (allocated(alpha_in%lam)) deallocate(alpha_in%lam)
       if (allocated(alpha_in%J)) deallocate(alpha_in%J)
       if (allocated(alpha_in%s2b)) deallocate(alpha_in%s2b)
       if (allocated(alpha_in%J2b)) deallocate(alpha_in%J2b)
       if (allocated(alpha_in%j_spect)) deallocate(alpha_in%j_spect)


       allocate(alpha_in%l(1:alpha_in%nchmax))
       allocate(alpha_in%lam(1:alpha_in%nchmax))
       allocate(alpha_in%J(1:alpha_in%nchmax))
       allocate(alpha_in%s2b(1:alpha_in%nchmax))
       allocate(alpha_in%J2b(1:alpha_in%nchmax))
       allocate(alpha_in%j_spect(1:alpha_in%nchmax))


      nch=1
      write(8,20)
      write(8,30)

        l=lbound
        do nJ2b=nint(2.*(l+s2b)),nint(2.*abs(l-s2b)),2
           j2b=nJ2b/2.0d0
         do lambda=0,lmax_cdcc
           do nJ_spect=nint(2*abs(lambda-s)),nint(2*(lambda+s)),2
            J_spect=nJ_spect/2.0d0
            jmin=max(jtmin,abs(J2b-J_spect))
            jmax=min(jtmax,J2b+J_spect)
      	    do nj=nint(2.*jmin),nint(2.*jmax),2
               j=nj/2.0d0
      	       alpha_in%l(nch)=l
               alpha_in%lam(nch)=lambda
               alpha_in%J(nch)=J
               alpha_in%s2b(nch)=s2b
               alpha_in%J2b(nch)=J2b
               alpha_in%j_spect(nch)=j_spect

               write(8,40)nch,0,l,s,s,s2b,J2b,
     +                       lambda,s,j_spect,j
      	           nch=nch+1
      	    end do
         end do
      end do
      end do


20    format('---For incoming channels the coupling coefficients are')
30    format(' a3b','|', 'a2b','|','(',' l ','(',' jb ',' jx ',')',
     +      ' s2b ',')',' J2b ', '(',' lam ',' jt ',')',' J3 ',',',
     +       ' Jtot ')
40    format(I4,I4,1x,I3,2x,f3.1,2x,f3.1,2x,f3.1,2x,
     +       f4.1,3x,I3,1x,f3.1,1x,f4.1,3x,f4.1)
      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine alpha_cdcc_in()
c     index of  | \alpha >_{in}
c     index of {(alpha_2b (lam jt)j_spect ; J}
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: l,lambda,lmax_cdcc,lbound
      real*8 :: jmin,jmax
      real*8 :: s2b,j2b,j_spect,J,s
      integer :: nJ,nJ2b,nJ_spect
      integer :: nch !index of {(l(jxjb)s2b)J2b (lam jt)J_spect ; J}
      real*8 :: jtmin, jtmax
      integer,dimension(10) :: l2b
      integer :: nl
      namelist /cdccwf/ l2b,lmax_cdcc,jtmin, jtmax, lbound
      l2b=-99
      open (unit=551,file='iav.in')
      rewind(551)
      read(551,nml=cdccwf)

      s2b=0.0d0 ! by assuming all spin zero particle, need to improved
      s=0.0d0 ! spin of the target


      incdcc%nchmax=0
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

         allocate(indexthox(1:exmax,1:n2ch%max))
         indexthoxmax=0
        do ie=1, exmax
          do acdcc=1, n2ch%max
             indexthoxmax=indexthoxmax+1
             indexthox(ie,acdcc)=indexthoxmax
          end do
        end do












      end subroutine





      end module
