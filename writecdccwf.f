      module writecdccwf
      use cdccchannels
      implicit none

C     complex*16, dimension(:,:,:),allocatable :: cdccwf_thox
      complex*16, dimension(:,:),allocatable :: cdccwf_thox_phi
      complex*16, dimension(:,:),allocatable :: cdccwf_thox_chi
C     complex*16, dimension(:,:,:),allocatable :: cdccwf_NEB


      contains

      subroutine initial_cdcc_wf_smoothie()
        use channels, only: jpiset,jpsets
        use parameters, only: maxeset
        use xcdcc,only:nquad,frad,rquad,nrcc,rvcc,parch
       implicit none

C      allocate(cdccwf_thox(1:nquad,1:nrcc,1:indexthoxmax))  ! note that rvcc(ir)=dble(ir-1)*hcm
C      allocate(cdccwf_NEB(1:nquad,1:nrcc,1:incdcc%nchmax))  ! note that rvcc(ir)=dble(ir-1)*hcm
       allocate(cdccwf_thox_phi(1:nquad,1:indexthoxmax))
       allocate(cdccwf_thox_chi(1:nrcc,1:indexthoxmax))
C      cdccwf_thox=0.0d0
C      cdccwf_NEB=0.0d0
       cdccwf_thox_phi=0.0d0
       cdccwf_thox_chi=0.0d0

      end subroutine








c *** -----------------------------------------------------------
c *** write CDCC wave function in each alpha n
c *** ------------------------------------------------------------
      subroutine cdcc_wf_thoxin(icc,incvec,ninc,nch,nr)
        use channels, only: jpiset,jpsets
        use parameters, only: maxeset
        use xcdcc,only:nquad,frad,rquad,nrcc,rvcc,parch,wfcdcc
      use channels, only: jptset,jpsets
      implicit none
      integer, parameter:: kwf=322
      integer :: ich,nch,nr,ir,icc,partot,ninc
      integer :: lsp,jsp,jc,parc,lam,iex,parp, lam_in
      real*8 :: jtot,exc,jp,jlp,jp_in 
      logical,dimension(1:nch) :: incvec
      complex*16 :: wf(nch,nr)

      integer :: jset,ne,nchan,iparity,nst,ii,ie
      real*8  :: xjtot
      integer :: alphacdcc,irbx,ira
      complex*16 :: phi, chi
      
      integer :: inc
      integer,dimension(ninc) :: inc_ich
      integer :: ain,nrank
      integer :: ijpset, njpsets
      integer,dimension(1:99) :: jpsetvalue
      
      
      do inc=1, ninc 
        do ich=1, nch
         if (incvec(ich)) then 
         inc_ich(inc) = ich 
         incvec(ich)=.false.
         exit
         else 
            cycle 
         end if 
         end do 
      end do 
      





      partot=jptset(icc)%partot
      jtot  =jptset(icc)%jtot



      write(kwf,350) icc,nch,jtot,partot
350   format(2x,"CC set:", i4, " Chans:",i4,
     &  " JTOT=",1f6.1, " Parity:",i2)

      do inc=1, ninc 
        jp_in= jptset(icc)%jp (inc_ich(inc)) 
        lam_in= jptset(icc)%l  (inc_ich(inc))
        
        do ain=1, alpha_in%nchmax
           if (nint(2.*jtot) == nint(2.*alpha_in%j(ain))
     +     .and.  nint(2.*jp_in) == nint(2.*alpha_in%j2b(ain))
     +     .and. lam_in == alpha_in%lam(ain)) exit
        end do 
        
        

      do ich=1,nch
       iex   =jptset(icc)%idx(ich)
       lam   =jptset(icc)%l  (ich)
       jp    =jptset(icc)%jp (ich)
       jlp   =jptset(icc)%jlp(ich)
       parp  = parch(iex)

       
       njpsets=0
       jpsetvalue=0
       do jset=1, jpsets
        ne=jpiset(jset)%nex
        nchan=jpiset(jset)%nchan
        xjtot=jpiset(jset)%jtot
        iparity=jpiset(jset)%partot
        if (nint(2.*xjtot)== nint(2.*jp) .and. iparity == parp) then 
           njpsets=njpsets+1
           jpsetvalue(njpsets) = jset
        end if 
       end do
       
       nrank=0
       do ijpset=1,njpsets
        nrank=jpiset(jpsetvalue(ijpset))%nex+nrank 
       end do 
C      if (ne /= nch) then 
C       write(*,*) "error in channel recouping"
C       write(*,*) "ne=",ne,"nch=",nch,"njpsets=",njpsets
C       stop
C      end if 
     
       ne=0
       do ijpset=1, njpsets    
       ne=jpiset(ijpset)%nex+ne 
C      if(ich>ne .and. njpsets>1) cycle !!!!!!?????????? do not understand why I need this ?????
       jset= jpsetvalue(ijpset)
       
       

       nchan=jpiset(jset)%nchan
       if(nchan/=1) stop "error in channel couplings"
       nst=0
       do ii=1, jset-1
         nst=nst+jpiset(ii)%nex
       end do



C      write(333,*)  "ich=",ich, "iex=",iex, "nst=",nst
C      if (iex>ne) cycle 
C      if (nst>=iex) cycle     
       if (iex>nst) ie=iex-nst
C      write(320,*) "inc=",inc,"ich=",ich,"ie=",ie,"iex=",iex,"nst=",nst
       
            


       do alphacdcc=1,incdcc%nchmax

        if (nint(2.*jtot) == nint(2.*incdcc%j(alphacdcc))
     +  .and.  nint(2.*jp) == nint(2.*incdcc%j2b(alphacdcc))
     +  .and. lam == incdcc%lam(alphacdcc)) exit
       end do
       
C      cycle 


       incdcc%n(alphacdcc)=nrank
C      write(319,*) "alphacdcc=",alphacdcc,"iex=",iex,"ie=",ie,
C    +              "nex=",jpiset(jset)%nex


       if (alphacdcc > incdcc%nchmax) stop "error in channels"




      write(kwf,354) ich,incdcc%n(alphacdcc),ie,jset,lam,
     +               jlp,jp,alphacdcc
354   format(2x," Channel:",i4, " n: ",i3, " IE:", i3, " jset:",I3
     &       " QN:",i3,1f5.1,1f5.1, " alphacdcc:", i3)





C      do irbx=1,nquad
C        do ira=1,nrcc
C
C        phi=frad(jset,ie,irbx,nchan)/rquad(irbx)
C        chi=wfcdcc(inc,ich,ira)
C
C        cdccwf_thox(irbx,ira,indexthox(ie,abar(ain,alphacdcc)))
C    +          =phi*chi
C        end do
C      end do
C
c-------new out put set
C      write(186,*) "indexthox=", indexthox(ie,abar(ain,alphacdcc)),
C    +      "n2ch=",abar(ain,alphacdcc)
       do irbx=1,nquad
         
        
         
         phi=frad(jset,ie,irbx,nchan)/rquad(irbx)
       if(njpsets>1) cdccwf_thox_phi(irbx,
     +               indexthox(iex,abar(ain,alphacdcc)))=phi
       if(njpsets==1) cdccwf_thox_phi(irbx,
     +               indexthox(ie,abar(ain,alphacdcc)))=phi
         

       end do

       do ira=1,nrcc
         chi=wfcdcc(inc,ich,ira)
         if(njpsets>1) cdccwf_thox_chi(ira,
     +                 indexthox(iex,abar(ain,alphacdcc)))=chi
         if(njpsets==1)cdccwf_thox_chi(ira,
     +                 indexthox(ie,abar(ain,alphacdcc)))=chi

       end do
       
       
       
       end do  ! ijpset



      enddo !ich
      
      end do ! inc 
      
      
      open (unit=222,file='rgrid.dat')
      do irbx=1, nquad
         write(222,*) rquad(irbx)
      end do 




      end subroutine
c *** -----------------------------------------------------------
      subroutine cdcc_wf_smoothie_out()
        use channels, only: jpiset,jpsets
        use parameters, only: maxeset
        use xcdcc,only:nquad,frad,rquad,nrcc,rvcc,parch
        implicit none
        integer :: ie, acdcc,ira,irbx
        
        
        
C       open (unit=551,file='cdccwf.dat')
C
C       do acdcc=1, incdcc%nchmax
C         write(551,*) acdcc
C        do ie=1, exmax
C
C
C
C         cdccwf_NEB(:,:,acdcc)= cdccwf_NEB(:,:,acdcc)
C    +                         + cdccwf_thox(:,:,indexthox(ie,acdcc))
C
C        end do
C        do ira=1,nrcc
C
C          do irbx=1,nquad
C            write(551,*) cdccwf_NEB(irbx,ira,acdcc)
C          end do
C        end do
C
C       end do
c-----------------new output
        open (unit=552,file='cdccwf_phi.dat')
        open (unit=553,file='cdccwf_chi.dat')
        open (unit=554,file='cdccwf_ranks.dat')
        do acdcc=1, n2ch%max
          write(552,*) acdcc
          write(553,*) acdcc
          write(554,*) incdcc%n(n2ch%a3b(acdcc))

         do ie=1, incdcc%n(n2ch%a3b(acdcc))
          do ira=1,nrcc
             write(553,*) cdccwf_thox_chi(ira,indexthox(ie,acdcc))
          end do
          do irbx=1,nquad
            write(552,*) cdccwf_thox_phi(irbx,indexthox(ie,acdcc))
          end do
         end do

        end do







      end subroutine


















      end module
