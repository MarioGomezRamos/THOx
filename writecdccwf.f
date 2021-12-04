      module writecdccwf
      use cdccchannels
      implicit none

      complex*16, dimension(:,:,:),allocatable :: cdccwf_thox
      complex*16, dimension(:,:),allocatable :: cdccwf_thox_phi
      complex*16, dimension(:,:),allocatable :: cdccwf_thox_chi
      complex*16, dimension(:,:,:),allocatable :: cdccwf_NEB


      contains

      subroutine initial_cdcc_wf_smoothie()
        use channels, only: jpiset,jpsets
        use parameters, only: maxeset
        use xcdcc,only:nquad,frad,rquad,nrcc,rvcc,parch
       implicit none

       allocate(cdccwf_thox(1:nquad,1:nrcc,1:indexthoxmax))  ! note that rvcc(ir)=dble(ir-1)*hcm
       allocate(cdccwf_NEB(1:nquad,1:nrcc,1:incdcc%nchmax))  ! note that rvcc(ir)=dble(ir-1)*hcm
       allocate(cdccwf_thox_phi(1:nquad,1:indexthoxmax))
       allocate(cdccwf_thox_chi(1:nrcc,1:indexthoxmax))
       cdccwf_thox=0.0d0
       cdccwf_NEB=0.0d0
       cdccwf_thox_phi=0.0d0
       cdccwf_thox_chi=0.0d0

      end subroutine








c *** -----------------------------------------------------------
c *** write CDCC wave function in each alpha n
c *** ------------------------------------------------------------
      subroutine cdcc_wf_thoxin(icc,wf,nch,nr)
        use channels, only: jpiset,jpsets
        use parameters, only: maxeset
        use xcdcc,only:nquad,frad,rquad,nrcc,rvcc,parch
      use channels, only: jptset,jpsets
      implicit none
      integer, parameter:: kwf=322
      integer :: ich,nch,nr,ir,icc,partot
      integer :: lsp,jsp,jc,parc,lam,iex,parp
      real*8 :: jtot,exc,jp,jlp
      complex*16 wf(nch,nr)

      integer :: jset,ne,nchan,iparity,nst,ii,ie
      real*8  :: xjtot
      integer :: alphacdcc,irbx,ira
      complex*16 :: phi, chi






      partot=jptset(icc)%partot
      jtot  =jptset(icc)%jtot



      write(kwf,350) icc,nch,jtot,partot
350   format(2x,"CC set:", i4, " Chans:",i4,
     &  " JTOT=",1f6.1, " Parity:",i2)

      do ich=1,nch
       iex   =jptset(icc)%idx(ich)
       lam   =jptset(icc)%l  (ich)
       jp    =jptset(icc)%jp (ich)
       jlp   =jptset(icc)%jlp(ich)
       parp  = parch(iex)


       do jset=1, jpsets
        ne=jpiset(jset)%nex
        nchan=jpiset(jset)%nchan
        xjtot=jpiset(jset)%jtot
        iparity=jpiset(jset)%partot
        if (nint(2.*xjtot)== nint(2.*jp) .and. iparity == parp) exit
       end do


       nchan=jpiset(jset)%nchan
       if(nchan/=1) stop "error in channel couplings"
       nst=0
       do ii=1, jset-1
         ne=jpiset(ii)%nex
         nst=nst+ne
       end do



       if (iex>nst) ie=iex-nst


       do alphacdcc=1,incdcc%nchmax

        if (nint(2.*jtot) == nint(2.*incdcc%j(alphacdcc))
     +  .and.  nint(2.*jp) == nint(2.*incdcc%j2b(alphacdcc))
     +  .and. lam == incdcc%lam(alphacdcc)) exit
       end do

       incdcc%n(alphacdcc)= jpiset(jset)%nex

C       write(322,*) "n=",incdcc%n(alphacdcc)



       if (alphacdcc > incdcc%nchmax) stop "error in channels"




      write(kwf,354) ich,incdcc%n(alphacdcc),ie,jset,lam,
     +               jlp,jp,alphacdcc
354   format(2x," Channel:",i4, " n: ",i3, " IE:", i3, " jset:",I3
     &       " QN:",i3,1f5.1,1f5.1, " alphacdcc:", i3)





       do irbx=1,nquad
         do ira=1,nrcc

         phi=frad(jset,ie,irbx,nchan)/rquad(irbx)
         chi=wf(ich,ira)

         cdccwf_thox(irbx,ira,indexthox(ie,alphacdcc))=phi*chi
         end do
       end do

c-------new out put set
       do irbx=1,nquad

         phi=frad(jset,ie,irbx,nchan)/rquad(irbx)
         cdccwf_thox_phi(irbx,indexthox(ie,alphacdcc))=phi

       end do

       do ira=1,nrcc
         chi=wf(ich,ira)
         cdccwf_thox_chi(ira,indexthox(ie,alphacdcc))=chi

       end do



      enddo !ich
      
      
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
        do acdcc=1, incdcc%nchmax
          write(552,*) acdcc
          write(553,*) acdcc
          write(554,*) incdcc%n(acdcc)

         do ie=1, incdcc%n(acdcc)
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
