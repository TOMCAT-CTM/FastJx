      program xsection
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----addX_FJX.f  -- with user supplied subroutine that supplies X-section x Q-yield
C----     generate fast-JX 18-bin X-sections
c---------revised and updated for (mprather,2/2013)
c
c     Edited 1/7/2020 by Martyn Chipperfield for processing FastJx files 
c         for TOMCAT or UKCA/UKESM 
c
c     Edited 1/9/2023 to solve Fortran issue of first call to XSGAS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer, parameter :: NB_ = 100
      integer, parameter :: NS_ = 40000
      integer, parameter :: NZ_ = 13550
      integer, parameter :: NJ_ = 18

      real*8 SRB(15,NJ_)
      real*8, dimension(NB_+1) :: WBIN
      real*8, dimension(NB_) :: FBIN, ABIN
      real*8, dimension(NJ_) :: FFBIN,AABIN
      integer IJX(NB_), ITT
      integer NB,I,J,J1,J2,K,K1,K2
      integer INIT
      real*8 W(NS_),F(NS_)
      integer IBINJ(NS_)
      real*8 WZ(NZ_),X(NZ_,3)
      real*8 W1,W2, TT,XP,XM, WW,XNEW, XDUM
      character*6 TITLNEW
      character*4 TITLET

      integer NQY
      parameter (NQY=1)
      integer IQY


c     Information about wavelength bins
      open (1, file='wavel-bins.dat', status='OLD')
        SRB(:,:) = 0.d0
        read(1,'(i5)') NB
        if (NB .gt. NB_) stop
        read(1,'(5x,f8.3)') (WBIN(I), I=1,NB+1)
        read(1,*)
        read(1,*)
        read(1,'(2x,15f5.1)') ((SRB(I,J),I=1,15),J=1,8)
        read(1,*)
        read(1,'(5x,i5)') (IJX(I),I=16,NB)
      close (1)

      open (2, file='solar-p05nm-UCI.dat', status='OLD')
        read(2,*)
        read(2,*)
        read(2,'(f10.4,e10.3)') (W(J),F(J), J=1,NS_)
      close (2)

      open (3, file='XO3-p05nm-UCI.dat', status='OLD')
        read(3,*)
        read(3,*)
        read(3,'(f10.4,3e10.3)') (WZ(J),X(J,1),X(J,2),X(J,3), J=1,NZ_)
      close (3)

c     Initialization call to user subroutine (IQY not set)
c     For this first call the 5th argument is just a dummy scalar variable.
      INIT = 0
      call XSGAS(WW,TT,XP,XM,XDUM,INIT,TITLNEW,IQY)

c     Synchronize with the O3 cross sections (whether done or not)
c     Will loop K=1:NZ_  (J = K - 1 + K1), wavel = WZ(K)
      do J=1,NS_
        if (WZ(1) .eq. W(J)) goto 10
      enddo
      write(6,*) 'Cannot synch the solar & Xsections'
      stop
c
   10 K1 = J
      K2 = min( NS_, NZ_+K1-1)
c     write(6,'(a,2f10.3,i10)') ' synch:',WZ(1),W(K1),K1
c     write(6,'(a,2f10.3,i10)') ' synch:',WZ(NZ_),W(K2),K2

c     Now assign bin #(I=1:77) to each p05nm microbin J (1:40000)
      IBINJ(:) = 0
      do I=1,NB
        W1 = WBIN(I)
        W2 = WBIN(I+1)
        do J=1,NS_
          if (W(J) .gt. W1) goto 11
        enddo
        J = NS_ + 1
   11   J1 = J
        do J=J1,NS_
          if (W(J) .gt. W2) goto 12
        enddo
        J = NS_ + 1
   12   J2 = J-1
        do J=J1,J2
          IBINJ(J) = I
        enddo
c       write(6,'(i5,2f9.3,2i9,2f9.3)') I, W1,W2, J1,J2,W(J1),W(J2)
      enddo
c     Be aware that this binning does not interpolate and is OK for large bins
c     It has 7% error in the very short wavel S-R bins of pratmo.
c     It should be fine for weighting cross sections!

!
!     This looping is set for temperature only
      XP = 999.d0
      XM = 2.46d19

c     Loop over different channels
      do 100 IQY=1,NQY

c     Major temperature-density loop
      do ITT =1,2

        if    (ITT.eq.1) then
          TITLET = ' 200'
          TT = 200.d0
        elseif(ITT.eq.2) then
          TITLET = ' 300'
          TT = 300.d0
        else
          stop 'T loop'
        endif
c
c---    Now ready to do any flux-weighted means over the 77-pratmo bins
        FBIN(:) = 0.d0
        ABIN(:) = 0.d0

c       Primary high-resolution wavelength loop
        do J = K1,K2
          K = J - K1 + 1
          I = IBINJ(J)
          if (I .gt. 0) then

cc          write(6,*) 'before call',K1,K2,I,W(J),IQY
            call XSGAS(W(J),TT,XP,XM,XNEW,INIT,TITLNEW,IQY)

            FBIN(I) = FBIN(I) + F(J)
            ABIN(I) = ABIN(I) + F(J)*XNEW
          endif
        enddo
c       Mean cross section (for this wavelength and temperature)
        do I=1,NB
          if (FBIN(I) .gt. 0.d0) ABIN(I) = ABIN(I)/FBIN(I)
        enddo
c---    Write out UCI std 77-bin data
c       write(6,'(a6,a4/(1p,8e10.3))') TITLNEW, TITLET, ABIN

c
c       Secondary sum 77-bin pratmo ==> 18-bin fast-JX
c       Combine fast-JX bins:
c       - non-SR bands (16:NB) are assigned a single JX bin
c       - SR bands are split (by Opacity Distrib Fn) into a range of JX bins
        FFBIN(:) = 0.d0
        AABIN(:) = 0.d0
        do I=16,NB
          J = IJX(I)
          FFBIN(J) = FFBIN(J) + FBIN(I)
          AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)
        enddo
        do I=1,15
        do J=1,NJ_
          FFBIN(J) = FFBIN(J) + FBIN(I)*SRB(I,J)
          AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)*SRB(I,J)
        enddo
        enddo
        do J=1,NJ_
          if (FFBIN(J) .gt. 0.d0) AABIN(J) = AABIN(J)/FFBIN(J)
        enddo

c       Save UCI fast-JX data bins
c       For BrONO2 write as two channels (0.85 (-> NO3) and 0.15 (-> NO2))
        write(6,'(a6,a4,1p,6e10.3/10x,6e10.3/10x,6e10.3)') TITLNEW, TITLET, AABIN

      enddo
c
c     End of loop over quantum yield
 100  continue
c
      stop
c
      end

c-------------sample subroutine for fast-JX Xsection generation---------
c-----------------------------------------------------------------------
      subroutine XSGAS(WW, TT, PP, MM, XXWT, INIT, TITLNEW, IQY)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     WW - wavelength  (nm)
c     TT - temperature (K)
c     IQY - Quantum yield index
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8, intent(in) :: WW, TT, PP, MM
      integer, intent(in):: IQY
      integer, intent(inout) :: INIT
      real*8, intent(out) :: XXWT
      character*6, intent(inout) :: TITLNEW

      character*80 FTBL,TABLE
      real*8 W(999), XW(999),XWU(999)
      real*8 WWL,XXW,FW
      real*8 QY, QY1, TTL
      integer NW, NB, N, I, IW
c
      save NW, W, XW

c     Cross section file
      FTBL = 'jpl_2015_so3.dat'

      if(INIT.eq.0) then
        open (3, file=FTBL, status='OLD')
        read(3,'(a)') TABLE
        write(6,'(2a/a)') ' openfile=',FTBL, TABLE
        read(3,*) NW
        do N=1,NW
          read(3,*) W(N),XW(N)
        enddo
        close(3)
        INIT = 1
      else

c       No temperature dependence.
        do N=1,NW
c         Copy of cross sections
          XWU(N) = XW(N) 
        enddo

c       Quantum yield = 1
        QY1=1.0
c       Value to use
        if    (IQY.eq.1) then
          TITLNEW = 'jso3 '
          QY=QY1
        else
          stop 'QY'
        endif

c       Interpolate X-section versus wavelength
        if(WW.lt.W(1).or.WW.gt.W(NW)) then
c
c         Return zero cross section
          XXWT=0.0
        else
c
c         Interpolate in range
          IW = 1
          do I=2,NW-1
            if (WW .gt. W(I)) IW = I
          enddo
          FW = (WW - W(IW))/(W(IW+1) - W(IW))
          FW = min(1.d0, max(0.d0, FW))
          XXW = XWU(IW) + FW*(XWU(IW+1)-XWU(IW))
c
c         Full cross section including quantum yield (cm2)
          XXWT = QY * XXW * 1.e-20
c         write (6,*) qy,xxw
c         stop

        endif

      endif

      return
      end
