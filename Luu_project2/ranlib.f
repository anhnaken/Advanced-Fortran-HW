C
C  function ran3
C
C  generates uniformly distributed random numbers 
C  between 0.0 and 1.1
C  taken from Numerical Recipes (fortran version) 
C
C  common block need to store inext, etc 
C
C INPUT: IDUM  integer "seed" 
C   note: best if IDUM is negative first call, in order 
C         to initialize
C OUTPUT: ran3 
C
      real FUNCTION RAN3(IDUM)
C         IMPLICIT REAL*4(M)
C         PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=2.5E-7)
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      common/randumb/inext,inextp,ma     ! needed on HP
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END


c**************************************************************
c      returns Gaussian random variable.                      *
C      Algorithm is that on p. 209-210 Computational Physics  *
C      by Koonin and Meridith (fortran version).              *
c      Note that the algorithm generates two variables,       *
C      here we randomly select between the two.               *
c**************************************************************
C  INPUT: iseed
C      integer seed; for best use, should initially be < 0 
C  OUTPUT: gaussvar, a gaussian distributed random number
C           with unit variance
C  FUNCTIONS CALLED: ran3
C
 
      real function gaussvar(iseed)
      implicit real (a-h,o-z)
      pi=2.00*asin(1.0)
      two=-2.0*log(1.0-ran3(iseed))
      radius=sqrt(two)
      theta=2.0*pi*ran3(iseed)
      gauss1=radius*cos(theta)
      gauss2=radius*sin(theta)
      xchoose=ran3(iseed)
      if(xchoose.le.0.50)then
         gaussvar=gauss1
      else
         gaussvar=gauss2
      end if
      return
      end
