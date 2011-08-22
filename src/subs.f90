! --Yipu Miao 05/09/2010
!********************************************************
! This file list common subroutines
! Subroutine List:
! upcase(string,iend)         : upcase string
! rdword(STRING,ISTART,ISTOP) : read word form spcified position
! rdnum                       : convert string to number
! rdinum                      : convert string to integer
! xnorm(a,i,j,k)              : Normalization
! xnewnorm
! priTim                      : Print Time
! priSym                      : Print Symmetry Matrix
! DIAG                        : Diagnalization
! EIGVEC
! DEGEN
! EIGVAL
! TRIDI                       : The above four subroutines are for eigen problem
! RANDOM                      : Random
! PriLab                      : elegent way to output atom list
! EffChar                     : check valid char
!
!********************************************************
! $Id: upcase.f90,v 1.1.1.1 2007/01/26 20:22:35 ayers Exp $
!--------------------------------------------------------
! Upcase(string,iend)
!--------------------------------------------------------
!
    SUBROUTINE UPCASE(STRING,IEND)

! CHANGES LOWER CASE CHARACTERS IN STRING(1:IEND) TO UPPER CASE.
! -- S. DIXON.

    CHARACTER STRING*(*),AUPP,ALOW,DUMMY
    INTEGER :: IALOW,IAUPP,IDUMMY
    DATA AUPP,ALOW /'A','a'/
    SAVE AUPP,ALOW
    IALOW = ICHAR(ALOW)
    IAUPP = ICHAR(AUPP)
    ISHIFT = IALOW - IAUPP
    DO I=1,IEND
        DUMMY = STRING(I:I)
        IDUMMY = ICHAR(DUMMY)
        IF(IDUMMY >= IALOW)THEN
            DUMMY = CHAR(IDUMMY - ISHIFT)
            STRING(I:I) = DUMMY
        ENDIF
    ENDDO
    RETURN
    end SUBROUTINE UPCASE
!********************************************************
!
!********************************************************
! rdword subroutine
!------------------------------------------------------------------------
! This was taken from the DivCon code by Steve Dixon.
! - Ed Brothers, February 18, 2002
!------------------------------------------------------------------------
!
    SUBROUTINE RDWORD(STRING,ISTART,ISTOP)

! LOCATES THE NEXT WORD IN STRING STARTING AT STRING(ISTART:ISTART).
! ON RETURN ISTART AND ISTOP WILL BE UPDATED AND THE WORD WILL BE
! IN STRING(ISTART:ISTOP).  IF THERE ARE NO MORE WORDS IN THE STRING,
! THEN BOTH ISTART AND ISTOP WILL BE RETURNED WITH VALUES OF ZERO.

    CHARACTER STRING*(*)
    LOGICAL :: INWORD

! GET DECLARED LENGTH OF STRING AND FIND THE NEXT CONTIGUOUS BLOCK
! OF NONBLANK CHARACTERS.

    IBEGIN = MAX(ISTART,1)
    IEND = LEN(STRING)
    IF(IEND < IBEGIN)THEN
        ISTART = 0
        ISTOP = 0
        RETURN
    ENDIF
    INWORD = .FALSE.
    IBEGIN = ISTART
    DO 100 I=IBEGIN,IEND
        IF(STRING(I:I) == ' ')THEN
            IF(INWORD)THEN
                ISTOP = I-1
                RETURN
            ENDIF
        ELSE
            IF( .NOT. INWORD)THEN
                INWORD = .TRUE.
                ISTART = I
            ENDIF
        ENDIF
    100 ENDDO

! IF WE GET HERE, THEN EITHER THE WORD FOUND EXTENDS ALL THE WAY TO
! THE END OF THE STRING, OR NO WORD WAS FOUND IN THE REMAINING
! PORTION OF THE STRING.

    IF(INWORD)THEN
        ISTOP = IEND
    ELSE
        ISTART = 0
        ISTOP = 0
    ENDIF
    RETURN
    end SUBROUTINE RDWORD

!********************************************************
! rdnum
! $Id: rdnum.f90,v 1.1.1.1 2007/01/26 20:22:35 ayers Exp $
!------------------------------------------------------------------------
! contains all the routines to convert strings to numbers
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    subroutine rdnum(string,istart,value,ierror)

! extracts a double precision floating point number from a character
! string.  the field of search starts at string(istart:istart) or
! after the first equals sign following the istart position.  the
! number is returned in value.  if an error is encountered, ierror
! is set to one.  this routine expects that there are no blank spaces
! embedded anywhere within the numerical field.

    implicit double precision (a-h,o-z)
    character string*(*),char*1,efield(4)*1
    data efield /'E','e','D','d'/
    save efield
    ierror = 0
    ibeg = istart
    istop = len(string)
    iend = istop
    do 10 i=istart,istop
        if(string(i:i) == ' ')then
            iend = i-1
            go to 20
        endif
    10 ENDDO
    20 if(iend < ibeg)then
        ierror = 1
        go to 1000
    endif
    ieq = index(string(ibeg:iend),'=')
    if(ieq /= 0) ibeg = ibeg + ieq
    call getnum(string,ibeg,iend,value,ierror)
    1000 return
    end subroutine rdnum



!c------------------------------------------------------------cc
!c------------------------------------------------------------cc

    subroutine getnum(string,ibeg,iend,value,ierror)
    implicit double precision(a-h,o-z)

    character string*(*),char*1,efield(4)*1
    data efield /'E','e','D','d'/
    save efield

!C-----------------------------------------------------CC

    value = 0.0d0

! check for algebraic sign.

    char = string(ibeg:ibeg)
    if(char == '-')then
        asign = -1.0d0
        ibeg = ibeg + 1
    elseif(char == '+')then
        asign = 1.0d0
        ibeg = ibeg + 1
    else
        asign = 1.0d0
    endif
    if(ibeg > iend)then
        ierror = 1
        go to 1000
    endif

! first determine the whole number equivalent of whatever is
! to the left of any decimal point.

    idecml = index(string(ibeg:iend),'.')
    if(idecml == 1)then
        if(ibeg == iend)then
        
        ! number is just a decimal point.  assume a value of zero.
        
            value = 0.0d0
            go to 1000
        endif
        xleft = 0.0d0
        ibeg = ibeg+1
    else
        i1 = ibeg
        if(idecml == 0)then
            i2 = iend
        else
            i2 = ibeg+idecml-2
        endif
        call whole(string,i1,i2,xleft,ierror)
        if(ierror /= 0) go to 1000
        value = xleft*asign
        if(idecml == 0 .OR. i2 == (iend-1)) go to 1000
        ibeg = i2+2
    endif

! determine the whole number equivalent of whatever is to the
! right of the decimal point.  account for e or d field format.

    do 30 i=1,4
        ie = index(string(ibeg:iend),efield(i))
        if(ie /= 0) go to 40
    30 ENDDO
    40 if(ie == 1)then
        value = xleft*asign
        ibeg = ibeg + 1
    else
        i1 = ibeg
        if(ie == 0)then
            i2 = iend
        else
            i2 = ibeg+ie-2
        endif
        call whole(string,i1,i2,xright,ierror)
        if(ierror /= 0) go to 1000
        xright = xright*10.0d0**(i1-i2-1)
        value = value + xright*asign
        if(ie == 0 .OR. i2 == (iend-1)) go to 1000
        ibeg = i2+2
    endif

! get the exponential portion.

    char = string(ibeg:ibeg)
    if(char == '-')then
        esign = -1.0d0
        ibeg = ibeg + 1
    elseif(char == '+')then
        esign = 1.0d0
        ibeg = ibeg + 1
    else
        esign = 1.0d0
    endif
    if(ibeg > iend) go to 1000
    i1 = ibeg
    i2 = iend
    call whole(string,i1,i2,expart,ierror)
    if(ierror /= 0) go to 1000
    value = value*10.0d0**(esign*expart)
    1000 return
    end subroutine getnum




!C------------------------------------------------------------CC
!C------------------------------------------------------------CC

    subroutine whole(string,ibeg,iend,value,ierror)

! returns the whole number in the field string(ibeg:iend).  only
! the numbers 0-9 are allowed to be present.

    implicit double precision (a-h,o-z)
    character string*(*)
    ierror = 0
    value = 0.0d0
    ichar0 = ichar('0')
    do 10 i=ibeg,iend
        idigit = ichar(string(i:i)) - ichar0
        if(idigit < 0 .OR. idigit > 9)then
            ierror = 1
            go to 1000
        endif
        value = 10.0d0*value + idigit
    10 ENDDO
    1000 return
    end subroutine whole


!********************************************************
! rdinum(string,istart,ivalue,ierror)
!-----------------------------------------------------------
! Read integer number
!------------------------------------------------------------
!

    subroutine rdinum(string,istart,ivalue,ierror)

    implicit double precision (a-h,o-z)
    character string*(*),char*1,efield(4)*1

    ierror = 0
    ibeg = istart
    istop = len(string)
    iend = istop
    do 10 i=istart,istop
        if(string(i:i) == ' ')then
            iend = i-1
            go to 20
        endif
    10 ENDDO
    20 if(iend < ibeg)then
        ierror = 1
        go to 1000
    endif
    ieq = index(string(ibeg:iend),'=')
    if(ieq /= 0) ibeg = ibeg + ieq
    call getinum(string,ibeg,iend,ivalue,ierror)
    1000 return
    end subroutine rdinum

!C------------------------------------------------------------CC
!C------------------------------------------------------------CC

    subroutine getinum(string,ib,ie,ivalue,ierror)

    implicit double precision (a-h,o-z)
    character string*(*),char*1,efield(4)*1

    ierror = 0
    ivalue = 0
    if (string(ib:ib) == '-') then
        ib = ib + 1
        isign = -1
    elseif (string(ib:ib) == '+') then
        ib = ib + 1
        isign = 1
    else
        isign = 1
    endif

    call iwhole(string,ib,ie,ivalue,ierror)
    ivalue = ivalue * isign

    return
    end subroutine getinum


!C------------------------------------------------------------CC
!C------------------------------------------------------------CC

    subroutine iwhole(string,ibeg,iend,ivalue,ierror)

! returns the whole number in the field string(ibeg:iend).  only
! the numbers 0-9 are allowed to be present.

    character string*(*)
    ivalue = 0
    ichar0 = ichar('0')
    do 10 i=ibeg,iend
        idigit = ichar(string(i:i)) - ichar0
        if(idigit < 0 .OR. idigit > 9)then
            ierror = 1
            go to 1000
        endif
        ivalue = 10*ivalue + idigit
    10 ENDDO
    1000 return
    end subroutine iwhole

!C------------------------------------------------------------CC
!C------------------------------------------------------------CC

    subroutine iatoi(str, istart, lstr, integ, ierror)

    implicit double precision(a-h,o-z)

    character str*(*), ch
    logical :: int, min

    integ = 0
    izero = ichar('0')
    nstr = len(str)
    do i=istart,nstr
        ch = str(i:i)
        call whatis2(ch, int, min)
        if ( .NOT. int) goto 20
    enddo
    20 lstr = i-1
    if (lstr == 0) return
    call getinum(str,istart,lstr,integ,ierror)

    end subroutine iatoi

!C------------------------------------------------------------CC
!C------------------------------------------------------------CC

    subroutine iatoimp(str, istart, lstr, integ, ierror)

    implicit double precision(a-h,o-z)

    character str*(*), ch
    logical :: int, min

    integ = 0
    izero = ichar('0')
    nstr = len(str)
    do i=istart,nstr
        ch = str(i:i)
        call whatis1i(ch, int)
        if ( .NOT. int) goto 20
    enddo
    20 lstr = i-1
    if (lstr == 0) return
    call getinum(str,istart,lstr,integ,ierror)

    end subroutine iatoimp
	
! $Id: whatis.f90,v 1.1.1.1 2007/01/26 20:22:35 ayers Exp $

!------------------------------------------------------------------------
    subroutine whatis1(this, float)

    implicit double precision(a-h,o-z)

    character this
    logical :: float

    float = .false.

    ithis = ichar(this)
    i0 = ichar('0')
    i9 = ichar('9')
    if ((ithis >= i0) .AND. (ithis <= i9)) then
        float = .true.
    elseif (this == '.') then
        float = .true.
    elseif (this == '-') then
        float = .true.
    elseif (this == '+') then
        float = .true.
    elseif (this == 'E') then
        float = .true.
    elseif (this == 'e') then
        float = .true.
    elseif (this == 'D') then
        float = .true.
    elseif (this == 'd') then
        float = .true.
    endif

    end subroutine whatis1

!C------------------------------------------------------------CC
!C------------------------------------------------------------CC

    subroutine whatis1i(this, int)

    implicit double precision(a-h,o-z)

    character this
    logical :: int

    int = .false.

    if (this == '-') then
        int = .true.
    elseif (this == '+') then
        int = .true.
    else
        ithis = ichar(this)
        i0 = ichar('0')
        i9 = ichar('9')
        if ((ithis >= i0) .AND. (ithis <= i9)) int = .TRUE. 
    endif

    end subroutine whatis1i

!C------------------------------------------------------------CC
!C------------------------------------------------------------CC

    subroutine whatis2(this, int, min)

    implicit double precision(a-h,o-z)

    character this
    logical :: int, min

    int = .false.
    min = .false.

    if (this == '-') then
        min = .true.
    else
        ithis = ichar(this)
        i0 = ichar('0')
        i9 = ichar('9')
        if ((ithis >= i0) .AND. (ithis <= i9)) int = .TRUE. 
    endif

    end subroutine whatis2

!C------------------------------------------------------------CC
!C------------------------------------------------------------CC

    subroutine whatis7(this,char,num,parl,parr,comma,eq,white)
    implicit double precision(a-h,o-z)

    character this
    logical :: char,num,parl,parr,comma,eq,white

!C--------------------------------------------CC

    char = .false.
    num = .false.
    parl = .false.
    parr = .false.
    comma = .false.
    eq = .false.
    white = .false.
    if (this == ' ') then
        white = .true.
    elseif (this == ',') then
        comma = .true.
    elseif (this == '=') then
        eq = .true.
    elseif (this == '(') then
        parl = .true.
    elseif (this == ')') then
        parr = .true.
    elseif (this == '/') then
        eq = .true.
    elseif (this == '+') then
        num = .true.
    elseif (this == '.') then
        num = .true.
    elseif (this == '-') then
        num = .true.
    else
        ithis = ichar(this)
        ia = ichar('a')
        iz = ichar('z')
        iaa = ichar('A')
        izz = ichar('Z')
        i0 = ichar('0')
        i9 = ichar('9')
        if (((ithis >= ia) .AND. (ithis <= iz)) .OR. &
        ((ithis >= iaa) .AND. (ithis <= izz))) then
            char = .true.
        elseif ((ithis >= i0) .AND. (ithis <= i9)) then
            num = .true.
        else
            white = .true.
        endif
    endif

    end subroutine whatis7


!********************************************************
! xnorm(a,i,j,k)
!-----------------------------------------------------------
! Normalization
!------------------------------------------------------------
! Ed Brothers. October 3, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    double precision function xnorm(a,i,j,k)
    use constants
    implicit double precision(a-h,o-z)

! The purpose of this function is to calculate the normalization
! constant of a gtf with exponent a and x,y,z exponents of i,j,k

! According to Mathematica, the simplied normalization constant is:

! 3/4 + (i + j + k)/2  3/4 + (i + j + k)/2
! 2                    a
! ut[3]= --------------------------------------------------------
! 1                  1                  1
! Sqrt[Gamma[- + i]] Sqrt[Gamma[- + j]] Sqrt[Gamma[- + k]]
! 2                  2                  2

! Constants:

!    pi=3.1415926535897932385
!    pito3half=5.568327996831707845284817982118835702014

    xnorm =   (2.d0*a)**(dble(i+j+k)/2.d0+.75d0)

! Before the Gamma function code, a quick note. All gamma functions are
! of the form:
! 1
! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
! 2

! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gaussians
! just requires a loop and multiplying by Pi^3/2

    gamma=1.d0
    DO L=1,i
        gamma = gamma * (dble(i - L) + .5d0)
    ENDDO
    DO L=1,j
        gamma = gamma * (dble(j - L) + .5d0)
    ENDDO
    DO L=1,k
        gamma = gamma * (dble(k - L) + .5d0)
    ENDDO

    xnorm = xnorm*pi**(-0.75d0)*gamma**(-.5d0)

    return
    END
!********************************************************
! xnewnorm(a,i,j,k)
!--------------------------------------------------------
! Ed Brothers. October 3, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    double precision function xnewnorm(i,j,k,nprim1,ddcoe,aaexp)
    use allmod
    implicit double precision(a-h,o-z)
    real*8 ddcoe(nprim1),aaexp(nprim)

    pi32 = pi**(1.5)
!    xnorm =   (2.d0*a)**(dble(i+j+k)/2.d0+.75d0)

        xnewnorm= 0.d0
        nxponent=-i-j-k
!        nxponent=-itype(1,Ibas)-itype(2,Ibas)-itype(3,Ibas)
        xponent=-1.5d0+dble(nxponent)
        DO Icon1 = 1,nprim1
            DO Icon2 = 1,nprim1
               xnewnorm = xnewnorm + ddcoe(Icon1)*ddcoe(Icon2) &
                *(aaexp(Icon1)+aaexp(Icon2))**xponent
!               print*,xnewnorm
            ENDDO
        ENDDO
!        print*,'*',xnewnorm
        gamma=1.d0
        DO L=1,i
            gamma = gamma * (dble(i - L) + .5d0)
        ENDDO
        DO L=1,j
            gamma = gamma * (dble(j - L) + .5d0)
        ENDDO
        DO L=1,k
            gamma = gamma * (dble(k - L) + .5d0)
        ENDDO
        xnewnorm = (xnewnorm*gamma*pi32)**(-.5d0)
!        print*,'**',gamma,pi32,xnewnorm
!        DO Icon1 = 1,ncontract(Ibas)
!            dcoeff(Icon1,Ibas) = dconew*dcoeff(Icon1,Ibas)
!        ENDDO

    return
    END


!********************************************************
! vett
!------------------------------------------------------------
! Alessandro GENONI 03/12/007
! Subroutine to build up the array kvett, whose elemts are kvett(i)=i*(i-1)/2
!
      Subroutine vett
       use ecpmod
       implicit double precision (a-h,o-z)
       do i=1,nbf12
         kvett(i)=i*(i-1)/2
       end do
       return
      end


!******  Print out total Time used by the job ******
!*     in parent program def. "real*8 Tim0,CPUTim"; "Tim0=CPUTim(0)" for initial time
      subroutine PrtTim(IOut,sec)
      Implicit Real*8(A-H,O-Z)

1000  Format(' Job cpu time:',I3,' days ',I2,' hours ',I2,' minutes ',F4.1,' seconds.')

      Time = sec
      NDays = (Time / (3600.0d0*24.0d0))
      Time = Time - (NDays*(3600.0d0*24.0d0))
      NHours = (Time / 3600.0d0)
      Time = Time - (NHours*3600.0d0)
      NMin = (Time / 60.0d0)
      Time = Time - (NMin*60.0d0)
      Write(IOut,1000) NDays, NHours, NMin, Time
      Return
      End
	  
!******* print action ****************
      subroutine PrtAct(io,line)
      implicit none
      integer io,L,leng,i
      character line*(*)

      leng=len(line)
	  L=0
	  write(io,*)
      do i=leng,1,-1
         if (line(i:i).ne.' ') then
            L=i
            exit
         endif
      enddo
      if (L.gt.0) then
         write(io,'(" @ ",a)') line(1:L)
      else
         write(io,*)
      endif
	  write(io,*)
	  call flush(io)
      return
      end
	  
!****** 2004.12.22 print column(c1:c2) of mat(n,n) to file(io) ******
      subroutine PriCol(io,n,mat,c1,c2,fm) ! format: f(x.y) x>7 sugg 12.5,12.7,14.9
      implicit none
      integer i,j,jj,n,io,c1,c2,n5,nf,nc,x,y,k
      real*8 mat(n,n)
      character fm*(*),ch,fm2*10
      character*40 fmt1,fmt2,fmt3,fmt4

      nc=c2-c1+1; n5=nc/5; nf=mod(nc,5)
      fm2=fm; ch=fm2(1:1); k=index(fm2,'.')
      read(fm2(2:k-1),*) x; read(fm2(k+1:10),*) y

      write(fmt1,101) ch,x,y; write(fmt2,102) nf,ch,x,y
 101  format('(i7,5',a1,i2,'.',i2,')')
 102  format('(i7,',i2,a1,i2,'.',i2,')')
      write(fmt3,103) x-7; write(fmt4,104) nf,x-7
 103  format('(3x,5(',i2,'x,i7))')
 104  format('(3x,',i2,'(',i2,'x,i7))')

      do jj=1,n5
         write(io,fmt3) (j,j=c1+(jj-1)*5,c1+jj*5-1)
         write(io,fmt1) (i,(mat(i,j),j=c1+(jj-1)*5,c1+jj*5-1),i=1,n)
!         if (jj.ne.n5.or.nf.ne.0) write(io,*)
      enddo

      if (nf.ne.0) then
         write(io,fmt4)(j,j=c1+n5*5,c2)
         write(io,fmt2) (i,(mat(i,j),j=c1+n5*5,c2),i=1,n)
      endif
      call flush(io)

      end
	  
!****** print symmetric mat(n,n) to file(io) ******
      subroutine PriSym(io,n,mat,fm) ! format: f(x.y) x>7 sugg 12.5,12.7,14.9
      implicit none
      integer j,jj,n,io,n5,nf,x,y,ini,ifi,k
      real*8 mat(n,n)
      character fm*(*),ch,fm2*10
      character*40 fmt1,fmt2,fmt3,fmt4

      n5=n/5
      nf=mod(n,5)
      fm2=fm
      ch=fm2(1:1)
      k=index(fm2,'.')
      read(fm2(2:k-1),*) x
      read(fm2(k+1:10),*) y

      write(fmt1,101) ch,x,y
      write(fmt2,102) nf,ch,x,y
 101  format('(i7,5',a1,i2,'.',i2,')')
 102  format('(i7,',i2,a1,i2,'.',i2,')')
      write(fmt3,103) x-7
      write(fmt4,104) nf,x-7
 103  format('(3x,5(',i2,'x,i7))')
 104  format('(3x,',i2,'(',i2,'x,i7))')

      do jj=1,n5
         ini=1+(jj-1)*5
         write(io,fmt3) (j,j=ini,jj*5)
         do k=1+(jj-1)*5,n
            ifi=min(jj*5,k)
            write(io,fmt1) k,(mat(k,j),j=ini,ifi)
         enddo
!         if (jj.ne.n5.or.nf.ne.0) write(io,*)
      enddo

      if (nf.ne.0) then
         ini=n-nf+1
         write(io,fmt4)(j,j=ini,n)
         do k=ini,n
            write(io,fmt2) k,(mat(k,j),j=ini,k)
         enddo
      endif
      call flush(io)

      end
!********************************************************
! Diag
!------------------------------------------------------------
! S. Dixon's diagonalization code from DivCon.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    SUBROUTINE DIAG(NDIM,A,NEVEC1,TOLERA,V,EVAL1,IDEGEN1,EVEC1, &
    IERROR)

! DRIVER ROUTINE FOR DIAGONALIZATION OF THE REAL, SYMMETRIC,
! MATRIX A.

! VARIABLES REQUIRED:
! ------------------

! NDIM = ORDER OF THE MATRIX A (I.E., A IS NDIM BY NDIM);

! A = REAL SYMMETRIC MATRIX TO BE DIAGONALIZED.  ONLY THE LOWER
! HALF OF A NEED BE FILLED WHEN CALLING.  A IS DESTROYED BY
! THIS ROUTINE.

! NEVEC1 = THE NUMBER OF EIGENVECTORS REQUIRED;

! TOLERA = TOLERANCE FACTOR USED IN THE QR ITERATION TO DETERMINE
! WHEN OFF-DIAGONAL ENTRIES ARE ESSENTIALLY ZERO.  (DEFAULT
! IS 1.0D-8).

! V = 3 BY NDIM WORKSPACE.


! VARIABLES RETURNED:
! ------------------

! EVAL1 = EIGENVALUES OF A (SORTED IN INCREASING ALGEBRAIC VALUE);

! IDEGEN1 = DEGENERACIES (NUMBER OF TIMES REPEATED) FOR EIGENVALUES;

! EVEC1 = EIGENVECTORS OF A (IN COLUMNS OF EVEC1);

! ERROR CODES:  IERROR=0 - SUCCESSFUL CALCULATION.
! IERROR=1 - NO CONVERGENCE IN QR ITERATION.


! PROGRAMMED BY S. L. DIXON, OCT., 1991.


!    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION A(NDIM,*),V(3,*),EVAL1(*),IDEGEN1(*),EVEC1(NDIM,*)
!    DIMENSION A(nbasis,nbasis),V(3,nbasis),EVAL1(nbasis), &
!    IDEGEN1(nbasis),EVEC1(nbasis,nbasis)

    DIMENSION A(NDIM,NDIM),V(3,NDIM),IDEGEN1(NDIM),EVAL1(NDIM)
	DIMENSION EVEC1(NDIM,NDIM)

! FLAG FOR WHETHER OR NOT TO COMPUTE EIGENVALUES:
		

    IF(NDIM == 1)THEN
        EVAL1(1) = A(1,1)
        IDEGEN1(1) = 1
        EVEC1(1,1) = 1.0D0
        RETURN
    ENDIF

! TRIDIAGONALIZE THE MATRIX A.  THIS WILL OVERWRITE THE DIAGONAL
! AND SUBDIAGONAL OF A WITH THE TRIDIAGONALIZED VERSION.  THE
! HOUSEHOLDER VECTORS ARE RETURNED IN THE ROWS ABOVE THE DIAGONAL,
! AND THE BETAHS ARE RETURNED BELOW THE SUBDIAGONAL.

    CALL TRIDI(NDIM,V,A)

! COMPUTE NORM OF TRIDIAGONAL MATRIX FROM THE "LARGEST" COLUMN.

    ANORM = ABS(A(1,1)) + ABS(A(2,1))
    IF(NDIM > 2)THEN
        DO 20 I=2,NDIM-1
            AICOL = ABS(A(I-1,I)) + ABS(A(I,I)) + ABS(A(I+1,I))
            ANORM = MAX(ANORM,AICOL)
        20 ENDDO
    ENDIF
    ANCOL = ABS(A(NDIM-1,NDIM)) + ABS(A(NDIM,NDIM))
    ANORM = MAX(ANORM,ANCOL)

! GET EIGENVALUES AND DEGENERACIES OF THE TRIDIAGONAL MATRIX A.
! IF THE CALLING ROUTINE HAS NOT SUPPLIED A TOLERANCE FACTOR FOR
! OFF-DIAGONAL ENTRIES IN THE QR ITERATION, A DEFAULT OF 1.0D-8
! WILL BE USED.

    TOLTMP = TOLERA
    IF(TOLTMP <= 0.0D0) TOLTMP = 1.0D-8
    CALL EIGVAL(NDIM,A,V,TOLTMP,ANORM,EVAL1,IERROR)
    IF(IERROR /= 0) RETURN

! DETERMINE DEGENERACIES OF EIGENVALUES.

    CALL DEGEN(NDIM,EVAL1,TOLTMP,ANORM,IDEGEN1)

! GET EIGENVECTORS OF TRIDIAGONALIZED VERSION OF A.

    IF(NEVEC1 <= 0) RETURN
    CALL EIGVEC(NDIM,NEVEC1,A,V,TOLTMP,ANORM,EVAL1,IDEGEN1,EVEC1)

! PREMULTIPLY EVEC1 BY THE HOUSEHOLDER MATRIX USED TO TRIDIAGONALIZE
! A.  THIS TRANSFORMS EIGENVECTORS OF THE TRIDIAGONAL MATRIX TO
! THOSE OF THE ORIGINAL MATRIX A.  SEE SUBROUTINE TRIDI FOR
! STORAGE OF HOUSEHOLDER TRANSFORMATION.

    IF(NDIM > 2)THEN
        DO 30 K=1,NDIM-2
            V(1,K) = A(K+2,K)
        30 ENDDO
    
    ! SWAP STORAGE SO THAT THE EXPENSIVE TRIPLE LOOP BELOW DOESN'T
    ! HAVE TO JUMP ACROSS COLUMNS TO GET ENTRIES OF A.
    
        DO 50 I=2,NDIM
            DO 40 J=1,I-1
                A(I,J) = A(J,I)
            40 ENDDO
        50 ENDDO
        DO 160 J=1,NEVEC1
            DO 140 M=2,NDIM-1
                K = NDIM - M
                BETAH = V(1,K)
!                IF(ABS(BETAH) < 1.0D-50) GO TO 140
                IF(ABS(BETAH) < 1.0D-50) CONTINUE
                SUM = 0.0D0
                DO 100 I=K+1,NDIM
                    SUM = SUM + A(I,K)*EVEC1(I,J)
                100 ENDDO
                BSUM = BETAH*SUM
                DO 120 I=K+1,NDIM
                    EVEC1(I,J) = EVEC1(I,J) - A(I,K)*BSUM
                120 ENDDO
            140 ENDDO
        160 ENDDO
    ENDIF
    RETURN
    end SUBROUTINE DIAG


!********************************************************
! TRIDI
!------------------------------------------------------------

    SUBROUTINE TRIDI(NDIM,V,A)

! TRIDIAGONALIZES A REAL, SYMMETRIC MATRIX A BY THE METHOD OF
! HOUSEHOLDER (J. H. WILKINSON, THE COMPUTER JOURNAL, VOL. 3,
! P. 23 (1960)).  NDIM IS THE ORDER OF A.  THE DIAGONAL AND
! SUBDIAGONAL OF A ARE OVERWRITTEN WITH THE TRIDIAGONALIZED
! VERSION OF A.  THE VECTORS USED IN EACH HOUSEHOLDER
! TRANSFORMATION ARE STORED ABOVE THE DIAGONAL IN THE FIRST
! NDIM-2 ROWS OF A.  THE BETAHS ARE RETURNED BELOW THE SUBDIAGONAL
! OF A.  V IS A WORKSPACE ARRAY.

! PROGRAMMED BY S. L. DIXON, OCT., 1991.


    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION A(NDIM,*),V(3,*)
    DIMENSION A(nbasis,nbasis),V(3,nbasis)

! THRESH WILL BE USED AS A THRESHOLD TO DETERMINE IF A VALUE SHOULD
! BE CONSIDERED TO BE ZERO.  THIS CAN BE CHANGED BY THE USER.

    THRESH = 1.0D-50

! IF A IS 2 BY 2 OR SMALLER, THEN IT IS ALREADY TRIDIAGONAL -- NO
! NEED TO CONTINUE.

    IF(NDIM <= 2) GO TO 1000
    DO 500 K=1,NDIM-2
    
    ! DETERMINE THE VECTOR V USED IN THE HOUSEHOLDER TRANSFORMATION P.
    ! FOR EACH VALUE OF K THE HOUSEHOLDER MATRIX P IS DEFINED AS:
    
    ! P = I - BETAH*V*V'
    
    
    ! CONSTRUCT A HOUSEHOLDER TRANSFORMATION ONLY IF THERE IS A NONZERO
    ! OFF-DIAGONAL ELEMENT BELOW A(K,K).
    
        ALPHA2 = 0.0D0
        DO 60 I=K+1,NDIM
            V(1,I) = A(I,K)
            ALPHA2 = ALPHA2 + V(1,I)**2
        60 ENDDO
        APTEMP = ALPHA2 - V(1,K+1)**2
        ALPHA = DSQRT(ALPHA2)
        IF(ALPHA >= THRESH)THEN
            BETAH = 1.0D0/(ALPHA*(ALPHA + ABS(V(1,K+1))))
            SGN = SIGN(1.0D0,V(1,K+1))
            V(1,K+1) = V(1,K+1) + SGN*ALPHA
        
        ! NOW OVERWRITE A WITH P'*A*P.  THE ENTRIES BELOW THE SUBDIAGONAL
        ! IN THE KTH COLUMN ARE ZEROED BY THE PREMULTIPLICATION BY P'.
        ! THESE ENTRIES WILL BE LEFT ALONE TO SAVE TIME.
        
            AKV = APTEMP + A(K+1,K)*V(1,K+1)
            S = BETAH*AKV
            A(K+1,K) = A(K+1,K) - S*V(1,K+1)
        
        ! NOW THE SUBMATRIX CONSISTING OF ROWS K+1,NDIM AND COLUMNS K+1,NDIM
        ! MUST BE OVERWRITTEN WITH THE TRANSFORMATION.
        
            DOT12 = 0.0D0
            BHALF = BETAH*0.5D0
            DO 220 I=K+1,NDIM
                SUM = 0.0D0
                DO 100 J=K+1,I
                    SUM = SUM + A(I,J)*V(1,J)
                100 ENDDO
                IF(I < NDIM)THEN
                    DO 180 J=I+1,NDIM
                    
                    ! AN UPPER TRIANGULAR ENTRY OF A WILL BE REQUIRED.  MUST USE
                    ! THE SYMMETRIC ENTRY IN THE LOWER TRIANGULAR PART OF A.
                    
                        SUM = SUM + A(J,I)*V(1,J)
                    180 ENDDO
                ENDIF
                V(2,I) = BETAH*SUM
                DOT12 = DOT12 + V(1,I)*V(2,I)
            220 ENDDO
            BH12 = BHALF*DOT12
            DO 300 I=K+1,NDIM
                V(2,I) = V(2,I) - BH12*V(1,I)
            300 ENDDO
            DO 350 J=K+1,NDIM
                DO 310 I=J,NDIM
                    A(I,J) = A(I,J) - V(1,I)*V(2,J) - V(2,I)*V(1,J)
                310 ENDDO
            
            ! STORE V(1,J) ABOVE THE DIAGONAL IN ROW K OF A
            
                A(K,J) = V(1,J)
            350 ENDDO
        
        ! STORE BETAH BELOW THE SUBDIAGONAL OF A.
        
            A(K+2,K) = BETAH
        ELSE
        
        ! NO HOUSEHOLDER TRANSFORMATION IS NECESSARY BECAUSE THE OFF-
        ! DIAGONALS ARE ALL ESSENTIALLY ZERO.
        
            A(K+2,K) = 0.0D0
            DO 460 J=K+1,NDIM
                A(K,J) = 0.0D0
            460 ENDDO
        ENDIF
    500 ENDDO
    1000 RETURN
    end SUBROUTINE TRIDI


!********************************************************
! EIGVAL
!------------------------------------------------------------

    SUBROUTINE EIGVAL(NDIM,A,BETAH,TOLERA,ANORM,EVAL1,IERROR)

! QR ROUTINE FOR THE DETERMINATION OF ALL THE EIGENVALUES
! OF THE NDIM BY NDIM SYMMETRIC, TRIDIAGONAL MATRIX A.

! INPUT:

! NDIM   = SIZE OF MATRIX A.
! A      = NDIM BY NDIM SYMMETRIC TRIDIAGONAL MATRIX.
! BETAH   = 3 BY NDIM WORKSPACE.
! TOLERA = SMALL NUMBER USED TO DETERMINE WHEN OFF-DIAGONAL
! ELEMENTS ARE ESSENTIALLY ZERO.
! ANORM  = ABSOLUTE COLUMN NORM OF TRIDIAGONAL MATRIX A.


! RETURNED:

! EVAL1   = EIGENVALUES OF A IN ASCENDING ORDER.
! IERROR = 1 IF QR ITERATION DID NOT CONVERGE; 0 OTHERWISE.

! PROGRAMMED BY S. L. DIXON.


    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION A(NDIM,*),BETAH(3,*),EVAL1(*)
    DIMENSION A(nbasis,nbasis),BETAH(3,nbasis),EVAL1(nbasis)
    IERROR = 0
! ITMAX = 20
    ITMAX = 200

! TOLERANCE FOR OFF-DIAGONAL ELEMENTS:

    EPSLON = TOLERA*ANORM

! COPY DIAGONAL ELEMENTS OF A TO EVAL1, AND SUBDIAGONAL ELEMENTS
! TO BETAH.

    EVAL1(1) = A(1,1)
    BETAH(1,1) = A(2,1)
    IF(NDIM > 2)THEN
        DO 50 I=2,NDIM-1
            EVAL1(I) = A(I,I)
            BETAH(1,I) = A(I+1,I)
        50 ENDDO
    ENDIF
    EVAL1(NDIM) = A(NDIM,NDIM)

! EACH QR ITERATION WILL OPERATE ON THE UNREDUCED TRIDIAGONAL
! SUBMATRIX WITH UPPER LEFT ELEMENT (L,L) AND LOWER RIGHT ELEMENT
! (N,N).

    L = 1
    N = NDIM
    ITER = 0

! FIND THE SMALLEST UNREDUCED SUBMATRIX WITH LOWER RIGHT CORNER AT
! (N,N).  I.E., SEARCH UPWARD FOR A BETAH THAT IS ZERO.

    80 KUPPER = N-L
    DO 100 K=1,KUPPER
        I = N-K
        IF(ABS(BETAH(1,I)) <= EPSLON)THEN
            L = I+1
            GO TO 150
        ENDIF
    100 ENDDO

! IF WE GET TO THE NEXT STATEMENT, THEN THERE ARE NO ZERO OFF-DIAGONALS
! FOR THE SUBMATRIX WITH UPPER LEFT A(L,L) AND LOWER RIGHT A(N,N).
! WE CAN STILL GET EIGENVALUES IF THE MATRIX IS 2 BY 2 OR 1 BY 1.
! OTHERWISE, DO ANOTHER QR ITERATION PROVIDED ITMAX CYCLES HAVE
! NOT OCCURRED.

    IF(L == N .OR. L == N-1)THEN
        GO TO 150
    ELSE
        IF(ITER == ITMAX)THEN
            IERROR = 1
            GO TO 1000
        ELSE
            GO TO 200
        ENDIF
    ENDIF

! IF WE GET TO 150 THEN A(L,L-1) IS ZERO AND THE UNREDUCED SUBMATRIX
! HAS UPPER LEFT AT A(L,L) AND LOWER RIGHT AT A(N,N).  WE CAN
! EXTRACT ONE EIGENVALUE IF THIS MATRIX IS 1 BY 1 AND 2 EIGENVALUES
! IF IT IS 2 BY 2.

    150 IF(L == N)THEN
    
    ! IT'S A 1 BY 1 AND EVAL1(N) IS AN EIGENVALUE.  IF L=2 OR 1 WE ARE
    ! DONE.  OTHERWISE, UPDATE N, RESET L AND ITER, AND REPEAT THE
    ! SEARCH.
    
        IF(L <= 2)THEN
            GO TO 500
        ELSE
            N = L-1
            L = 1
            ITER = 0
            GO TO 80
        ENDIF
    ELSEIF(L == N-1)THEN
    
    ! THE UNREDUCED SUBMATRIX IS A 2 BY 2.  OVERWRITE EVAL1(N-1)
    ! AND EVAL1(N) WITH THE EIGENVALUES OF THE LOWER RIGHT 2 BY 2.
    
        BTERM = EVAL1(N-1) + EVAL1(N)
        ROOT1 = BTERM*0.5D0
        ROOT2 = ROOT1
        DISCR = BTERM**2 - 4.0D0*(EVAL1(N-1)*EVAL1(N)-BETAH(1,N-1)**2)
        IF(DISCR > 0.0D0)THEN
            D = DSQRT(DISCR)*0.5D0
            ROOT1 = ROOT1 - D
            ROOT2 = ROOT2 + D
        ENDIF
        EVAL1(N-1) = ROOT1
        EVAL1(N) = ROOT2
    
    ! SEE IF WE ARE DONE.  IF NOT, RESET N, L, AND ITER AND LOOK
    ! FOR NEXT UNREDUCED SUBMATRIX.
    
        IF(L <= 2)THEN
            GO TO 500
        ELSE
            N = L-1
            L = 1
            ITER = 0
            GO TO 80
        ENDIF
    ELSE
    
    ! AN EIGENVALUE WAS FOUND AND THE NEW UNREDUCED MATRIX LIMITS
    ! N AND L ARE SET.  DO A QR ITERATION ON NEW MATRIX.
    
        ITER = 0
        GO TO 200
    ENDIF

! QR ITERATION BEGINS HERE.

    200 ITER = ITER + 1

! USE EIGENVALUES OF THE LOWER RIGHT 2 BY 2 TO COMPUTE SHIFT.  SHIFT
! BY THE EIGENVALUE CLOSEST TO EVAL1(N).

    D = (EVAL1(N-1) - EVAL1(N))*0.5D0
    SIGND = 1.0D0
    IF(D < 0.0D0) SIGND = -1.0D0
    SHIFT = EVAL1(N) + D - SIGND*DSQRT(D*D + BETAH(1,N-1)**2)
    P = EVAL1(L) - SHIFT
    R = BETAH(1,L)
    T = EVAL1(L)
    W = BETAH(1,L)

! OVERWRITE A WITH Q'*A*Q.

    DO 250 K=L,N-1
        D = DSQRT(P*P + R*R)
        C = P/D
        S = R/D
        IF(K /= L) BETAH(1,K-1) = D
        CC = C*C
        SS = 1.0D0 - CC
        CS = C*S
        CSW = 2.0D0*CS*W
        AK1 = EVAL1(K+1)
        EVAL1(K) = CC*T + CSW + SS*AK1
        P = (CC - SS)*W + CS*(AK1 - T)
        T = SS*T - CSW + CC*AK1
        R = S*BETAH(1,K+1)
        W = C*BETAH(1,K+1)
    250 ENDDO
    BETAH(1,N-1) = P
    EVAL1(N) = T

! GO BACK AND SEE IF L AND N NEED TO BE UPDATED.

    GO TO 80

! SORT EIGENVALUES IN ASCENDING ALGEBRAIC ORDER.

    500 DO 600 I=2,NDIM
        JMAX = NDIM-I+1
        ISORT = 0
        DO 550 J=1,JMAX
            IF(EVAL1(J) > EVAL1(J+1))THEN
                ETEMP = EVAL1(J)
                EVAL1(J) = EVAL1(J+1)
                EVAL1(J+1) = ETEMP
                ISORT = 1
            ENDIF
        550 ENDDO
        IF(ISORT == 0) GO TO 1000
    600 ENDDO
    1000 RETURN
    end SUBROUTINE EIGVAL

!********************************************************
! DEGEN
!------------------------------------------------------------
    SUBROUTINE DEGEN(NDIM,EVAL1,TOLERA,ANORM,IDEGEN1)

! DETERMINES DEGENERACIES OF THE EIGENVALUES.

! INPUT:

! NDIM   = SIZE OF MATRIX BEING DIAGONALIZED.
! EVAL1   = SORTED EIGENVALUES (INCREASING VALUE).
! TOLERA = SAME TOLERANCE USED TO DETERMINE EIGENVALUES.
! ANORM  = ABSOLUTE COLUMN NORM OF TRIDIAGONAL MATRIX.


! RETURNED:

! IDEGEN1 = DEGENERACIES OF EIGENVALUES.

    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION EVAL1(*),IDEGEN1(*)
    DIMENSION EVAL1(nbasis),IDEGEN1(nbasis)

! DETERMINE DEGENERACIES OF EIGENVALUES.  ADJACENT EIGENVALUES
! WILL BE CONSIDERED TO BE DEGENERATE WHEN THEY DIFFER BY LESS
! THAN DTOLER.

    DTOLER = MAX(ANORM*DSQRT(TOLERA),1.0D-8)
    NSAME = 1
    DO 200 I=2,NDIM
        DIFF = ABS(EVAL1(I-1) - EVAL1(I))
        IF(DIFF <= DTOLER)THEN
        
        ! EIGENVALUES I-1 AND I ARE DEGENERATE.
        
            NSAME = NSAME + 1
            IF(I == NDIM)THEN
            
            ! WE'VE COME TO THE LAST REQUESTED EIGENVALUE, AND IT'S TIME
            ! TO ASSIGN DEGENERACIES FOR THE BLOCK ENDING WITH THE ITH
            ! EIGENVALUE.
            
                DO 100 J=I-NSAME+1,I
                    IDEGEN1(J) = NSAME
                100 ENDDO
            ENDIF
        
        ! GO TO THE NEXT EIGENVALUE (IF THERE ARE ANY LEFT) AND SEE IF
        ! IT'S DEGENERATE WITH THE NSAME EIGENVALUES WE'VE ALREADY
        ! FOUND.
        
            GO TO 200
        ELSE
        
        ! EITHER EIGENVALUE I-1 IS NONDEGENERATE OR IT'S THE LAST
        ! EIGENVALUE IN A DEGENERATE BLOCK.  CORRESPONDINGLY, ASSIGN THE
        ! PROPER DEGENERACY TO I-1 OR TO EACH EIGENVALUE IN THE BLOCK.
        
            DO 150 J=I-NSAME,I-1
                IDEGEN1(J) = NSAME
            150 ENDDO
            NSAME = 1
        
        ! IF I=NDIM THEN IT MUST BE THE CASE THAT THIS LAST EIGENVALUE
        ! IS NONDEGENERATE.
        
            IF(I == NDIM) IDEGEN1(I) = 1
        ENDIF
    200 ENDDO
    RETURN
    end SUBROUTINE DEGEN


!********************************************************
! EIGVEC
!------------------------------------------------------------
    SUBROUTINE EIGVEC(NDIM,NEVEC1,A,AWORK,TOLERA,ANORM,EVAL1,IDEGEN1, &
    EVEC1)

! INVERSE ITERATION ROUTINE FOR EIGENVECTOR DETERMINATION.
! CALCULATES THE EIGENVECTORS OF AN NDIM BY NDIM SYMMETRIC,
! TRIDIAGONAL MATRIX A.

! INPUT:

! NDIM   = SIZE OF MATRIX A.
! NEVEC1  = NUMBER OF EIGENVECTORS REQUIRED.
! A      = NDIM BY NDIM TRIDIAGONAL MATRIX.
! TOLERA = SAME TOLERANCE USED TO DETERMINE EIGENVALUES.
! ANORM  = ABSOLUTE COLUMN NORM OF TRIDIAGONAL MATRIX A.
! EVAL1   = SORTED EIGENVALUES OF A.
! IDEGEN1 = DEGENERACIES OF EIGENVALUES.


! RETURNED:

! EVEC1   = EIGENVECTORS OF TRIDIAGONAL MATRIX (IN COLUMNS).

! PROGRAMMED BY S. L. DIXON, OCT., 1991.


    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION A(NDIM,*),AWORK(3,*),EVAL1(*),IDEGEN1(*),EVEC1(NDIM,*)
    DIMENSION A(nbasis,nbasis),AWORK(3,nbasis),EVAL1(nbasis), &
    IDEGEN1(nbasis),EVEC1(nbasis,nbasis)
    LOGICAL :: ORTH
    IRAND = 13876532

! COMPUTE THRESHOLD EPSLON WHICH WILL BE USED IF THE INVERSE ITERATION
! MATRIX IS SINGULAR.

    EPSLON = ANORM*TOLERA

! WHEN DEGENERACIES OCCUR, THERE ARE RARE INSTANCES WHEN THE
! DEGENERATE BLOCK OF EIGENVECTORS ARE NOT LINEARLY INDEPENDENT.
! IN THESE CASES, AN ADDITIONAL PASS THROUGH THE INVERSE ITERATION
! (WITH A NEW SET OF RANDOM NUMBERS) IS CARRIED OUT, I.E., CONTROL
! PASSES TO STATEMENT 40.  IVECT WILL KEEP TRACK OF THE CURRENT
! STARTING EIGENVECTOR WHEN ADDITIONAL PASSES ARE NECESSARY.

    IVECT = 1
    NFAIL = 0
    NPRTRB = 0

! DO ONE ITERATION FOR EACH EIGENVECTOR.

    40 ORTH = .TRUE.
    NDEGEN = 0
    MSTART = IVECT
    DO 380 M=MSTART,NEVEC1
        IF(IDEGEN1(M) > 1) NDEGEN = NDEGEN + 1
        Z = EVAL1(M)
    
    ! IF THE INVERSE ITERATION HAS FAILED TWICE DUE TO NON-ORTHOGONALITY
    ! OF DEGENERATE EIGENVECTORS, PERTURB THE EIGENVALUE BY A SMALL
    ! AMOUNT.
    
        IF(NFAIL >= 2)THEN
            NPRTRB = NPRTRB + 1
            Z = Z + 0.001D0*dble(NPRTRB)*TOLERA
        ENDIF
    
    ! STORE THE TRIDIAGONAL ENTRIES OF THE INVERSE ITERATION MATRIX IN
    ! THE 3 BY NDIM WORKSPACE AWORK.
    
        AWORK(1,1) = 0.0D0
        AWORK(2,1) = A(1,1) - Z
        AWORK(3,1) = A(2,1)
        IF(NDIM > 2)THEN
            DO 80 I=2,NDIM-1
                AWORK(1,I) = A(I,I-1)
                AWORK(2,I) = A(I,I) - Z
                AWORK(3,I) = A(I+1,I)
            80 ENDDO
        ENDIF
        AWORK(1,NDIM) = A(NDIM,NDIM-1)
        AWORK(2,NDIM) = A(NDIM,NDIM) - Z
        AWORK(3,NDIM) = 0.0D0
    
    ! ASSIGN INVERSE ITERATION VECTOR FROM RANDOM NUMBERS.
    
        DO 120 I=1,NDIM
            CALL RANDOM(IRAND,RNDOM)
            RNDOM = 2.0D0*(RNDOM - 0.5D0)
            EVEC1(I,M) = RNDOM*TOLERA
        120 ENDDO
    
    ! CARRY OUT FORWARD GAUSSIAN ELIMINATION WITH ROW PIVOTING
    ! ON THE INVERSE ITERATION MATRIX.
    
        DO 160 K=1,NDIM-1
            ADIAG = ABS(AWORK(2,K))
            ASUB = ABS(AWORK(1,K+1))
            IF(ADIAG >= ASUB)THEN
            
            ! USE PIVOTAL ELEMENT FROM ROW K.
            
                IF(AWORK(2,K) == 0.0D0) AWORK(2,K) = EPSLON
                T = AWORK(1,K+1)/AWORK(2,K)
                AWORK(1,K+1) = 0.0D0
                AWORK(2,K+1) = AWORK(2,K+1) - T*AWORK(3,K)
            
            ! LEFT-JUSTIFY EQUATION K SO THAT DIAGONAL ENTRY IS STORED
            ! IN AWORK(1,K).
            
                AWORK(1,K) = AWORK(2,K)
                AWORK(2,K) = AWORK(3,K)
                AWORK(3,K) = 0.0D0
            
            ! OPERATE ON VECTOR AS WELL.
            
                EVEC1(K+1,M) = EVEC1(K+1,M) - T*EVEC1(K,M)
            ELSE
            
            ! USE PIVOTAL ELEMENT FROM ROW K+1 AND SWAP ROWS K AND K+1.
            
                IF(AWORK(1,K+1) == 0.0D0) AWORK(1,K+1) = EPSLON
                T = AWORK(2,K)/AWORK(1,K+1)
                ATEMP = AWORK(3,K) - T*AWORK(2,K+1)
                AWORK(1,K) = AWORK(1,K+1)
                AWORK(2,K) = AWORK(2,K+1)
                AWORK(3,K) = AWORK(3,K+1)
                AWORK(1,K+1) = 0.0D0
                AWORK(2,K+1) = ATEMP
                AWORK(3,K+1) = -T*AWORK(3,K+1)
            
            ! OPERATE ON VECTOR AND SWAP ENTRIES.
            
                ETEMP = EVEC1(K+1,M)
                EVEC1(K+1,M) = EVEC1(K,M) - ETEMP*T
                EVEC1(K,M) = ETEMP
            ENDIF
        160 ENDDO
    
    ! FORWARD ELIMINATION COMPLETE.  BACK SUBSTITUTE TO GET SOLUTION.
    ! OVERWRITE COLUMN M OF EVEC1 WITH SOLUTION.
    
        IF(AWORK(2,NDIM) == 0.0D0) AWORK(2,NDIM) = EPSLON
        EVEC1(NDIM,M) = EVEC1(NDIM,M)/AWORK(2,NDIM)
        ETEMP = EVEC1(NDIM-1,M) - AWORK(2,NDIM-1)*EVEC1(NDIM,M)
        EVEC1(NDIM-1,M) = ETEMP/AWORK(1,NDIM-1)
        ENORM = EVEC1(NDIM,M)**2 + EVEC1(NDIM-1,M)**2
        IF(NDIM > 2)THEN
        
        ! CAUTION: PROBLEM LOOP FOR SOME IBM RS/6000 COMPILERS.  VALUE
        ! OF K CAN GET LOST WHEN OPTIMIZE FLAG IS USED.
        
            DO 200 L=1,NDIM-2
                K = NDIM-L-1
                ETEMP = EVEC1(K,M) - AWORK(2,K)*EVEC1(K+1,M) &
                - AWORK(3,K)*EVEC1(K+2,M)
                EVEC1(K,M) = ETEMP/AWORK(1,K)
                ENORM = ENORM + EVEC1(K,M)**2
            200 ENDDO
        ENDIF
        EINV = 1.0D0/DSQRT(ENORM)
    
    ! NORMALIZE EIGENVECTOR.
    
        DO 240 I=1,NDIM
            EVEC1(I,M) = EVEC1(I,M)*EINV
        240 ENDDO
    
    ! IF WE HAVE COME TO THE END OF A DEGENERATE BLOCK OF EIGENVECTORS,
    ! ORTHOGONALIZE THE BLOCK.
    
        IF(NDEGEN > 1)THEN
            IF(NDEGEN == IDEGEN1(M) .OR. M == NEVEC1)THEN
                JSTART = M-NDEGEN+1
                CALL ORTHOG(NDIM,NDEGEN,JSTART,EVEC1,ORTH)
                IF(ORTH)THEN
                    NFAIL = 0
                    NPRTRB = 0
                
                ! THE DEGENERATE VECTORS WERE LINEARLY INDEPENDENT AND WERE
                ! SUCCESSFULLY ORTHOGONALIZED.
                
                    IVECT = IVECT + NDEGEN
                    NDEGEN = 0
                ELSE
                
                ! THE BLOCK IS APPARENTLY NOT LINEARLY INDEPENDENT.  GO BACK
                ! AND REPEAT THE INVERSE ITERATION FOR THESE VECTORS.  AFTER
                ! AN INDEPENDENT SET HAS BEEN FOUND, ANY ADDITIONAL EIGENVECTORS
                ! WILL BE DETERMINED.
                
                    NFAIL = NFAIL + 1
                    GO TO 40
                ENDIF
            ENDIF
        ENDIF
    
    ! THE CURRENT EIGENVECTOR SHOULD BE OKAY IF IT IS NONDEGENERATE.
    
        IF(IDEGEN1(M) == 1) IVECT = IVECT + 1
    380 ENDDO
    RETURN
    end SUBROUTINE EIGVEC

!********************************************************
! ORTHOG
!------------------------------------------------------------

    SUBROUTINE ORTHOG(NDIM,NVECT,JSTART,VECT,ORTH)

! CONSTRUCTS A SET OF ORTHONORMAL VECTORS FROM THE NVECT LINEARLY
! INDEPENDENT, NORMALIZED VECTORS IN THE ARRAY VECT.  THE VECTORS
! SHOULD BE STORED COLUMNWISE, STARTING IN COLUMN JSTART.  VECT IS
! OVERWRITTEN WITH THE ORTHONORMAL SET.  ALL VECTORS ARE NDIM BY 1.
! ORTH IS RETURNED WITH A VALUE OF .TRUE. IF THE SET WAS LINEARLY
! INDEPENDENT AND .FALSE. OTHERWISE.

! PROGRAMMED BY S. L. DIXON.


    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION VECT(NDIM,*)
    DIMENSION VECT(nbasis,nbasis)
    LOGICAL :: ORTH

    ORTH = .TRUE.
    ORTEST = 1.0D-8

! BEGIN ORTHOGONALIZATION.

    JSTOP = JSTART + NVECT - 1
    DO 120 J=JSTART,JSTOP
        IF(J > JSTART)THEN
        
        ! SUBTRACT OFF COMPONENTS OF PREVIOUSLY DETERMINED ORTHOGONAL
        ! VECTORS FROM THE VECTOR IN COLUMN J.
        
            DO 60 JPREV=JSTART,J-1
                DOT = 0.0D0
                DO 20 I=1,NDIM
                    DOT = DOT + VECT(I,JPREV)*VECT(I,J)
                20 ENDDO
                DO 40 I=1,NDIM
                    VECT(I,J) = VECT(I,J) - DOT*VECT(I,JPREV)
                40 ENDDO
            60 ENDDO
        ENDIF
    
    ! NORMALIZE COLUMN J.
    
        VJNORM = 0.0D0
        DO 80 I=1,NDIM
            VJNORM = VJNORM + VECT(I,J)**2
        80 ENDDO
        VJNORM = DSQRT(VJNORM)
    
    ! IF THE NORM OF THIS VECTOR IS TOO SMALL THEN THE VECTORS ARE
    ! NOT LINEARLY INDEPENDENT.
    
        IF(VJNORM < ORTEST)THEN
            ORTH = .FALSE.
            GO TO 1000
        ENDIF
        DO 100 I=1,NDIM
            VECT(I,J) = VECT(I,J)/VJNORM
        100 ENDDO
    120 ENDDO
    1000 RETURN
    end SUBROUTINE ORTHOG

!********************************************************
! RANDOM
!------------------------------------------------------------

    SUBROUTINE RANDOM(IX,X)

! GENERATES INTEGER*4 AND DOUBLE PRECISION RANDOM NUMBERS ON
! THE INTERVALS:

! 0 < IX < (2**31)-1
! AND      0.0 < X < 1.0

! ON THE FIRST CALL IX SHOULD SATISFY THE TOP INEQUALITY.

! NUMBERS ARE GENERATED USING THE RELATION,

! IX = IX*IC (MODULO (2**31)-1),  WHERE IC = 7**5


    double precision :: X


! INITIALIZE:     I15 = 2**15
! I16 = 2**16
! I31_1 = 2**31 - 1
! IC = 7**5

    DATA I15 /32768/
    DATA I16 /65536/
    DATA I31_1 /2147483647/
    DATA IC /16807/

    SAVE I15,I16,I31_1,IC

    IX16 = IX/I16
    I16RMD = IX - IX16*I16

! NOTE THAT IX = IX16*I16 + I16RMD    (I16RMD = 16-BIT REMAINDER)

    IX16IC = IX16*IC
    IXIC31 = IX16IC/I15

! NOTE:   IX*IC = (IX16*I16 + I16RMD)*IC
! = (IX16*IC*I16) + (I16RMD*IC)
! = (IX16IC*I16 ) + (I16RMD*IC)
! = (  TERM 1   ) + ( TERM 2  )
! AND,

! IX16IC = ((IX16IC/I15)*I15)  +  (IX16IC - (IX16IC/I15)*I15))
! = (  IXIC31  )*I15)  +  (IX16IC - (  IXIC31  )*I15 )
! (     15-BIT REMAINDER     )

! THEREFORE,  TERM 1 = ((IXIC31*I15) + (IX16IC - IXIC31*I15))*I16

! THEN,

! (   TERM A     )   (        TERM B          )   ( TERM C  )
! IX*IC = ((IXIC31*I16*I15) + (IX16IC - IXIC31*I15)*I16) + (I16RMD*IC)
! = (                  TERM 1                    ) + ( TERM 2  )


! NOTE THAT TERM B AND TERM C ARE BOTH LESS THAN 2**31 - 1.  ONLY
! TERM A HAS THE POSSIBILITY OF EXCEEDING 2**31 - 1.  BUT SINCE
! I16*I15 = 2**31, THE FACTOR IXIC31 INDICATES EXACTLY HOW MANY TIMES
! TERM A "WRAPS" AROUND THE INTERVAL (0,2**31 - 1).  THUS, FOR THE
! MODULO OPERATION, TERM A MAY BE REPLACED BY IXIC31.  THE SUM OF
! TERM A AND TERM B MIGHT EXCEED 2**31 - 1, BUT WE CAN SUBSTRACT
! 2**31 - 1 FROM ONE OF THEM TO PREVENT THIS FROM HAPPENING.

    IX = IXIC31 + ((IX16IC-IXIC31*I15)*I16 - I31_1) + I16RMD*IC

! ADD I31_1 BACK IN IF THE SUBTRACTION MADE IX NEGATIVE.

    IF(IX < 0) IX = IX + I31_1

! MAKE X RANDOM ON (0.0,1.0) BY MULTIPLYING IX BY 1.0/I31_1

    X = dble(IX)*4.6566128752458D-10
    RETURN
    end SUBROUTINE RANDOM
	
!
!PrtLab(line,kk,ktmp) written by Li Wei	
!****** 2005.01.09 write Ktmp(kk) into line: eg. "1,3,5-8,2*0,9"
!
      subroutine PrtLab(line,kk,Ktmp)
      integer Ktmp(kk)
      character ch,ch2,line*(*),line1*100,line2*100
      parameter (ch=',',ch2='-')

      line=' '; ini=1; fc=1

      write(line1,*) Ktmp(1)
      call EffChar(line1,1,100,k1,k2)
      line(ini:ini+k2-k1)=line1(k1:k2)
      ini=ini+k2-k1+1
      nz=0
      if (Ktmp(1)==0.and.Ktmp(2)==0) then
         fc=0; nz=1; ini=1
      endif

      do 110 i=2,kk
         write(line1,*) Ktmp(i)
         call EffChar(line1,1,100,k1,k2)
         if (Ktmp(i)-Ktmp(i-1)==1) then
            if (i==kk.or.Ktmp(i+1)-Ktmp(i).ne.1) then
               line(ini:ini+k2-k1+1)=ch2//line1(k1:k2)
               ini=ini+k2-k1+2
            endif
         elseif (Ktmp(i)==0) then
            nz=nz+1
            if (i==kk) then
               if (nz==1) then
                  line(ini:ini+1)=ch//'0'; ini=ini+2
               else
                  write(line2,*) nz; call EffChar(line2,1,100,k3,k4)
                  if (fc==1) then  ! 2005.01.09 add
                     line(ini:ini+k4-k3+3)=ch//line2(k3:k4)//'*0'
                     ini=ini+k4-k3+4; exit
                  else
                     line(ini:ini+k4-k3+2)=line2(k3:k4)//'*0'
                     ini=ini+k4-k3+3; fc=1; exit
                  endif
               endif
            else
               if (Ktmp(i+1).ne.0.and.nz==1) then
                  line(ini:ini+1)=ch//'0'; ini=ini+2; nz=0
               elseif (Ktmp(i+1).ne.0.and.nz.ne.1) then
                  write(line2,*) nz; call EffChar(line2,1,100,k3,k4)
                  if (fc==1) then  ! 2005.01.09 add
                     line(ini:ini+k4-k3+3)=ch//line2(k3:k4)//'*0'
                     ini=ini+k4-k3+4; nz=0
                  else
                     line(ini:ini+k4-k3+2)=line2(k3:k4)//'*0'
                     ini=ini+k4-k3+3; nz=0; fc=1
                  endif
               endif
            endif
         else
            line(ini:ini+k2-k1+1)=ch//line1(k1:k2)
            ini=ini+k2-k1+2
         endif
 110  enddo

      end

!****** 2005.01.07 move blank of two sides in a line ******
      subroutine EffChar(line,ini,ifi,k1,k2)
      implicit none
      integer ini,ifi,k1,k2,i,j
      character line*(*)

      do i=ini,ifi
         if (line(i:i).ne.' ') then
            k1=i; exit
         endif
      enddo

      do i=ifi,ini,-1
         if (line(i:i).ne.' ') then
            k2=i; exit
         endif
      enddo

      end

! Ed Brothers. October 12, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    double precision function ekinetic(a,b,i,j,k,ii,jj,kk,Ax,Ay,Az,Bx, &
    By,Bz)
    implicit double precision(a-h,o-z)
    double precision :: kinetic

! The purpose of this subroutine is to calculate the kinetic energy
! of an electron  distributed between gtfs with orbital exponents a
! and b on A and B with angular momentums defined by i,j,k (a's x, y
! and z exponents, respectively) and ii,jj,and kk on B.

! The first step is to see if this function is zero due to symmetry.
! If it is not, reset kinetic to 0.

    kinetic = (1+(-1)**(i+ii))*(1+(-1)**(j+jj))*(1+(-1)**(k+kk)) &
    +(Ax-Bx)**2 + (Ay-By)**2 + (Az-Bz)**2
    IF (kinetic == 0.d0) goto 100
    kinetic=0.d0

! Kinetic energy is the integral of an orbital times the second derivative
! over space of the other orbital.  For GTFs, this means that it is just a
! sum of various overlap integrals with the powers adjusted.

    xi = dble(i)
    xj = dble(j)
    xk = dble(k)
    kinetic = kinetic + (-1.d0+xi)*xi*overlap(a,b,i-2,j,k,ii,jj,kk, &
    Ax,Ay,Az,Bx,By,Bz) &
    - 2.d0*a*(1.d0+2.d0*xi)*overlap(a,b,i,j,k,ii,jj,kk, &
    Ax,Ay,Az,Bx,By,Bz) &
    + 4.d0*(a**2.d0)*overlap(a,b,i+2,j,k,ii,jj,kk, &
    Ax,Ay,Az,Bx,By,Bz)
    kinetic = kinetic + (-1.d0+xj)*xj*overlap(a,b,i,j-2,k,ii,jj,kk, &
    Ax,Ay,Az,Bx,By,Bz) &
    - 2.d0*a*(1.d0+2.d0*xj)*overlap(a,b,i,j,k,ii,jj,kk, &
    Ax,Ay,Az,Bx,By,Bz) &
    + 4.d0*(a**2.d0)*overlap(a,b,i,j+2,k,ii,jj,kk, &
    Ax,Ay,Az,Bx,By,Bz)
    kinetic = kinetic + (-1.d0+xk)*xk*overlap(a,b,i,j,k-2,ii,jj,kk, &
    Ax,Ay,Az,Bx,By,Bz) &
    - 2.d0*a*(1.d0+2.d0*xk)*overlap(a,b,i,j,k,ii,jj,kk, &
    Ax,Ay,Az,Bx,By,Bz) &
    + 4.d0*(a**2.d0)*overlap(a,b,i,j,k+2,ii,jj,kk, &
    Ax,Ay,Az,Bx,By,Bz)

    100 continue
    ekinetic = kinetic/(-2.d0)
    return
    end function ekinetic
	
	    SUBROUTINE LSOLVE(N,ISIZE,A,B,W,THRESH,X,IERROR)

! USES GAUSSIAN ELIMINATION WITH ROW PIVOTING TO SOLVE THE ORDER N
! LINEAR SYSTEM:

! A*X = B.

! W IS A WORK VECTOR OF LENGTH N, AND THRESH IS A THRESHOLD FOR
! ZERO PIVOTAL ELEMENTS OF A.

! ERROR CODES:  IERROR = 0 - SOLUTION FOUND SUCCESSFULLY;
! IERROR = 1 - ZERO PIVOT ENCOUNTERED, SINGULAR MATRIX.


! This code was written by Steve Dixon for use in Divcon, and has been
! slightly adapted for use in this code. -Ed Brothers.

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION A(ISIZE,ISIZE),B(ISIZE),W(ISIZE),X(ISIZE)
    IERROR = 0
    IF(THRESH <= 0.0D0) THRESH = 1.0D-12
    IF(N == 1)THEN
        IF(ABS(A(1,1)) < THRESH)THEN
            IERROR = 1
        ELSE
            X(1) = B(1)/A(1,1)
        ENDIF
        RETURN
    ENDIF

! COMPUTE NORM OF A AS THE AVERAGE ABSOLUTE VALUE.

    AVE = 0.0D0
    DO 20 I=1,N
        DO 10 J=1,N
            AVE = AVE + ABS(A(I,J))
        10 ENDDO
    20 ENDDO
    AVE = AVE/(N*N)
    AMIN = AVE*THRESH

! BEGIN GAUSSIAN ELIMINATION.

    DO 100 K=1,N-1
    
    ! CHECK ENTRIES K THROUGH N OF THE KTH COLUMN OF A TO FIND THE
    ! LARGEST VALUE.
    
        AIKMAX = 0.0D0
        IROW = K
        DO 30 I=K,N
            AIK = ABS(A(I,K))
            IF(AIK > AIKMAX)THEN
                AIKMAX = AIK
                IROW = I
            ENDIF
        30 ENDDO
    
    ! AIKMAX IS THE ABSOLUTE VALUE OF THE PIVOTAL ELEMENT.  IF
    ! AIKMAX IS SMALLER THAN THE THRESHOLD AMIN THEN THE LEADING
    ! SUBMATRIX IS NEARLY SINGULAR.  IN THIS CASE, RETURN TO CALLING
    ! ROUTINE WITH AN ERROR.
    
        IF(AIKMAX < AMIN)THEN
            IERROR = 1
            RETURN
        ENDIF
    
    ! IF IROW IS NOT EQUAL TO K THEN SWAP ROWS K AND IROW OF THE
    ! MATRIX A AND SWAP ENTRIES K AND IROW OF THE VECTOR B.
    
        IF(IROW /= K)THEN
            DO 60 J=1,N
                W(J) = A(IROW,J)
                A(IROW,J) = A(K,J)
                A(K,J) = W(J)
            60 ENDDO
            BSWAP = B(IROW)
            B(IROW) = B(K)
            B(K) = BSWAP
        ELSE
            DO 70 J=K+1,N
                W(J) = A(K,J)
            70 ENDDO
        ENDIF
    
    ! FOR J GREATER THAN OR EQUAL TO I, OVERWRITE A(I,J) WITH
    ! U(I,J); FOR I GREATER THAN J OVERWRITE A(I,J) WITH L(I,J).
    
        DO 90 I=K+1,N
            T = A(I,K)/A(K,K)
            A(I,K) = T
            DO 80 J=K+1,N
                A(I,J) = A(I,J) - T*W(J)
            80 ENDDO
        90 ENDDO
    100 ENDDO
    IF(ABS(A(N,N)) < AMIN)THEN
        IERROR = 1
        RETURN
    ENDIF

! WE NOW HAVE STORED IN A THE L-U DECOMPOSITION OF P*A, WHERE
! P IS A PERMUTATION MATRIX OF THE N BY N IDENTITY MATRIX.  IN
! THE VECTOR B WE HAVE P*B.  WE NOW SOLVE L*U*X = P*B FOR THE
! VECTOR X.  FIRST OVERWRITE B WITH THE SOLUTION TO L*Y = B
! VIA FORWARD ELIMINATION.

    DO 160 I=2,N
        DO 150 J=1,I-1
            B(I) = B(I) - A(I,J)*B(J)
        150 ENDDO
    160 ENDDO

! NOW SOLVE U*X = B FOR X VIA BACK SUBSTITUTION.

    X(N) = B(N)/A(N,N)
    DO 200 K=2,N
        I = N+1-K
        X(I) = B(I)
        DO 190 J=I+1,N
            X(I) = X(I) - A(I,J)*X(J)
        190 ENDDO
        X(I) = X(I)/A(I,I)
    200 ENDDO
    RETURN
    end SUBROUTINE LSOLVE
! Ed Brothers. January 22, 2002
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    double precision function ssw(gridx,gridy,gridz,iparent)
    use allmod
    implicit double precision(a-h,o-z)

! This subroutie calculates the Scuseria-Stratmann wieghts.  There are
! two conditions that cause the weights to be unity: If there is only
! one atom:

    IF (natom == 1)  THEN
        ssw=1.d0
        return
    ENDIF

! Another time the weight is unity is r(iparent,g)<.5*(1-a)*R(i,n)
! where r(iparent,g) is the distance from the parent atom to the grid
! point, a is a parameter (=.64) and R(i,n) is the distance from the
! parent atom to it's nearest neighbor.

    xparent=xyz(1,iparent)
    yparent=xyz(2,iparent)
    zparent=xyz(3,iparent)

    rig=(gridx-xparent)*(gridx-xparent)
    rig=rig+(gridy-yparent)*(gridy-yparent)
    rig=rig+(gridz-zparent)*(gridz-zparent)
    rig=Dsqrt(rig)

    IF (rig < 0.18d0*distnbor(iparent)) THEN
        ssw=1.d0
        return
    ENDIF

! If neither of those are the case, we have to actually calculate the
! weight.  First we must calculate the unnormalized wieght of the grid point
! with respect to the parent atom.

! Step one of calculating the unnormalized weight is finding the confocal
! elliptical coordinate between each cell.  This it the mu with subscripted
! i and j in the paper:
! Stratmann, Scuseria, and Frisch, Chem. Phys. Lett., v 257,
! 1996, pg 213-223.

    wofparent=1.d0
    Jatm=1
    DO WHILE (Jatm.ne.iparent.and.wofparent.ne.0.d0)
        xJatm=xyz(1,Jatm)
        yJatm=xyz(2,Jatm)
        zJatm=xyz(3,Jatm)
        rjg=(gridx-xJatm)*(gridx-xJatm)
        rjg=rjg+(gridy-yJatm)*(gridy-yJatm)
        rjg=rjg+(gridz-zJatm)*(gridz-zJatm)
        rjg=Dsqrt(rjg)
        Rij=(xparent-xJatm)*(xparent-xJatm)
        Rij=Rij+(yparent-yJatm)*(yparent-yJatm)
        Rij=Rij+(zparent-zJatm)*(zparent-zJatm)
        Rij=Dsqrt(Rij)
        confocal=(rig-rjg)/Rij
        IF (confocal >= 0.64d0) THEN
        ! gofconfocal=1.d0
        ! wofparent=wofparent*.5d0*(1.d0-1.d0)
            wofparent=0.d0
        ELSEIF(confocal >= -0.64d0) THEN
            frctn=confocal/0.64d0
            frctnto3=frctn*frctn*frctn
            frctnto5=frctnto3*frctn*frctn
            frctnto7=frctnto5*frctn*frctn
            gofconfocal=(35.d0*frctn-35.d0*frctnto3+21.d0*frctnto5 &
            -5.d0*frctnto7)/16.d0
            wofparent=wofparent*.5d0*(1.d0-gofconfocal)
        ELSE
        ! gofconfocal=-1.d0
        ! wofparent=wofparent*.5d0*(1.d0-(-1.d0))
        ! wofparent=wofparent
            continue
        ENDIF
        Jatm=Jatm+1
    ENDDO

    Jatm=iparent+1
    DO WHILE (Jatm.le.natom.and.wofparent.ne.0.d0)
        xJatm=xyz(1,Jatm)
        yJatm=xyz(2,Jatm)
        zJatm=xyz(3,Jatm)
        rjg=(gridx-xJatm)*(gridx-xJatm)
        rjg=rjg+(gridy-yJatm)*(gridy-yJatm)
        rjg=rjg+(gridz-zJatm)*(gridz-zJatm)
        rjg=Dsqrt(rjg)
        Rij=(xparent-xJatm)*(xparent-xJatm)
        Rij=Rij+(yparent-yJatm)*(yparent-yJatm)
        Rij=Rij+(zparent-zJatm)*(zparent-zJatm)
        Rij=Dsqrt(Rij)
        confocal=(rig-rjg)/Rij
        IF (confocal >= 0.64d0) THEN
        ! gofconfocal=1.d0
        ! wofparent=wofparent*.5d0*(1.d0-1.d0)
            wofparent=0.d0
        ELSEIF(confocal >= -0.64d0) THEN
            frctn=confocal/0.64d0
            frctnto3=frctn*frctn*frctn
            frctnto5=frctnto3*frctn*frctn
            frctnto7=frctnto5*frctn*frctn
            gofconfocal=(35.d0*frctn-35.d0*frctnto3+21.d0*frctnto5 &
            -5.d0*frctnto7)/16.d0
            wofparent=wofparent*.5d0*(1.d0-gofconfocal)
        ELSE
        ! gofconfocal=-1.d0
        ! wofparent=wofparent*.5d0*(1.d0-(-1.d0))
        ! wofparent=wofparent
            continue
        ENDIF
        Jatm=Jatm+1
    ENDDO

    totalw=wofparent
    IF (wofparent == 0.d0) THEN
        ssw=0.d0
        return
    ENDIF

! Now we have the unnormalized weight of the grid point with regard to the
! parent atom.  Now we have to do this for all other atom pairs to
! normalize the grid weight.


    DO Iatm=1,natom
        IF (iatm == iparent) goto 50
        xIatm=xyz(1,Iatm)
        yIatm=xyz(2,Iatm)
        zIatm=xyz(3,Iatm)
        wofiatm=1.d0
        rig=(gridx-xIatm)*(gridx-xIatm)
        rig=rig+(gridy-yIatm)*(gridy-yIatm)
        rig=rig+(gridz-zIatm)*(gridz-zIatm)
        rig=Dsqrt(rig)
        Jatm=1
        DO WHILE (Jatm.ne.Iatm.and.wofiatm.ne.0.d0)
            rjg=(gridx-xyz(1,Jatm))*(gridx-xyz(1,Jatm))
            rjg=rjg+(gridy-xyz(2,Jatm))*(gridy-xyz(2,Jatm))
            rjg=rjg+(gridz-xyz(3,Jatm))*(gridz-xyz(3,Jatm))
            rjg=Dsqrt(rjg)
            Rij=(xIatm-xyz(1,Jatm))*(xIatm-xyz(1,Jatm))
            Rij=Rij+(yIatm-xyz(2,Jatm))*(yIatm-xyz(2,Jatm))
            Rij=Rij+(zIatm-xyz(3,Jatm))*(zIatm-xyz(3,Jatm))
            Rij=Dsqrt(Rij)
            confocal=(rig-rjg)/Rij
            IF (confocal >= 0.64d0) THEN
            ! gofconfocal=1.d0
            ! wofiatm=wofiatm*.5d0*(1.d0-1.d0)
                wofiatm=0.d0
            ELSEIF(confocal >= -0.64d0) THEN
                frctn=confocal/0.64d0
                frctnto3=frctn*frctn*frctn
                frctnto5=frctnto3*frctn*frctn
                frctnto7=frctnto5*frctn*frctn
                gofconfocal=(35.d0*frctn-35.d0*frctnto3+21.d0*frctnto5 &
                -5.d0*frctnto7)/16.d0
                wofiatm=wofiatm*.5d0*(1.d0-gofconfocal)
            ELSE
            ! gofconfocal=-1.d0
            ! wofiatm=wofiatm*.5d0*(1.d0-(-1.d0))
            ! wofiatm=wofiatm
                continue
            ENDIF
            Jatm=Jatm+1
        ENDDO

        Jatm=Iatm+1
        DO WHILE (Jatm.le.natom.and.wofiatm.ne.0.d0)
            rjg=(gridx-xyz(1,Jatm))*(gridx-xyz(1,Jatm))
            rjg=rjg+(gridy-xyz(2,Jatm))*(gridy-xyz(2,Jatm))
            rjg=rjg+(gridz-xyz(3,Jatm))*(gridz-xyz(3,Jatm))
            rjg=Dsqrt(rjg)
            Rij=(xIatm-xyz(1,Jatm))*(xIatm-xyz(1,Jatm))
            Rij=Rij+(yIatm-xyz(2,Jatm))*(yIatm-xyz(2,Jatm))
            Rij=Rij+(zIatm-xyz(3,Jatm))*(zIatm-xyz(3,Jatm))
            Rij=Dsqrt(Rij)
            confocal=(rig-rjg)/Rij
            IF (confocal >= 0.64d0) THEN
            ! gofconfocal=1.d0
            ! wofiatm=wofiatm*.5d0*(1.d0-1.d0)
                wofiatm=0.d0
            ELSEIF(confocal >= -0.64d0) THEN
                frctn=confocal/0.64d0
                frctnto3=frctn*frctn*frctn
                frctnto5=frctnto3*frctn*frctn
                frctnto7=frctnto5*frctn*frctn
                gofconfocal=(35.d0*frctn-35.d0*frctnto3+21.d0*frctnto5 &
                -5.d0*frctnto7)/16.d0
                wofiatm=wofiatm*.5d0*(1.d0-gofconfocal)
            ELSE
            ! gofconfocal=-1.d0
            ! wofiatm=wofiatm*.5d0*(1.d0-(-1.d0))
            ! wofiatm=wofiatm
                continue
            ENDIF
            Jatm=Jatm+1
        ENDDO

        totalw = totalw+wofiatm
        50 continue
    ENDDO

    ssw=wofparent/totalw

    END function ssw

