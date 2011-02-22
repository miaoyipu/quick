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
! RANdoM                      : Random
! PriLab                      : elegent way to output atom list
! EffChar                     : check valid char
! bndang                      : return band angel
! dihedr                      : return dihedral angle
! wrtrestart                  : rewrite input file
! copySym                     : complete symmetry of matrix
! sum2Mat                     : Sum over 2 matrix
!
!********************************************************


!--------------------------------------------------------
! Upcase(string,iend)
!--------------------------------------------------------
! CHANGES LOWER CASE CHARACTERS IN STRING(1:IEND) TO UPPER CASE.
! -- S. DIXON.
!--------------------------------------------------------
SUBROUTINE UPCASE(STRING,IEND)
  CHARACTER STRING*(*),AUPP,ALOW,DUMMY
  INTEGER :: IALOW,IAUPP,IDUMMY
  DATA AUPP,ALOW /'A','a'/
  SAVE AUPP,ALOW
  IALOW = ICHAR(ALOW)
  IAUPP = ICHAR(AUPP)
  ISHifT = IALOW - IAUPP
  do I=1,IEND
     DUMMY = STRING(I:I)
     IDUMMY = ICHAR(DUMMY)
     if(IDUMMY >= IALOW)then
        DUMMY = CHAR(IDUMMY - ISHifT)
        STRING(I:I) = DUMMY
     endif
  enddo
  RETURN
end SUBROUTINE UPCASE





!--------------------------------------------------------
! Upcase(string,iend)
!------------------------------------------------------------------------
! This was taken from the DivCon code by Steve Dixon.
! - Ed Brothers, February 18, 2002
!------------------------------------------------------------------------
SUBROUTINE RDWORD(STRING,ISTART,ISTOP)

  ! LOCATES THE NEXT WORD IN STRING STARTING AT STRING(ISTART:ISTART).
  ! ON RETURN ISTART AND ISTOP WILL BE UPDATED AND THE WORD WILL BE
  ! IN STRING(ISTART:ISTOP).  if THERE ARE NO MORE WORDS IN THE STRING,
  ! then BOTH ISTART AND ISTOP WILL BE RETURNED WITH VALUES OF ZERO.

  CHARACTER STRING*(*)
  LOGICAL :: INWORD

  ! GET DECLARED LENGTH OF STRING AND FIND THE NEXT CONTIGUOUS BLOCK
  ! OF NONBLANK CHARACTERS.

  IBEGIN = MAX(ISTART,1)
  IEND = LEN(STRING)
  if(IEND < IBEGIN)then
     ISTART = 0
     ISTOP = 0
     RETURN
  endif
  INWORD = .FALSE.
  IBEGIN = ISTART
  do 100 I=IBEGIN,IEND
     if(STRING(I:I) == ' ')then
        if(INWORD)then
           ISTOP = I-1
           RETURN
        endif
     else
        if( .NOT. INWORD)then
           INWORD = .TRUE.
           ISTART = I
        endif
     endif
100 enddo

  ! if WE GET HERE, then EITHER THE WORD FOUND EXTENDS ALL THE WAY TO
  ! THE END OF THE STRING, OR NO WORD WAS FOUND IN THE REMAINING
  ! PORTION OF THE STRING.

  if(INWORD)then
     ISTOP = IEND
  else
     ISTART = 0
     ISTOP = 0
  endif
  RETURN
end SUBROUTINE RDWORD








!--------------------------------------------------------
! rdnum
!--------------------------------------------------------
! contains all the routines to convert strings to numbers
!--------------------------------------------------------
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
10 enddo
20 if(iend < ibeg)then
     ierror = 1
     go to 1000
  endif
  ieq = index(string(ibeg:iend),'=')
  if(ieq /= 0) ibeg = ibeg + ieq
  call getnum(string,ibeg,iend,value,ierror)
1000 return
end subroutine rdnum



!--------------------------------------------------------
! getnum
!--------------------------------------------------------
! routine to read number from string
!--------------------------------------------------------
subroutine getnum(string,ibeg,iend,value,ierror)
  implicit double precision(a-h,o-z)

  character string*(*),char*1,efield(4)*1
  data efield /'E','e','D','d'/
  save efield

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
30 enddo
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




!------------------------------------------------------------
!  whole
!------------------------------------------------------------
! returns the whole number in the field string(ibeg:iend).  only
! the numbers 0-9 are allowed to be present.
!------------------------------------------------------------
subroutine whole(string,ibeg,iend,value,ierror)
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
10 enddo
1000 return
end subroutine whole





!-----------------------------------------------------------
! rdinum(string,istart,ivalue,ierror)
!-----------------------------------------------------------
! Read integer number
!-----------------------------------------------------------

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
10 enddo
20 if(iend < ibeg)then
     ierror = 1
     go to 1000
  endif
  ieq = index(string(ibeg:iend),'=')
  if(ieq /= 0) ibeg = ibeg + ieq
  call getinum(string,ibeg,iend,ivalue,ierror)
1000 return
end subroutine rdinum



!-----------------------------------------------------------
! getinum(string,istart,ivalue,ierror)
!-----------------------------------------------------------
! read integer number
!-----------------------------------------------------------
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




!-----------------------------------------------------------
! iwhole
!-----------------------------------------------------------
! returns the whole number in the field string(ibeg:iend).  only
! the numbers 0-9 are allowed to be present.
!-----------------------------------------------------------
subroutine iwhole(string,ibeg,iend,ivalue,ierror)

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
10 enddo
1000 return
end subroutine iwhole







!------------------------------------------------------------
! iatioi
!------------------------------------------------------------

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



!------------------------------------------------------------
! iatoimp
!------------------------------------------------------------

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

!-----------------------------------------------------------
! whatis1
!-----------------------------------------------------------

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
! whatis1i
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
! whatis2
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
! whatis7
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






!-----------------------------------------------------------
! xnorm(a,i,j,k)
!-----------------------------------------------------------
double precision function xnorm(a,i,j,k)
  use quick_constants_module
  implicit double precision(a-h,o-z)
  !-----------------------------------------------------------
  ! Ed Brothers. October 3, 2001
  ! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
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
  !-----------------------------------------------------------

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
  do L=1,i
     gamma = gamma * (dble(i - L) + .5d0)
  enddo
  do L=1,j
     gamma = gamma * (dble(j - L) + .5d0)
  enddo
  do L=1,k
     gamma = gamma * (dble(k - L) + .5d0)
  enddo

  xnorm = xnorm*pi**(-0.75d0)*gamma**(-.5d0)

  return
END function xnorm



!-----------------------------------------------------------
! xnewnorm(a,i,j,k)
!--------------------------------------------------------
! Ed Brothers. October 3, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
!-----------------------------------------------------------
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
  do Icon1 = 1,nprim1
     do Icon2 = 1,nprim1
        xnewnorm = xnewnorm + ddcoe(Icon1)*ddcoe(Icon2) &
             *(aaexp(Icon1)+aaexp(Icon2))**xponent
        !               print*,xnewnorm
     enddo
  enddo
  !        print*,'*',xnewnorm
  gamma=1.d0
  do L=1,i
     gamma = gamma * (dble(i - L) + .5d0)
  enddo
  do L=1,j
     gamma = gamma * (dble(j - L) + .5d0)
  enddo
  do L=1,k
     gamma = gamma * (dble(k - L) + .5d0)
  enddo
  xnewnorm = (xnewnorm*gamma*pi32)**(-.5d0)
  !        print*,'**',gamma,pi32,xnewnorm
  !        do Icon1 = 1,ncontract(Ibas)
  !            dcoeff(Icon1,Ibas) = dconew*dcoeff(Icon1,Ibas)
  !        enddo

  return
END function xnewnorm


!-----------------------------------------------------------
! vett
!------------------------------------------------------------
! Alessandro GENONI 03/12/007
! Subroutine to build up the array kvett, whose elemts are kvett(i)=i*(i-1)/2
!-----------------------------------------------------------
Subroutine vett
  use quick_ecp_module
  implicit double precision (a-h,o-z)
  do i=1,nbf12
     kvett(i)=i*(i-1)/2
  end do
  return
end Subroutine vett


!-----------------------------------------------------------
! PrtTime
!-----------------------------------------------------------
!  Print out total Time used by the job
!*     in parent program def. "real*8 Tim0,CPUTim"; "Tim0=CPUTim(0)" for initial time
!-----------------------------------------------------------
subroutine PrtTim(IOut,sec)
  Implicit Real*8(A-H,O-Z)

1000 Format(' Job cpu time:',I3,' days ',I2,' hours ',I2,' minutes ',F4.1,' seconds.')

  Time = sec
  NDays = (Time / (3600.0d0*24.0d0))
  Time = Time - (NDays*(3600.0d0*24.0d0))
  NHours = (Time / 3600.0d0)
  Time = Time - (NHours*3600.0d0)
  NMin = (Time / 60.0d0)
  Time = Time - (NMin*60.0d0)
  Write(IOut,1000) NDays, NHours, NMin, Time
  Return
End subroutine PrtTim





!-----------------------------------------------------------
! PrtAct
!-----------------------------------------------------------
! print action
!-----------------------------------------------------------

subroutine PrtAct(io,line)
  implicit none
  integer io,L,leng,i
  character line*(*)

  leng=len(line)
  L=0
  write(io,'(a)')
  write(io,'(" @ ",a)') line
  write(io,'(a)')
  call flush(io)
  return
end subroutine PrtAct




!-----------------------------------------------------------
! PriCol
!-----------------------------------------------------------
! 2004.12.22 print column(c1:c2) of mat(n,n) to file(io)
!-----------------------------------------------------------
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
101 format('(i7,5',a1,i2,'.',i2,')')
102 format('(i7,',i2,a1,i2,'.',i2,')')
  write(fmt3,103) x-7; write(fmt4,104) nf,x-7
103 format('(3x,5(',i2,'x,i7))')
104 format('(3x,',i2,'(',i2,'x,i7))')

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

end subroutine PriCol



!-----------------------------------------------------------      
! PriSym
!-----------------------------------------------------------
! print symmetric mat(n,n) to file(io)
!-----------------------------------------------------------
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
101 format('(i7,5',a1,i2,'.',i2,')')
102 format('(i7,',i2,a1,i2,'.',i2,')')
  write(fmt3,103) x-7
  write(fmt4,104) nf,x-7
103 format('(3x,5(',i2,'x,i7))')
104 format('(3x,',i2,'(',i2,'x,i7))')

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

end subroutine PriSym






!-----------------------------------------------------------
! Diag
!------------------------------------------------------------
! S. Dixon's diagonalization code from DivCon.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
!-----------------------------------------------------------

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
  IMPLICIT doUBLE PRECISION (A-H,O-Z)
  ! DIMENSION A(NDIM,*),V(3,*),EVAL1(*),IDEGEN1(*),EVEC1(NDIM,*)
  !    DIMENSION A(nbasis,nbasis),V(3,nbasis),EVAL1(nbasis), &
  !    IDEGEN1(nbasis),EVEC1(nbasis,nbasis)

  DIMENSION A(NDIM,NDIM),V(3,NDIM),IDEGEN1(NDIM),EVAL1(NDIM)
  DIMENSION EVEC1(NDIM,NDIM)

  ! FLAG FOR WHETHER OR NOT TO COMPUTE EIGENVALUES:


  if(NDIM == 1)then
     EVAL1(1) = A(1,1)
     IDEGEN1(1) = 1
     EVEC1(1,1) = 1.0D0
     RETURN
  endif

  ! TRIDIAGONALIZE THE MATRIX A.  THIS WILL OVERWRITE THE DIAGONAL
  ! AND SUBDIAGONAL OF A WITH THE TRIDIAGONALIZED VERSION.  THE
  ! HOUSEHOLDER VECTORS ARE RETURNED IN THE ROWS ABOVE THE DIAGONAL,
  ! AND THE BETAHS ARE RETURNED BELOW THE SUBDIAGONAL.

  CALL TRIDI(NDIM,V,A)

  ! COMPUTE NORM OF TRIDIAGONAL MATRIX FROM THE "LARGEST" COLUMN.

  ANORM = ABS(A(1,1)) + ABS(A(2,1))
  if(NDIM > 2)then
     do 20 I=2,NDIM-1
        AICOL = ABS(A(I-1,I)) + ABS(A(I,I)) + ABS(A(I+1,I))
        ANORM = MAX(ANORM,AICOL)
20   enddo
  endif
  ANCOL = ABS(A(NDIM-1,NDIM)) + ABS(A(NDIM,NDIM))
  ANORM = MAX(ANORM,ANCOL)

  ! GET EIGENVALUES AND DEGENERACIES OF THE TRIDIAGONAL MATRIX A.
  ! if THE CALLING ROUTINE HAS NOT SUPPLIED A TOLERANCE FACTOR FOR
  ! OFF-DIAGONAL ENTRIES IN THE QR ITERATION, A DEFAULT OF 1.0D-8
  ! WILL BE USED.

  TOLTMP = TOLERA
  if(TOLTMP <= 0.0D0) TOLTMP = 1.0D-8
  CALL EIGVAL(NDIM,A,V,TOLTMP,ANORM,EVAL1,IERROR)
  if(IERROR /= 0) RETURN

  ! DETERMINE DEGENERACIES OF EIGENVALUES.

  CALL DEGEN(NDIM,EVAL1,TOLTMP,ANORM,IDEGEN1)

  ! GET EIGENVECTORS OF TRIDIAGONALIZED VERSION OF A.

  if(NEVEC1 <= 0) RETURN
  CALL EIGVEC(NDIM,NEVEC1,A,V,TOLTMP,ANORM,EVAL1,IDEGEN1,EVEC1)

  ! PREMULTIPLY EVEC1 BY THE HOUSEHOLDER MATRIX USED TO TRIDIAGONALIZE
  ! A.  THIS TRANSFORMS EIGENVECTORS OF THE TRIDIAGONAL MATRIX TO
  ! THOSE OF THE ORIGINAL MATRIX A.  SEE SUBROUTINE TRIDI FOR
  ! STORAGE OF HOUSEHOLDER TRANSFORMATION.

  if(NDIM > 2)then
     do 30 K=1,NDIM-2
        V(1,K) = A(K+2,K)
30   enddo

     ! SWAP STORAGE SO THAT THE EXPENSIVE TRIPLE LOOP BELOW doESN'T
     ! HAVE TO JUMP ACROSS COLUMNS TO GET ENTRIES OF A.

     do 50 I=2,NDIM
        do 40 J=1,I-1
           A(I,J) = A(J,I)
40      enddo
50   enddo
     do 160 J=1,NEVEC1
        do 140 M=2,NDIM-1
           K = NDIM - M
           BETAH = V(1,K)
           !                if(ABS(BETAH) < 1.0D-50) GO TO 140
           if(ABS(BETAH) < 1.0D-50) CONTINUE
           SUM = 0.0D0
           do 100 I=K+1,NDIM
              SUM = SUM + A(I,K)*EVEC1(I,J)
100        enddo
           BSUM = BETAH*SUM
           do 120 I=K+1,NDIM
              EVEC1(I,J) = EVEC1(I,J) - A(I,K)*BSUM
120        enddo
140     enddo
160  enddo
  endif
  RETURN
end SUBROUTINE DIAG






!-----------------------------------------------------------
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
  IMPLICIT doUBLE PRECISION (A-H,O-Z)
  ! DIMENSION A(NDIM,*),V(3,*)
  DIMENSION A(nbasis,nbasis),V(3,nbasis)

  ! THRESH WILL BE USED AS A THRESHOLD TO DETERMINE if A VALUE SHOULD
  ! BE CONSIDERED TO BE ZERO.  THIS CAN BE CHANGED BY THE USER.

  THRESH = 1.0D-50

  ! if A IS 2 BY 2 OR SMALLER, then IT IS ALREADY TRIDIAGONAL -- NO
  ! NEED TO CONTINUE.

  if(NDIM <= 2) GO TO 1000
  do 500 K=1,NDIM-2

     ! DETERMINE THE VECTOR V USED IN THE HOUSEHOLDER TRANSFORMATION P.
     ! FOR EACH VALUE OF K THE HOUSEHOLDER MATRIX P IS DEFINED AS:

     ! P = I - BETAH*V*V'


     ! CONSTRUCT A HOUSEHOLDER TRANSFORMATION ONLY if THERE IS A NONZERO
     ! OFF-DIAGONAL ELEMENT BELOW A(K,K).

     ALPHA2 = 0.0D0
     do 60 I=K+1,NDIM
        V(1,I) = A(I,K)
        ALPHA2 = ALPHA2 + V(1,I)**2
60   enddo
     APTEMP = ALPHA2 - V(1,K+1)**2
     ALPHA = DSQRT(ALPHA2)
     if(ALPHA >= THRESH)then
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

        doT12 = 0.0D0
        BHALF = BETAH*0.5D0
        do 220 I=K+1,NDIM
           SUM = 0.0D0
           do 100 J=K+1,I
              SUM = SUM + A(I,J)*V(1,J)
100        enddo
           if(I < NDIM)then
              do 180 J=I+1,NDIM

                 ! AN UPPER TRIANGULAR ENTRY OF A WILL BE REQUIRED.  MUST USE
                 ! THE SYMMETRIC ENTRY IN THE LOWER TRIANGULAR PART OF A.

                 SUM = SUM + A(J,I)*V(1,J)
180           enddo
           endif
           V(2,I) = BETAH*SUM
           doT12 = doT12 + V(1,I)*V(2,I)
220     enddo
        BH12 = BHALF*doT12
        do 300 I=K+1,NDIM
           V(2,I) = V(2,I) - BH12*V(1,I)
300     enddo
        do 350 J=K+1,NDIM
           do 310 I=J,NDIM
              A(I,J) = A(I,J) - V(1,I)*V(2,J) - V(2,I)*V(1,J)
310        enddo

           ! STORE V(1,J) ABOVE THE DIAGONAL IN ROW K OF A

           A(K,J) = V(1,J)
350     enddo

        ! STORE BETAH BELOW THE SUBDIAGONAL OF A.

        A(K+2,K) = BETAH
     else

        ! NO HOUSEHOLDER TRANSFORMATION IS NECESSARY BECAUSE THE OFF-
        ! DIAGONALS ARE ALL ESSENTIALLY ZERO.

        A(K+2,K) = 0.0D0
        do 460 J=K+1,NDIM
           A(K,J) = 0.0D0
460     enddo
     endif
500 enddo
1000 RETURN
end SUBROUTINE TRIDI






!-----------------------------------------------------------
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
  ! IERROR = 1 if QR ITERATION DID NOT CONVERGE; 0 OTHERWISE.

  ! PROGRAMMED BY S. L. DIXON.


  use allmod
  IMPLICIT doUBLE PRECISION (A-H,O-Z)
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
  if(NDIM > 2)then
     do 50 I=2,NDIM-1
        EVAL1(I) = A(I,I)
        BETAH(1,I) = A(I+1,I)
50   enddo
  endif
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
  do 100 K=1,KUPPER
     I = N-K
     if(ABS(BETAH(1,I)) <= EPSLON)then
        L = I+1
        GO TO 150
     endif
100 enddo

  ! if WE GET TO THE NEXT STATEMENT, then THERE ARE NO ZERO OFF-DIAGONALS
  ! FOR THE SUBMATRIX WITH UPPER LEFT A(L,L) AND LOWER RIGHT A(N,N).
  ! WE CAN STILL GET EIGENVALUES if THE MATRIX IS 2 BY 2 OR 1 BY 1.
  ! OTHERWISE, do ANOTHER QR ITERATION PROVIDED ITMAX CYCLES HAVE
  ! NOT OCCURRED.

  if(L == N .OR. L == N-1)then
     GO TO 150
  else
     if(ITER == ITMAX)then
        IERROR = 1
        GO TO 1000
     else
        GO TO 200
     endif
  endif

  ! if WE GET TO 150 then A(L,L-1) IS ZERO AND THE UNREDUCED SUBMATRIX
  ! HAS UPPER LEFT AT A(L,L) AND LOWER RIGHT AT A(N,N).  WE CAN
  ! EXTRACT ONE EIGENVALUE if THIS MATRIX IS 1 BY 1 AND 2 EIGENVALUES
  ! if IT IS 2 BY 2.

150 if(L == N)then

     ! IT'S A 1 BY 1 AND EVAL1(N) IS AN EIGENVALUE.  if L=2 OR 1 WE ARE
     ! doNE.  OTHERWISE, UPDATE N, RESET L AND ITER, AND REPEAT THE
     ! SEARCH.

     if(L <= 2)then
        GO TO 500
     else
        N = L-1
        L = 1
        ITER = 0
        GO TO 80
     endif
  elseif(L == N-1)then

     ! THE UNREDUCED SUBMATRIX IS A 2 BY 2.  OVERWRITE EVAL1(N-1)
     ! AND EVAL1(N) WITH THE EIGENVALUES OF THE LOWER RIGHT 2 BY 2.

     BTERM = EVAL1(N-1) + EVAL1(N)
     ROOT1 = BTERM*0.5D0
     ROOT2 = ROOT1
     DISCR = BTERM**2 - 4.0D0*(EVAL1(N-1)*EVAL1(N)-BETAH(1,N-1)**2)
     if(DISCR > 0.0D0)then
        D = DSQRT(DISCR)*0.5D0
        ROOT1 = ROOT1 - D
        ROOT2 = ROOT2 + D
     endif
     EVAL1(N-1) = ROOT1
     EVAL1(N) = ROOT2

     ! SEE if WE ARE doNE.  if NOT, RESET N, L, AND ITER AND LOOK
     ! FOR NEXT UNREDUCED SUBMATRIX.

     if(L <= 2)then
        GO TO 500
     else
        N = L-1
        L = 1
        ITER = 0
        GO TO 80
     endif
  else

     ! AN EIGENVALUE WAS FOUND AND THE NEW UNREDUCED MATRIX LIMITS
     ! N AND L ARE SET.  do A QR ITERATION ON NEW MATRIX.

     ITER = 0
     GO TO 200
  endif

  ! QR ITERATION BEGINS HERE.

200 ITER = ITER + 1

  ! USE EIGENVALUES OF THE LOWER RIGHT 2 BY 2 TO COMPUTE SHifT.  SHifT
  ! BY THE EIGENVALUE CLOSEST TO EVAL1(N).

  D = (EVAL1(N-1) - EVAL1(N))*0.5D0
  SIGND = 1.0D0
  if(D < 0.0D0) SIGND = -1.0D0
  SHifT = EVAL1(N) + D - SIGND*DSQRT(D*D + BETAH(1,N-1)**2)
  P = EVAL1(L) - SHifT
  R = BETAH(1,L)
  T = EVAL1(L)
  W = BETAH(1,L)

  ! OVERWRITE A WITH Q'*A*Q.

  do 250 K=L,N-1
     D = DSQRT(P*P + R*R)
     C = P/D
     S = R/D
     if(K /= L) BETAH(1,K-1) = D
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
250 enddo
  BETAH(1,N-1) = P
  EVAL1(N) = T

  ! GO BACK AND SEE if L AND N NEED TO BE UPDATED.

  GO TO 80

  ! SORT EIGENVALUES IN ASCENDING ALGEBRAIC ORDER.

500 do 600 I=2,NDIM
     JMAX = NDIM-I+1
     ISORT = 0
     do 550 J=1,JMAX
        if(EVAL1(J) > EVAL1(J+1))then
           ETEMP = EVAL1(J)
           EVAL1(J) = EVAL1(J+1)
           EVAL1(J+1) = ETEMP
           ISORT = 1
        endif
550  enddo
     if(ISORT == 0) GO TO 1000
600 enddo
1000 RETURN
end SUBROUTINE EIGVAL



!-----------------------------------------------------------
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
  IMPLICIT doUBLE PRECISION (A-H,O-Z)
  ! DIMENSION EVAL1(*),IDEGEN1(*)
  DIMENSION EVAL1(nbasis),IDEGEN1(nbasis)

  ! DETERMINE DEGENERACIES OF EIGENVALUES.  ADJACENT EIGENVALUES
  ! WILL BE CONSIDERED TO BE DEGENERATE WHEN THEY DifFER BY LESS
  ! THAN DTOLER.

  DTOLER = MAX(ANORM*DSQRT(TOLERA),1.0D-8)
  NSAME = 1
  do 200 I=2,NDIM
     DifF = ABS(EVAL1(I-1) - EVAL1(I))
     if(DifF <= DTOLER)then

        ! EIGENVALUES I-1 AND I ARE DEGENERATE.

        NSAME = NSAME + 1
        if(I == NDIM)then

           ! WE'VE COME TO THE LAST REQUESTED EIGENVALUE, AND IT'S TIME
           ! TO ASSIGN DEGENERACIES FOR THE BLOCK ENDING WITH THE ITH
           ! EIGENVALUE.

           do 100 J=I-NSAME+1,I
              IDEGEN1(J) = NSAME
100        enddo
        endif

        ! GO TO THE NEXT EIGENVALUE (if THERE ARE ANY LEFT) AND SEE if
        ! IT'S DEGENERATE WITH THE NSAME EIGENVALUES WE'VE ALREADY
        ! FOUND.

        GO TO 200
     else

        ! EITHER EIGENVALUE I-1 IS NONDEGENERATE OR IT'S THE LAST
        ! EIGENVALUE IN A DEGENERATE BLOCK.  CORRESPONDINGLY, ASSIGN THE
        ! PROPER DEGENERACY TO I-1 OR TO EACH EIGENVALUE IN THE BLOCK.

        do 150 J=I-NSAME,I-1
           IDEGEN1(J) = NSAME
150     enddo
        NSAME = 1

        ! if I=NDIM then IT MUST BE THE CASE THAT THIS LAST EIGENVALUE
        ! IS NONDEGENERATE.

        if(I == NDIM) IDEGEN1(I) = 1
     endif
200 enddo
  RETURN
end SUBROUTINE DEGEN









!-----------------------------------------------------------
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
  IMPLICIT doUBLE PRECISION (A-H,O-Z)
  ! DIMENSION A(NDIM,*),AWORK(3,*),EVAL1(*),IDEGEN1(*),EVEC1(NDIM,*)
  DIMENSION A(nbasis,nbasis),AWORK(3,nbasis),EVAL1(nbasis), &
       IDEGEN1(nbasis),EVEC1(nbasis,nbasis)
  LOGICAL :: ORTH
  IRAND = 13876532

  ! COMPUTE THRESHOLD EPSLON WHICH WILL BE USED if THE INVERSE ITERATION
  ! MATRIX IS SINGULAR.

  EPSLON = ANORM*TOLERA

  ! WHEN DEGENERACIES OCCUR, THERE ARE RARE INSTANCES WHEN THE
  ! DEGENERATE BLOCK OF EIGENVECTORS ARE NOT LINEARLY INDEPENDENT.
  ! IN THESE CASES, AN ADDITIONAL PASS THROUGH THE INVERSE ITERATION
  ! (WITH A NEW SET OF RANdoM NUMBERS) IS CARRIED OUT, I.E., CONTROL
  ! PASSES TO STATEMENT 40.  IVECT WILL KEEP TRACK OF THE CURRENT
  ! STARTING EIGENVECTOR WHEN ADDITIONAL PASSES ARE NECESSARY.

  IVECT = 1
  NFAIL = 0
  NPRTRB = 0

  ! do ONE ITERATION FOR EACH EIGENVECTOR.

40 ORTH = .TRUE.
  NDEGEN = 0
  MSTART = IVECT
  do 380 M=MSTART,NEVEC1
     if(IDEGEN1(M) > 1) NDEGEN = NDEGEN + 1
     Z = EVAL1(M)

     ! if THE INVERSE ITERATION HAS FAILED TWICE DUE TO NON-ORTHOGONALITY
     ! OF DEGENERATE EIGENVECTORS, PERTURB THE EIGENVALUE BY A SMALL
     ! AMOUNT.

     if(NFAIL >= 2)then
        NPRTRB = NPRTRB + 1
        Z = Z + 0.001D0*dble(NPRTRB)*TOLERA
     endif

     ! STORE THE TRIDIAGONAL ENTRIES OF THE INVERSE ITERATION MATRIX IN
     ! THE 3 BY NDIM WORKSPACE AWORK.

     AWORK(1,1) = 0.0D0
     AWORK(2,1) = A(1,1) - Z
     AWORK(3,1) = A(2,1)
     if(NDIM > 2)then
        do 80 I=2,NDIM-1
           AWORK(1,I) = A(I,I-1)
           AWORK(2,I) = A(I,I) - Z
           AWORK(3,I) = A(I+1,I)
80      enddo
     endif
     AWORK(1,NDIM) = A(NDIM,NDIM-1)
     AWORK(2,NDIM) = A(NDIM,NDIM) - Z
     AWORK(3,NDIM) = 0.0D0

     ! ASSIGN INVERSE ITERATION VECTOR FROM RANdoM NUMBERS.

     do 120 I=1,NDIM
        CALL RANdoM(IRAND,RNdoM)
        RNdoM = 2.0D0*(RNdoM - 0.5D0)
        EVEC1(I,M) = RNdoM*TOLERA
120  enddo

     ! CARRY OUT FORWARD GAUSSIAN ELIMINATION WITH ROW PIVOTING
     ! ON THE INVERSE ITERATION MATRIX.

     do 160 K=1,NDIM-1
        ADIAG = ABS(AWORK(2,K))
        ASUB = ABS(AWORK(1,K+1))
        if(ADIAG >= ASUB)then

           ! USE PIVOTAL ELEMENT FROM ROW K.

           if(AWORK(2,K) == 0.0D0) AWORK(2,K) = EPSLON
           T = AWORK(1,K+1)/AWORK(2,K)
           AWORK(1,K+1) = 0.0D0
           AWORK(2,K+1) = AWORK(2,K+1) - T*AWORK(3,K)

           ! LEFT-JUSTifY EQUATION K SO THAT DIAGONAL ENTRY IS STORED
           ! IN AWORK(1,K).

           AWORK(1,K) = AWORK(2,K)
           AWORK(2,K) = AWORK(3,K)
           AWORK(3,K) = 0.0D0

           ! OPERATE ON VECTOR AS WELL.

           EVEC1(K+1,M) = EVEC1(K+1,M) - T*EVEC1(K,M)
        else

           ! USE PIVOTAL ELEMENT FROM ROW K+1 AND SWAP ROWS K AND K+1.

           if(AWORK(1,K+1) == 0.0D0) AWORK(1,K+1) = EPSLON
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
        endif
160  enddo

     ! FORWARD ELIMINATION COMPLETE.  BACK SUBSTITUTE TO GET SOLUTION.
     ! OVERWRITE COLUMN M OF EVEC1 WITH SOLUTION.

     if(AWORK(2,NDIM) == 0.0D0) AWORK(2,NDIM) = EPSLON
     EVEC1(NDIM,M) = EVEC1(NDIM,M)/AWORK(2,NDIM)
     ETEMP = EVEC1(NDIM-1,M) - AWORK(2,NDIM-1)*EVEC1(NDIM,M)
     EVEC1(NDIM-1,M) = ETEMP/AWORK(1,NDIM-1)
     ENORM = EVEC1(NDIM,M)**2 + EVEC1(NDIM-1,M)**2
     if(NDIM > 2)then

        ! CAUTION: PROBLEM LOOP FOR SOME IBM RS/6000 COMPILERS.  VALUE
        ! OF K CAN GET LOST WHEN OPTIMIZE FLAG IS USED.

        do 200 L=1,NDIM-2
           K = NDIM-L-1
           ETEMP = EVEC1(K,M) - AWORK(2,K)*EVEC1(K+1,M) &
                - AWORK(3,K)*EVEC1(K+2,M)
           EVEC1(K,M) = ETEMP/AWORK(1,K)
           ENORM = ENORM + EVEC1(K,M)**2
200     enddo
     endif
     EINV = 1.0D0/DSQRT(ENORM)

     ! NORMALIZE EIGENVECTOR.

     do 240 I=1,NDIM
        EVEC1(I,M) = EVEC1(I,M)*EINV
240  enddo

     ! if WE HAVE COME TO THE END OF A DEGENERATE BLOCK OF EIGENVECTORS,
     ! ORTHOGONALIZE THE BLOCK.

     if(NDEGEN > 1)then
        if(NDEGEN == IDEGEN1(M) .OR. M == NEVEC1)then
           JSTART = M-NDEGEN+1
           CALL ORTHOG(NDIM,NDEGEN,JSTART,EVEC1,ORTH)
           if(ORTH)then
              NFAIL = 0
              NPRTRB = 0

              ! THE DEGENERATE VECTORS WERE LINEARLY INDEPENDENT AND WERE
              ! SUCCESSFULLY ORTHOGONALIZED.

              IVECT = IVECT + NDEGEN
              NDEGEN = 0
           else

              ! THE BLOCK IS APPARENTLY NOT LINEARLY INDEPENDENT.  GO BACK
              ! AND REPEAT THE INVERSE ITERATION FOR THESE VECTORS.  AFTER
              ! AN INDEPENDENT SET HAS BEEN FOUND, ANY ADDITIONAL EIGENVECTORS
              ! WILL BE DETERMINED.

              NFAIL = NFAIL + 1
              GO TO 40
           endif
        endif
     endif

     ! THE CURRENT EIGENVECTOR SHOULD BE OKAY if IT IS NONDEGENERATE.

     if(IDEGEN1(M) == 1) IVECT = IVECT + 1
380 enddo
  RETURN
end SUBROUTINE EIGVEC




!-----------------------------------------------------------
! ORTHOG
!------------------------------------------------------------
SUBROUTINE ORTHOG(NDIM,NVECT,JSTART,VECT,ORTH)

  ! CONSTRUCTS A SET OF ORTHONORMAL VECTORS FROM THE NVECT LINEARLY
  ! INDEPENDENT, NORMALIZED VECTORS IN THE ARRAY VECT.  THE VECTORS
  ! SHOULD BE STORED COLUMNWISE, STARTING IN COLUMN JSTART.  VECT IS
  ! OVERWRITTEN WITH THE ORTHONORMAL SET.  ALL VECTORS ARE NDIM BY 1.
  ! ORTH IS RETURNED WITH A VALUE OF .TRUE. if THE SET WAS LINEARLY
  ! INDEPENDENT AND .FALSE. OTHERWISE.

  ! PROGRAMMED BY S. L. DIXON.


  use allmod
  IMPLICIT doUBLE PRECISION (A-H,O-Z)
  ! DIMENSION VECT(NDIM,*)
  DIMENSION VECT(nbasis,nbasis)
  LOGICAL :: ORTH

  ORTH = .TRUE.
  ORTEST = 1.0D-8

  ! BEGIN ORTHOGONALIZATION.

  JSTOP = JSTART + NVECT - 1
  do 120 J=JSTART,JSTOP
     if(J > JSTART)then

        ! SUBTRACT OFF COMPONENTS OF PREVIOUSLY DETERMINED ORTHOGONAL
        ! VECTORS FROM THE VECTOR IN COLUMN J.

        do 60 JPREV=JSTART,J-1
           doT = 0.0D0
           do 20 I=1,NDIM
              doT = doT + VECT(I,JPREV)*VECT(I,J)
20         enddo
           do 40 I=1,NDIM
              VECT(I,J) = VECT(I,J) - doT*VECT(I,JPREV)
40         enddo
60      enddo
     endif

     ! NORMALIZE COLUMN J.

     VJNORM = 0.0D0
     do 80 I=1,NDIM
        VJNORM = VJNORM + VECT(I,J)**2
80   enddo
     VJNORM = DSQRT(VJNORM)

     ! if THE NORM OF THIS VECTOR IS TOO SMALL then THE VECTORS ARE
     ! NOT LINEARLY INDEPENDENT.

     if(VJNORM < ORTEST)then
        ORTH = .FALSE.
        GO TO 1000
     endif
     do 100 I=1,NDIM
        VECT(I,J) = VECT(I,J)/VJNORM
100  enddo
120 enddo
1000 RETURN
end SUBROUTINE ORTHOG






!-----------------------------------------------------------
! RANDOM
!------------------------------------------------------------

SUBROUTINE RANDOM(IX,X)

  ! GENERATES INTEGER*4 AND doUBLE PRECISION RANdoM NUMBERS ON
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

  ! then,

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

  ! ADD I31_1 BACK IN if THE SUBTRACTION MADE IX NEGATIVE.

  if(IX < 0) IX = IX + I31_1

  ! MAKE X RANdoM ON (0.0,1.0) BY MULTIPLYING IX BY 1.0/I31_1

  X = dble(IX)*4.6566128752458D-10
  RETURN
end SUBROUTINE RANdoM







!-----------------------------------------------------------
! PrtLab
!-----------------------------------------------------------
!Li Wei write Ktmp(kk) into line: eg. "1,3,5-8,2*0,9"
!-----------------------------------------------------------
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
110 enddo

end subroutine PrtLab





!-----------------------------------------------------------
! EffChar
!-----------------------------------------------------------
! 2005.01.07 move blank of two sides in a line
!-----------------------------------------------------------
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

end subroutine EffChar




!-----------------------------------------------------------
! ekinetic
!-----------------------------------------------------------
! Ed Brothers. October 12, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
!-----------------------------------------------------------
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
  if (kinetic == 0.d0) goto 100
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



!-----------------------------------------------------------
! LSOLVE
!-----------------------------------------------------------
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

  IMPLICIT doUBLE PRECISION (A-H,O-Z)
  DIMENSION A(ISIZE,ISIZE),B(ISIZE),W(ISIZE),X(ISIZE)
  IERROR = 0
  if(THRESH <= 0.0D0) THRESH = 1.0D-12
  if(N == 1)then
     if(ABS(A(1,1)) < THRESH)then
        IERROR = 1
     else
        X(1) = B(1)/A(1,1)
     endif
     RETURN
  endif

  ! COMPUTE NORM OF A AS THE AVERAGE ABSOLUTE VALUE.

  AVE = 0.0D0
  do 20 I=1,N
     do 10 J=1,N
        AVE = AVE + ABS(A(I,J))
10   enddo
20 enddo
  AVE = AVE/(N*N)
  AMIN = AVE*THRESH

  ! BEGIN GAUSSIAN ELIMINATION.

  do 100 K=1,N-1

     ! CHECK ENTRIES K THROUGH N OF THE KTH COLUMN OF A TO FIND THE
     ! LARGEST VALUE.

     AIKMAX = 0.0D0
     IROW = K
     do 30 I=K,N
        AIK = ABS(A(I,K))
        if(AIK > AIKMAX)then
           AIKMAX = AIK
           IROW = I
        endif
30   enddo

     ! AIKMAX IS THE ABSOLUTE VALUE OF THE PIVOTAL ELEMENT.  if
     ! AIKMAX IS SMALLER THAN THE THRESHOLD AMIN then THE LEADING
     ! SUBMATRIX IS NEARLY SINGULAR.  IN THIS CASE, RETURN TO CALLING
     ! ROUTINE WITH AN ERROR.

     if(AIKMAX < AMIN)then
        IERROR = 1
        RETURN
     endif

     ! if IROW IS NOT EQUAL TO K then SWAP ROWS K AND IROW OF THE
     ! MATRIX A AND SWAP ENTRIES K AND IROW OF THE VECTOR B.

     if(IROW /= K)then
        do 60 J=1,N
           W(J) = A(IROW,J)
           A(IROW,J) = A(K,J)
           A(K,J) = W(J)
60      enddo
        BSWAP = B(IROW)
        B(IROW) = B(K)
        B(K) = BSWAP
     else
        do 70 J=K+1,N
           W(J) = A(K,J)
70      enddo
     endif

     ! FOR J GREATER THAN OR EQUAL TO I, OVERWRITE A(I,J) WITH
     ! U(I,J); FOR I GREATER THAN J OVERWRITE A(I,J) WITH L(I,J).

     do 90 I=K+1,N
        T = A(I,K)/A(K,K)
        A(I,K) = T
        do 80 J=K+1,N
           A(I,J) = A(I,J) - T*W(J)
80      enddo
90   enddo
100 enddo
  if(ABS(A(N,N)) < AMIN)then
     IERROR = 1
     RETURN
  endif

  ! WE NOW HAVE STORED IN A THE L-U DECOMPOSITION OF P*A, WHERE
  ! P IS A PERMUTATION MATRIX OF THE N BY N IDENTITY MATRIX.  IN
  ! THE VECTOR B WE HAVE P*B.  WE NOW SOLVE L*U*X = P*B FOR THE
  ! VECTOR X.  FIRST OVERWRITE B WITH THE SOLUTION TO L*Y = B
  ! VIA FORWARD ELIMINATION.

  do 160 I=2,N
     do 150 J=1,I-1
        B(I) = B(I) - A(I,J)*B(J)
150  enddo
160 enddo

  ! NOW SOLVE U*X = B FOR X VIA BACK SUBSTITUTION.

  X(N) = B(N)/A(N,N)
  do 200 K=2,N
     I = N+1-K
     X(I) = B(I)
     do 190 J=I+1,N
        X(I) = X(I) - A(I,J)*X(J)
190  enddo
     X(I) = X(I)/A(I,I)
200 enddo
  RETURN
end SUBROUTINE LSOLVE








!-----------------------------------------------------------
! ssw
!-----------------------------------------------------------
! Ed Brothers. January 22, 2002
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
!-----------------------------------------------------------

double precision function ssw(gridx,gridy,gridz,iparent)
  use allmod
  implicit double precision(a-h,o-z)

  ! This subroutie calculates the Scuseria-Stratmann wieghts.  There are
  ! two conditions that cause the weights to be unity: If there is only
  ! one atom:

  if (natom == 1)  then
     ssw=1.d0
     return
  endif

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

  if (rig < 0.18d0*distnbor(iparent)) then
     ssw=1.d0
     return
  endif

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
  do while (Jatm.ne.iparent.and.wofparent.ne.0.d0)
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
     if (confocal >= 0.64d0) then
        ! gofconfocal=1.d0
        ! wofparent=wofparent*.5d0*(1.d0-1.d0)
        wofparent=0.d0
     elseif(confocal >= -0.64d0) then
        frctn=confocal/0.64d0
        frctnto3=frctn*frctn*frctn
        frctnto5=frctnto3*frctn*frctn
        frctnto7=frctnto5*frctn*frctn
        gofconfocal=(35.d0*frctn-35.d0*frctnto3+21.d0*frctnto5 &
             -5.d0*frctnto7)/16.d0
        wofparent=wofparent*.5d0*(1.d0-gofconfocal)
     else
        ! gofconfocal=-1.d0
        ! wofparent=wofparent*.5d0*(1.d0-(-1.d0))
        ! wofparent=wofparent
        continue
     endif
     Jatm=Jatm+1
  enddo

  Jatm=iparent+1
  do while (Jatm.le.natom.and.wofparent.ne.0.d0)
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
     if (confocal >= 0.64d0) then
        ! gofconfocal=1.d0
        ! wofparent=wofparent*.5d0*(1.d0-1.d0)
        wofparent=0.d0
     elseif(confocal >= -0.64d0) then
        frctn=confocal/0.64d0
        frctnto3=frctn*frctn*frctn
        frctnto5=frctnto3*frctn*frctn
        frctnto7=frctnto5*frctn*frctn
        gofconfocal=(35.d0*frctn-35.d0*frctnto3+21.d0*frctnto5 &
             -5.d0*frctnto7)/16.d0
        wofparent=wofparent*.5d0*(1.d0-gofconfocal)
     else
        ! gofconfocal=-1.d0
        ! wofparent=wofparent*.5d0*(1.d0-(-1.d0))
        ! wofparent=wofparent
        continue
     endif
     Jatm=Jatm+1
  enddo

  totalw=wofparent
  if (wofparent == 0.d0) then
     ssw=0.d0
     return
  endif

  ! Now we have the unnormalized weight of the grid point with regard to the
  ! parent atom.  Now we have to do this for all other atom pairs to
  ! normalize the grid weight.


  do Iatm=1,natom
     if (iatm == iparent) goto 50
     xIatm=xyz(1,Iatm)
     yIatm=xyz(2,Iatm)
     zIatm=xyz(3,Iatm)
     wofiatm=1.d0
     rig=(gridx-xIatm)*(gridx-xIatm)
     rig=rig+(gridy-yIatm)*(gridy-yIatm)
     rig=rig+(gridz-zIatm)*(gridz-zIatm)
     rig=Dsqrt(rig)
     Jatm=1
     do while (Jatm.ne.Iatm.and.wofiatm.ne.0.d0)
        rjg=(gridx-xyz(1,Jatm))*(gridx-xyz(1,Jatm))
        rjg=rjg+(gridy-xyz(2,Jatm))*(gridy-xyz(2,Jatm))
        rjg=rjg+(gridz-xyz(3,Jatm))*(gridz-xyz(3,Jatm))
        rjg=Dsqrt(rjg)
        Rij=(xIatm-xyz(1,Jatm))*(xIatm-xyz(1,Jatm))
        Rij=Rij+(yIatm-xyz(2,Jatm))*(yIatm-xyz(2,Jatm))
        Rij=Rij+(zIatm-xyz(3,Jatm))*(zIatm-xyz(3,Jatm))
        Rij=Dsqrt(Rij)
        confocal=(rig-rjg)/Rij
        if (confocal >= 0.64d0) then
           ! gofconfocal=1.d0
           ! wofiatm=wofiatm*.5d0*(1.d0-1.d0)
           wofiatm=0.d0
        elseif(confocal >= -0.64d0) then
           frctn=confocal/0.64d0
           frctnto3=frctn*frctn*frctn
           frctnto5=frctnto3*frctn*frctn
           frctnto7=frctnto5*frctn*frctn
           gofconfocal=(35.d0*frctn-35.d0*frctnto3+21.d0*frctnto5 &
                -5.d0*frctnto7)/16.d0
           wofiatm=wofiatm*.5d0*(1.d0-gofconfocal)
        else
           ! gofconfocal=-1.d0
           ! wofiatm=wofiatm*.5d0*(1.d0-(-1.d0))
           ! wofiatm=wofiatm
           continue
        endif
        Jatm=Jatm+1
     enddo

     Jatm=Iatm+1
     do while (Jatm.le.natom.and.wofiatm.ne.0.d0)
        rjg=(gridx-xyz(1,Jatm))*(gridx-xyz(1,Jatm))
        rjg=rjg+(gridy-xyz(2,Jatm))*(gridy-xyz(2,Jatm))
        rjg=rjg+(gridz-xyz(3,Jatm))*(gridz-xyz(3,Jatm))
        rjg=Dsqrt(rjg)
        Rij=(xIatm-xyz(1,Jatm))*(xIatm-xyz(1,Jatm))
        Rij=Rij+(yIatm-xyz(2,Jatm))*(yIatm-xyz(2,Jatm))
        Rij=Rij+(zIatm-xyz(3,Jatm))*(zIatm-xyz(3,Jatm))
        Rij=Dsqrt(Rij)
        confocal=(rig-rjg)/Rij
        if (confocal >= 0.64d0) then
           ! gofconfocal=1.d0
           ! wofiatm=wofiatm*.5d0*(1.d0-1.d0)
           wofiatm=0.d0
        elseif(confocal >= -0.64d0) then
           frctn=confocal/0.64d0
           frctnto3=frctn*frctn*frctn
           frctnto5=frctnto3*frctn*frctn
           frctnto7=frctnto5*frctn*frctn
           gofconfocal=(35.d0*frctn-35.d0*frctnto3+21.d0*frctnto5 &
                -5.d0*frctnto7)/16.d0
           wofiatm=wofiatm*.5d0*(1.d0-gofconfocal)
        else
           ! gofconfocal=-1.d0
           ! wofiatm=wofiatm*.5d0*(1.d0-(-1.d0))
           ! wofiatm=wofiatm
           continue
        endif
        Jatm=Jatm+1
     enddo

     totalw = totalw+wofiatm
50   continue
  enddo

  ssw=wofparent/totalw

END function ssw


!-----------------------------------------------------------
! IOrder
!-----------------------------------------------------------
! 2004.12.23 find order for the first n1 a dimension arr(n)
!     incr=1: increasing order; incr=-1: decreasing order
!-----------------------------------------------------------

subroutine IOrder(n,arr)
  implicit none
  integer n,incr,i,j,l,k1,k2
  integer arr(n),PP,k

  do i=1,n-1
     do j=i+1,n
        if (arr(i).ge.arr(j)) then
           k=arr(i);arr(i)=arr(j);arr(j)=k
        endif
     enddo
  enddo

end subroutine IOrder


!-----------------------------------------------------------
! PrtDate
!-----------------------------------------------------------
! 2004.12.24 Print current time for system
!-----------------------------------------------------------

subroutine PrtDate(io,note)
  implicit none
  integer io,i
  character datim*26,note*(*)

  i=len(note)
  call GDate(datim)
  write (io,*) note(1:i)//' '//datim(1:24)
  call flush(io)
end subroutine PrtDate


!-----------------------------------------------------------
! GDate
!-----------------------------------------------------------
!*Deck GDate
Subroutine GDate(Date1)
  Implicit Integer(A-Z)
  !C
  !C     This wrapper routine either calls FDate (on bsd systems) or
  !C     builds the 24-character date in some other way.
  !C
  Character*(*) Date1
  !C
  !C#ifdef IBM_RS6K
  !C#define GDATE_doNE
  !C      Character*26 LDate
  !C      LDate = ' '
  !C      Junk = GCTime(LDate)
  !C      Date1 = LDate
  !C#endif
  !C#ifndef GDATE_doNE
  Call FDate(Date1)
  !C#endif
  If(Len(Date1).gt.24) Date1(25:) = ' '
  Return
end Subroutine GDate



!-----------------------------------------------------------
! CopyDMat
!-----------------------------------------------------------
! 2010.10.26 Copy matrix(double precision)
!-----------------------------------------------------------
subroutine CopyDMat(fromMat,toMat,n)
  implicit none
  integer n,i,j
  double precision fromMat(n,n),toMat(n,n)

  do i=1,n
     do j=1,n
        toMat(j,i)=fromMat(j,i)
     enddo
  enddo

end subroutine CopyDMat


!-----------------------------------------------------------
! BNDANG
!-----------------------------------------------------------
! Return bond ang
!-----------------------------------------------------------
SUBROUTINE BNDANG(I,IA,IB,ANGLE)
  use allmod
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)


  ! COMPUTES THE ANGLE IN RADIANS FORMED BY ATOMS I-IA-IB.  RETURNED
  ! IN ANGLE.


  ! FORM THE BOND VECTORS WITH VERTEX IA.

  X1 = XYZ(1,I) - XYZ(1,IA)
  Y1 = XYZ(2,I) - XYZ(2,IA)
  Z1 = XYZ(3,I) - XYZ(3,IA)
  X2 = XYZ(1,IB) - XYZ(1,IA)
  Y2 = XYZ(2,IB) - XYZ(2,IA)
  Z2 = XYZ(3,IB) - XYZ(3,IA)
  VNORM1 = DSQRT(X1*X1 + Y1*Y1 + Z1*Z1)
  VNORM2 = DSQRT(X2*X2 + Y2*Y2 + Z2*Z2)
  COSINE = (X1*X2 + Y1*Y2 + Z1*Z2)/(VNORM1*VNORM2)
  IF(ABS(COSINE) > 1.0D0) COSINE = SIGN(1.0D0,COSINE)
  ANGLE = ACOS(COSINE)
  RETURN
end SUBROUTINE BNDANG




!-----------------------------------------------------------
! DIHEDR
!-----------------------------------------------------------
! Return dihedral angel
! $Id: dihedr.f90,v 1.1.1.1 2007/01/26 20:22:34 ayers Exp $
!------------------------------------------------------------------------
SUBROUTINE DIHEDR(XYZ,I,IA,IB,IC,DIH)

  ! DETERMINES THE I-IA-IB-IC DIHEDRAL ANGLE IN RADIANS.  THE
  ! ANGLE DIH IS POSITIVE IF IC IS LOCATED CLOCKWISE FROM I WHEN
  ! VIEWING FROM IA THROUGH IB.

  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DIMENSION XYZ(3,*)

  ! SHIFT IA-I AND IB-IC BOND VECTORS TO A COMMON ORIGIN.

  AIIX = XYZ(1,I) - XYZ(1,IA)
  AIIY = XYZ(2,I) - XYZ(2,IA)
  AIIZ = XYZ(3,I) - XYZ(3,IA)
  BCX = XYZ(1,IC) - XYZ(1,IB)
  BCY = XYZ(2,IC) - XYZ(2,IB)
  BCZ = XYZ(3,IC) - XYZ(3,IB)

  ! FORM THE IA-IB BOND AXIS VECTOR.

  ABX = XYZ(1,IB) - XYZ(1,IA)
  ABY = XYZ(2,IB) - XYZ(2,IA)
  ABZ = XYZ(3,IB) - XYZ(3,IA)

  ! REMOVE FROM (AIIX,AIIY,AIIZ) AND (BCX,BCY,BCZ) ANY PROJECTION ALONG
  ! THE (ABX,ABY,ABZ) AXIS.

  DOT1 = AIIX*ABX + AIIY*ABY + AIIZ*ABZ
  ABSQR = ABX**2 + ABY**2 + ABZ**2
  PROJ1 = DOT1/ABSQR
  AIIX = AIIX - PROJ1*ABX
  AIIY = AIIY - PROJ1*ABY
  AIIZ = AIIZ - PROJ1*ABZ
  DOT2 = BCX*ABX + BCY*ABY + BCZ*ABZ
  PROJ2 = DOT2/ABSQR
  BCX = BCX - PROJ2*ABX
  BCY = BCY - PROJ2*ABY
  BCZ = BCZ - PROJ2*ABZ

  ! COMPUTE THE CROSS-PRODUCT (AIIX,AIIY,AIIZ) X (BCX,BCY,BCZ).  STORE
  ! IT IN THE VECTOR (AIBCX,AIBCY,AIBCZ).

  AIBCX = AIIY*BCZ - AIIZ*BCY
  AIBCY = AIIZ*BCX - AIIX*BCZ
  AIBCZ = AIIX*BCY - AIIY*BCX

  ! IF (AIBCX,AIBCY,AIBCZ) POINTS IN THE SAME DIRECTION AS
  ! (ABX,ABY,ABZ) THEN IC IS LOCATED CLOCKWISE FROM I WHEN
  ! VIEWED FROM IA TOWARD IB.  THUS, IN MOVING ALONG THE PATH
  ! I-IA-IB-IC, A CLOCKWISE OR POSITIVE ANGLE IS OBSERVED.
  ! TO DETERMINE WHETHER THESE VECTORS POINT IN THE SAME OR
  ! OPPOSITE DIRECTIONS, COMPUTE THEIR DOT PRODUCT.

  DOT3 = AIBCX*ABX + AIBCY*ABY + AIBCZ*ABZ
  DIREC = SIGN(1.0D0,DOT3)

  ! COMPUTE THE DIHEDRAL ANGLE DIH.

  DOT4 = AIIX*BCX + AIIY*BCY + AIIZ*BCZ
  AILENG = SQRT(AIIX**2 + AIIY**2 + AIIZ**2)
  BCLENG = SQRT(BCX**2 + BCY**2 + BCZ**2)
  IF(ABS(AILENG*BCLENG) < 1.0D-5)THEN
     COSINE = 1.0D0
  ELSE
     COSINE = DOT4/(AILENG*BCLENG)
  ENDIF
  COSINE = MAX(-1.0D0,COSINE)
  COSINE = MIN(1.0D0,COSINE)
  DIH = ACOS(COSINE)*DIREC
  RETURN
end SUBROUTINE DIHEDR


!-----------------------------------------------------------
! wrtrestart
!-----------------------------------------------------------
! Ed Brothers. August 18, 2002
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine wrtrestart
  use allmod
  implicit double precision(a-h,o-z)

  !    logical :: present
  character(len=80) keywd
  !    character(len=20) tempstring

  ! The purpose of this routine is to write out an input file based
  ! on the result of a previous calculation.  Thus the restart file
  ! should be the same as the input file, except gifferent geometry.

  open(infile,file=infilename,status='old')
  open(irstfile,file=rstfilename,status='unknown')
  write (ioutfile,'(/" WROTE A RESTART FILE. ")')

  istart = 1
  ifinal = 80
  read (infile,'(A80)') keywd
  call upcase(keywd,80)

  DO WHILE(istart.ne.0.and.ifinal.ne.0)
     write(irstfile,'(A80)') keywd
     read (infile,'(A80)') keywd
     call upcase(keywd,80)
     call rdword(keywd,ISTART,IFINAL)
  ENDDO

  ! After copying the keywords, put in a blank space.

  write(irstfile,'("  ")')

  ! Now we can write the molecular geometry.

  DO I=1,natom
     Write (irstfile,'(A2,6x,F14.9,3x,F14.9,3x,F14.9)') &
          symbol(iattype(I)),xyz(1,I)*0.529177249d0, &
          xyz(2,I)*0.529177249d0,xyz(3,I)*0.529177249d0
  ENDDO

  ! Now another blank line.

  write(irstfile,'("  ")')

  ! If this is DFT calculation, write the grid specification out.
  ! This has to be copied from the input file, as otherwise the radial
  ! number would get smaller as you restarted multiple times.  For
  ! instance, if you specify 75 radial points, the code will chop off
  ! the insignifigant ones, giving you 72 points.  If that is written to
  ! the restart file, when you use the restart file you will actually
  ! have a coarser grid.

  IF (quick_method%DFT .OR. quick_method%SEDFT) THEN
     ! skip all geometry specification and the terminating blank line.
     DO I=1,natom+1
        read (infile,'(A80)') keywd
     ENDDO

     ! Now read in the grid specifiaction and write it out.  Not that the
     ! number of lines is the number of regions+2.

     DO I=1,iregion+2
        read (infile,'(A80)') keywd
        write(irstfile,'(A80)') keywd
     ENDDO

  ENDIF

  close(infile)
  close(irstfile)

  return
end subroutine wrtrestart


!-----------------------------------------------------------
! Complete the symmetry of matrix
!-----------------------------------------------------------
! if we have symmetry matrix but only get right-top part infos
! we can simply copy it to left-bottom
!-----------------------------------------------------------
subroutine copySym(O,n)
implicit none
integer n,i,j
double precision O(n,n)

do i=1,n
    do j=i+1,n
        O(i,j)=O(j,i)
    enddo
enddo

end subroutine copySym


!-----------------------------------------------------------
! sum2Mat
!-----------------------------------------------------------
! sum2Mat=sigma(i,j) Mat1(i,j)*Mat2(j,i)
!-----------------------------------------------------------
function sum2Mat(Mat1,Mat2,n)
implicit none
integer n,i,j
double precision Mat1(n,n),Mat2(n,n),sum2Mat

sum2Mat=0.0d0
do j=1,n
    do i=1,n
        sum2Mat=sum2Mat+Mat1(i,j)*Mat2(j,i)
    enddo
enddo
end function sum2Mat


!-----------------------------------------------------------
! IMatMul
!-----------------------------------------------------------
!multiplication for two matrix m3=m1xm2 [integer]
!-----------------------------------------------------------
      subroutine IMatMul(n,m1,m2,m3)
      implicit none
      integer n,i,j,k,m1(n,n),m2(n,n),m3(n,n),P

      do i=1,n; do j=1,n
         P=0
         do k=1,n
            P=P+m1(i,k)*m2(k,j)
         enddo
         m3(i,j)=P
      enddo; enddo

      end

!-----------------------------------------------------------
! DMatMul
!-----------------------------------------------------------
!multiplication for two matrix m3=m1xm2 double precision
!-----------------------------------------------------------

      subroutine DMatMul(n,m1,m2,m3)
      implicit none
      integer n,i,j,k
      double precision m1(n,n),m2(n,n),m3(n,n),P

      do j=1,n
        do i=1,n
            P=0d0
            do k=1,n
                P=P+m1(i,k)*m2(k,j)
            enddo
            m3(i,j)=P
        enddo
      enddo
      return
      end subroutine DMatMul


!-----------------------------------------------------------
! greedy_distrubutie
!-----------------------------------------------------------
! Yipu Miao 11/19/2010
! greedy algrithm to obtain optimized distrubution
!-----------------------------------------------------------    
    subroutine greedy_distrubute(j,n,node,node_jn,node_j)
    implicit none
    integer n               ! no. of stuff
    integer j(n)            ! Value of n stuff
    integer jorder(n)       ! value list after sort
    integer node            ! no. of nodes
    integer node_jn(0:node-1)   ! total no. of stuff one node takes
    integer node_j(0:node-1,n)  ! the nth stuff one node takes
    integer i,jn,k,jj,ii,maxnode,minnode
    integer max_val,min_val
    integer node_jtot(0:node-1) ! total value of one node
    logical isUsed(n)

    do i=0,node-1
        node_jn(i)=0
        node_jtot(i)=0
    enddo
    do i=1,n
        isUsed(i)=.false.
    enddo
  
    minnode=0
    ! The basis idea is to put the largest element to the group with fewest 
    ! value
    do i=1,n
        ! first find the most valuable element
        max_val=0
        do jj=1,n
           if((j(jj).ge.max_val).and.(.not.isUsed(jj))) then
              ii=jj
              max_val=j(jj)
           endif
        enddo
        isUsed(ii)=.true.

        ! then put it to the group with minimum value
        node_jn(minnode)=node_jn(minnode)+1
        node_j(minnode,node_jn(minnode))=ii
        node_jtot(minnode)=node_jtot(minnode)+j(ii)
        
        ! find now which group is the most valueless
        do jn=0,node-1
            if (node_jtot(jn).lt.node_jtot(minnode)) minnode=jn
        enddo
    enddo

    
!    do i=0,node-1
!       write(*,*) "I=",i,node_jtot(i),node_jn(i)
!       do jn=1,node_jn(i)
!           write(*,*) i,jn,node_j(i,jn)
!       enddo
!   enddo
    end subroutine greedy_distrubute
    
    