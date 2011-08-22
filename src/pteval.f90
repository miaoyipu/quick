! Ed Brothers. January 17, 2002
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine pteval(gridx,gridy,gridz,phi,dphidx,dphidy,dphidz,Iphi)
    use allmod
    implicit double precision(a-h,o-z)

! Given a point in space, this function calculates the value of basis
! function I and the value of its cartesian derivatives in all three
! derivatives.

    x1=(gridx-xyz(1,ncenter(Iphi)))
    y1=(gridy-xyz(2,ncenter(Iphi)))
    z1=(gridz-xyz(3,ncenter(Iphi)))
    rsquared=x1*x1+y1*y1+z1*z1


    phi=0.d0
    dphidx=0.d0
    dphidy=0.d0
    dphidz=0.d0
    IF (rsquared > sigrad2(Iphi)) THEN
        continue
    ELSE
        IF (itype(1,Iphi) == 0) THEN
            x1imin1=0.d0
            x1i=1.d0
            x1iplus1=x1
        ELSE
            x1imin1=x1**(itype(1,Iphi)-1)
            x1i=x1imin1*x1
            x1iplus1=x1i*x1
        ENDIF

        IF (itype(2,Iphi) == 0) THEN
            y1imin1=0.d0
            y1i=1.d0
            y1iplus1=y1
        ELSE
            y1imin1=y1**(itype(2,Iphi)-1)
            y1i=y1imin1*y1
            y1iplus1=y1i*y1
        ENDIF

        IF (itype(3,Iphi) == 0) THEN
            z1imin1=0.d0
            z1i=1.d0
            z1iplus1=z1
        ELSE
            z1imin1=z1**(itype(3,Iphi)-1)
            z1i=z1imin1*z1
            z1iplus1=z1i*z1
        ENDIF

        xtype=dble(itype(1,Iphi))
        ytype=dble(itype(2,Iphi))
        ztype=dble(itype(3,Iphi))

        temp=0.d0
        DO Icon=1,ncontract(IPhi)
            temp = dcoeff(Icon,IPhi)*DExp((-aexp(Icon,IPhi))*rsquared)
            Phi=Phi+temp
            dphidx=dphidx+temp*(-2.d0*(aexp(Icon,IPhi))*x1iplus1 &
            +xtype*x1imin1)
            dphidy=dphidy+temp*(-2.d0*(aexp(Icon,IPhi))*y1iplus1 &
            +ytype*y1imin1)
            dphidz=dphidz+temp*(-2.d0*(aexp(Icon,IPhi))*z1iplus1 &
            +ztype*z1imin1)
        ENDDO
        Phi=phi*x1i*y1i*z1i
 

        dphidx=dphidx*y1i*z1i
        dphidy=dphidy*x1i*z1i
        dphidz=dphidz*x1i*y1i
    ENDIF

    return
    end subroutine pteval



