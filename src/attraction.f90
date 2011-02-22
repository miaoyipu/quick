! Ed Brothers. October 23, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    double precision function attraction(a,b,i,j,k,ii,jj,kk,Ax,Ay,Az, &
    Bx,By,Bz,Cx,Cy,Cz,Z)
    use quick_constants_module
    implicit double precision(a-h,o-z)
    dimension aux(0:20)

! Variables needed later:
!    pi=3.1415926535897932385

    g = a+b
    Px = (a*Ax + b*Bx)/g
    Py = (a*Ay + b*By)/g
    Pz = (a*Az + b*Bz)/g

    PCsquare = (Px-Cx)**2 + (Py -Cy)**2 + (Pz -Cz)**2

! The purpose of this subroutine is to calculate the nuclear attraction
! of an electron  distributed between gtfs with orbital exponents a
! and b on A and B with angular momentums defined by i,j,k (a's x, y
! and z exponents, respectively) and ii,jj,k and kk on B with the core at
! (Cx,Cy,Cz) with charge Z. m is the "order" of the integral which
! arises from the recusion relationship.

! The this is taken from the recursive relation found in Obara and Saika,
! J. Chem. Phys. 84 (7) 1986, 3963.

! The first step is generating all the necessary auxillary integrals.
! These are (0|1/rc|0)^(m) = 2 Sqrt (g/Pi) (0||0) Fm(g(Rpc)^2)
! The values of m range from 0 to i+j+k+ii+jj+kk.

    U = g* PCsquare
    Maxm = i+j+k+ii+jj+kk
    call FmT(Maxm,U,aux)
    constant = overlap(a,b,0,0,0,0,0,0,Ax,Ay,Az,Bx,By,Bz) &
    * 2.d0 * sqrt(g/Pi)
    do L = 0,maxm
        aux(L) = aux(L)*constant
    enddo

! At this point all the auxillary integrals have been calculated.
! It is now time to decompase the attraction integral to it's
! auxillary integrals through the recursion scheme.  To do this we use
! a recursive function.

    attraction = attrecurse(i,j,k,ii,jj,kk,0,aux,Ax,Ay,Az,Bx,By,Bz, &
    Cx,Cy,Cz,Px,Py,Pz,g)
    attraction = attraction*(-1.d0)* Z
    return
    end function attraction



! Ed Brothers. October 23, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    double precision recursive function &
    attrecurse(i,j,k,ii,jj,kk,m,aux,Ax,Ay, &
    Az,Bx,By,Bz,Cx,Cy,Cz,Px,Py,Pz,g) &
    result(attrec)
    implicit double precision(a-h,o-z)
    dimension iexponents(6),center(12),aux(0:20)

! The this is taken from the recursive relation found in Obara and Saika,
! J. Chem. Phys. 84 (7) 1986, 3963.

! If this is one of the auxillary integrals (s||s)^(m), assign value and
! return.

    if (i+j+k+ii+jj+kk == 0) then
        attrec=aux(m)

    ! Otherwise, use the recusion relation from Obara and Saika.  The first
    ! step is to find the lowest nonzero angular momentum exponent.  This is
    ! because the more exponents equal zero the fewer terms need to be
    ! calculated, and each recursive loop reduces the angular momentum
    ! exponents. This therefore reorders the atoms and sets the exponent
    ! to be reduced.

    else
        iexponents(1) = i
        iexponents(2) = j
        iexponents(3) = k
        iexponents(4) = ii
        iexponents(5) = jj
        iexponents(6) = kk
        center(7) = Cx
        center(8) = Cy
        center(9) = Cz
        center(10)= Px
        center(11)= Py
        center(12)= Pz
        ilownum=300
        ilowex=300
        do L=1,6
            if (iexponents(L) < ilowex .AND. iexponents(L) /= 0) then
                ilowex=iexponents(L)
                ilownum=L
            endif
        enddo
        if (ilownum <= 3) then
            center(1)=Ax
            center(2)=Ay
            center(3)=Az
            center(4)=Bx
            center(5)=By
            center(6)=Bz
        else
            center(4)=Ax
            center(5)=Ay
            center(6)=Az
            center(1)=Bx
            center(2)=By
            center(3)=Bz
            iexponents(4) = i
            iexponents(5) = j
            iexponents(6) = k
            iexponents(1) = ii
            iexponents(2) = jj
            iexponents(3) = kk
            ilownum = ilownum - 3
        endif

    ! The first step is lowering the orbital exponent by one.

        iexponents(ilownum) = iexponents(ilownum)-1

    ! At this point, calculate the first two terms of the recusion
    ! relation.

        attrec=0.d0
        PA = center(9+ilownum)-center(ilownum)
        if (PA /= 0) attrec = attrec + PA * &
        attrecurse(iexponents(1),iexponents(2), &
        iexponents(3),iexponents(4), &
        iexponents(5),iexponents(6), &
        m,aux, &
        center(1),center(2),center(3), &
        center(4),center(5),center(6), &
        center(7),center(8),center(9), &
        center(10),center(11),center(12),g)

        PC = center(9+ilownum)-center(6+ilownum)
        if (PC /= 0) attrec = attrec - PC * &
        attrecurse(iexponents(1),iexponents(2), &
        iexponents(3),iexponents(4), &
        iexponents(5),iexponents(6), &
        m+1,aux, &
        center(1),center(2),center(3), &
        center(4),center(5),center(6), &
        center(7),center(8),center(9), &
        center(10),center(11),center(12),g)

    ! The next two terms only arise is the angual momentum of the dimension
    ! of A that has already been lowered is not zero.  In other words, if a
    ! (px|1/rc|px) was passed to this subroutine, we are now considering
    ! (s|1/rc|px), and the following term does not arise, as the x expoent
    ! on A is zero.

        if (iexponents(ilownum) /= 0) then
            coeff = dble(iexponents(ilownum))/(2.d0*g)
            iexponents(ilownum) = iexponents(ilownum)-1
            attrec = attrec + coeff*( &
            attrecurse(iexponents(1),iexponents(2), &
            iexponents(3),iexponents(4), &
            iexponents(5),iexponents(6), &
            m,aux, &
            center(1),center(2),center(3), &
            center(4),center(5),center(6), &
            center(7),center(8),center(9), &
            center(10),center(11),center(12),g) &
            -attrecurse(iexponents(1),iexponents(2), &
            iexponents(3),iexponents(4), &
            iexponents(5),iexponents(6), &
            m+1,aux, &
            center(1),center(2),center(3), &
            center(4),center(5),center(6), &
            center(7),center(8),center(9), &
            center(10),center(11),center(12),g) &
            )
            iexponents(ilownum) = iexponents(ilownum)+1
        endif

    ! The next two terms only arise is the angual momentum of the dimension
    ! of A that has already been lowered is not zero in B.  If a
    ! (px|1/rc|px) was passed to this subroutine, we are now considering
    ! (s|1/rc|px), and the following term does arise, as the x exponent on
    ! B is 1.

        if (iexponents(ilownum+3) /= 0) then
            coeff = dble(iexponents(ilownum+3))/(2.d0*g)
            iexponents(ilownum+3) = iexponents(ilownum+3)-1
            attrec = attrec + coeff*( &
            attrecurse(iexponents(1),iexponents(2), &
            iexponents(3),iexponents(4), &
            iexponents(5),iexponents(6), &
            m,aux, &
            center(1),center(2),center(3), &
            center(4),center(5),center(6), &
            center(7),center(8),center(9), &
            center(10),center(11),center(12),g) &
            -attrecurse(iexponents(1),iexponents(2), &
            iexponents(3),iexponents(4), &
            iexponents(5),iexponents(6), &
            m+1,aux, &
            center(1),center(2),center(3), &
            center(4),center(5),center(6), &
            center(7),center(8),center(9), &
            center(10),center(11),center(12),g) &
            )
        endif
    endif

    return
    end