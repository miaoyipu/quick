! Ed Brothers. October 3, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    function overlap(a,b,i,j,k,ii,jj,kk,Ax,Ay,Az,Bx, &
    By,Bz)
    use quick_constants_module
    implicit double precision(a-h,o-z)

! The purpose of this subroutine is to calculate the overlap between
! two normalized gaussians. i,j and k are the x,y,
! and z exponents for the gaussian with exponent a, and ii,jj, and kk
! have the same order for b.

! Constants:

!    pi=3.1415926535897932385
!    pito3half=5.568327996831707845284817982118835702014
!    pito3half = pi**(1.5)

! The first step is to see if this function is zero due to symmetry.
! If it is not, reset overlap to 0.

    overlap = (1+(-1)**(i+ii))*(1+(-1)**(j+jj))*(1+(-1)**(k+kk)) &
    +(Ax-Bx)**2 + (Ay-By)**2 + (Az-Bz)**2
    if (overlap == 0.d0) goto 100
    overlap=0.d0


! If it is not zero, construct P and g values.  The gaussian product
! theory states the product of two s gaussians on centers A and B
! with exponents a and b forms a new s gaussian on P with exponent
! g.  (g comes from gamma, as is "alpha,beta, gamma" and P comes
! from "Product." Also needed are the PA differences.

    g = a+b
    Px = (a*Ax + b*Bx)/g
    Py = (a*Ay + b*By)/g
    Pz = (a*Az + b*Bz)/g

    PAx= Px-Ax
    PAy= Py-Ay
    PAz= Pz-Az
    PBx= Px-Bx
    PBy= Py-By
    PBz= Pz-Bz

! There is also a few factorials that are needed in the integral many
! times.  Calculate these as well.

    xnumfact=fact(i)*fact(ii)*fact(j)*fact(jj)*fact(k)*fact(kk)

! Now start looping over i,ii,j,jj,k,kk to form all the required elements.

    do iloop=0,i
        do iiloop=0,ii
            do jloop=0,j
                do jjloop=0,jj
                    do kloop=0,k
                        do kkloop=0,kk
                            ix=iloop+iiloop
                            jy=jloop+jjloop
                            kz=kloop+kkloop

                        ! Check to see if this element is zero.

                            element=(1+(-1)**(ix))*(1+(-1)**(jy))*(1+(-1)**(kz))/8
                            if (element == 0) goto 50

                        ! Continue calculating the elements.  The next elements arise from the
                        ! different angular momentums portion of the GPT.

                            element=PAx**(i-iloop) &
                            *PBx**(ii-iiloop) &
                            *PAy**(j-jloop) &
                            *PBy**(jj-jjloop) &
                            *PAz**(k-kloop) &
                            *PBz**(kk-kkloop) &
                            *xnumfact &
                            /(fact(iloop)*fact(iiloop)* &
                            fact(jloop)*fact(jjloop)* &
                            fact(kloop)*fact(kkloop)* &
                            fact(i-iloop)*fact(ii-iiloop)* &
                            fact(j-jloop)*fact(jj-jjloop)* &
                            fact(k-kloop)*fact(kk-kkloop))

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element*g**(dble(-3 - ix - jy - kz)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,ix/2
                                element = element * (dble(ix)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,jy/2
                                element = element * (dble(jy)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,kz/2
                                element = element * (dble(kz)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                            50 continue
                            overlap = overlap + element
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

! The final step is multiplying in the K factor (from the gpt)

    overlap = overlap* Exp(-((a*b*((Ax - Bx)**2.d0 + (Ay - By)**2.d0 &
    + (Az - Bz)**2.d0))/(a + b)))

    100 continue
    return
    end function overlap

    function ssoverlap(a,b,Ax,Ay,Az,Bx,By,Bz)
    use quick_constants_module
    implicit double precision(a-h,o-z)

!    pito3half=5.568327996831707845284817982118835702014
!    pito3half = pi**1.5
    g = a+b

    ssoverlap = pito3half*1.d0*g**(-3.d0/2.d0)

    ssoverlap = ssoverlap* Exp(-((a*b*((Ax - Bx)**2.d0 + &
    (Ay - By)**2.d0 &
    + (Az - Bz)**2.d0))/(a + b)))

    end function ssoverlap


! Ed Brothers. October 3, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine overlapone(a,b,i,j,k,ii,jj,kk,Ax,Ay,Az,Bx, &
    By,Bz,fmmtemparray)
    use quick_constants_module
    implicit double precision(a-h,o-z)

           real*8 Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Px,Py,Pz,g
           real*8 AA(3),BB(3),CC(3),PP(3)
           real*8 fmmtemparray(0:2,0:2,1:2)
!           common /xiaofmm/fmmonearray,AA,BB,CC,PP,g

! The purpose of this subroutine is to calculate the overlap between
! two normalized gaussians. i,j and k are the x,y,
! and z exponents for the gaussian with exponent a, and ii,jj, and kk
! have the same order for b.

! Constants:

!    pi=3.1415926535897932385
!    pito3half=5.568327996831707845284817982118835702014
!    pito3half = pi**(1.5)

! The first step is to see if this function is zero due to symmetry.
! If it is not, reset overlap to 0.

! If it is not zero, construct P and g values.  The gaussian product
! theory states the product of two s gaussians on centers A and B
! with exponents a and b forms a new s gaussian on P with exponent
! g.  (g comes from gamma, as is "alpha,beta, gamma" and P comes
! from "Product." Also needed are the PA differences.

    g = a+b
    Px = (a*Ax + b*Bx)/g
    Py = (a*Ay + b*By)/g
    Pz = (a*Az + b*Bz)/g

    PAx= Px-Ax
    PAy= Py-Ay
    PAz= Pz-Az
    PBx= Px-Bx
    PBy= Py-By
    PBz= Pz-Bz

    Y1=0.0d0
    Y2=0.0d0
    Y3=0.0d0
    Y4=0.0d0

! There is also a few factorials that are needed in the integral many
! times.  Calculate these as well.

    xnumfact=fact(i)*fact(ii)*fact(j)*fact(jj)*fact(k)*fact(kk)

! Now start looping over i,ii,j,jj,k,kk to form all the required elements.

    do iloop=0,i
        do iiloop=0,ii
            do jloop=0,j
                do jjloop=0,jj
                    do kloop=0,k
                        do kkloop=0,kk
                            ix=iloop+iiloop
                            jy=jloop+jjloop
                            kz=kloop+kkloop

                        ! Continue calculating the elements.  The next elements arise from the
                        ! different angular momentums portion of the GPT.

                            element00=PAx**(i-iloop) &
                            *PBx**(ii-iiloop) &
                            *PAy**(j-jloop) &
                            *PBy**(jj-jjloop) &
                            *PAz**(k-kloop) &
                            *PBz**(kk-kkloop) &
                            *xnumfact &
                            /(fact(iloop)*fact(iiloop)* &
                            fact(jloop)*fact(jjloop)* &
                            fact(kloop)*fact(kkloop)* &
                            fact(i-iloop)*fact(ii-iiloop)* &
                            fact(j-jloop)*fact(jj-jjloop)* &
                            fact(k-kloop)*fact(kk-kkloop))

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix))*(1+(-1)**(jy))*(1+(-1)**(kz))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 50
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix - jy - kz)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,ix/2
                                element = element * (dble(ix)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,jy/2
                                element = element * (dble(jy)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,kz/2
                                element = element * (dble(kz)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                            Y1=Y1+element
                            50 continue
!                            fmmtemparray(0,0,1) = fmmtemparray(0,0,1) + element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix))*(1+(-1)**(jy))*(1+(-1)**(kz+1))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 60
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix - jy - kz-1)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,ix/2
                                element = element * (dble(ix)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,jy/2
                                element = element * (dble(jy)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,(kz+1)/2
                                element = element * (dble(kz+1)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                            Y2=Y2+element
                            60 continue
!                            fmmtemparray(1,0,1) = fmmtemparray(1,0,1) + element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix+1))*(1+(-1)**(jy))*(1+(-1)**(kz))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 70
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix -1 - jy - kz)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,(ix+1)/2
                                element = element * (dble(ix+1)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,jy/2
                                element = element * (dble(jy)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,kz/2
                                element = element * (dble(kz)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.
     
                            Y3=Y3+element
                            70 continue
!                            fmmtemparray(1,1,1) = fmmtemparray(1,1,1) + element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix))*(1+(-1)**(jy+1))*(1+(-1)**(kz))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 80
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix - jy -1 - kz)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,ix/2
                                element = element * (dble(ix)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,(jy+1)/2
                                element = element * (dble(jy+1)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,kz/2
                                element = element * (dble(kz)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                            Y4=Y4+element
                            80 continue
!                            fmmtemparray(1,1,2) = fmmtemparray(1,1,2) + element


                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

! The final step is multiplying in the K factor (from the gpt)

!    overlap = overlap* Exp(-((a*b*((Ax - Bx)**2.d0 + (Ay - By)**2.d0 &
!    + (Az - Bz)**2.d0))/(a + b)))

    fmmtemparray(0,0,1)=Y1
    fmmtemparray(1,0,1)=Y2
    fmmtemparray(1,1,1)=Y3
    fmmtemparray(1,1,2)=Y4

    100 continue
    return
    end 

! Ed Brothers. October 3, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine overlaptwo(a,b,i,j,k,ii,jj,kk,Ax,Ay,Az,Bx, &
    By,Bz,fmmonearray)
    use quick_constants_module
    implicit double precision(a-h,o-z)

           real*8 Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Px,Py,Pz,g
           real*8 AA(3),BB(3),CC(3),PP(3)
           real*8 fmmonearray(0:2,0:2,1:2)
!           common /xiaofmm/fmmonearray,AA,BB,CC,PP,g

! The purpose of this subroutine is to calculate the overlap between
! two normalized gaussians. i,j and k are the x,y,
! and z exponents for the gaussian with exponent a, and ii,jj, and kk
! have the same order for b.

! Constants:

!    pi=3.1415926535897932385
!    pito3half=5.568327996831707845284817982118835702014
!    pito3half = pi**(1.5)

! The first step is to see if this function is zero due to symmetry.
! If it is not, reset overlap to 0.

! If it is not zero, construct P and g values.  The gaussian product
! theory states the product of two s gaussians on centers A and B
! with exponents a and b forms a new s gaussian on P with exponent
! g.  (g comes from gamma, as is "alpha,beta, gamma" and P comes
! from "Product." Also needed are the PA differences.

    g = a+b
    Px = (a*Ax + b*Bx)/g
    Py = (a*Ay + b*By)/g
    Pz = (a*Az + b*Bz)/g

    PAx= Px-Ax
    PAy= Py-Ay
    PAz= Pz-Az
    PBx= Px-Bx
    PBy= Py-By
    PBz= Pz-Bz

    Y1=0.0d0
    Y2=0.0d0
    Y3=0.0d0
    Y4=0.0d0
    Y5=0.0d0
    Y6=0.0d0
    Y7=0.0d0
    Y8=0.0d0
    Y9=0.0d0

! There is also a few factorials that are needed in the integral many
! times.  Calculate these as well.

    xnumfact=fact(i)*fact(ii)*fact(j)*fact(jj)*fact(k)*fact(kk)

! Now start looping over i,ii,j,jj,k,kk to form all the required elements.

    do iloop=0,i
        do iiloop=0,ii
            do jloop=0,j
                do jjloop=0,jj
                    do kloop=0,k
                        do kkloop=0,kk
                            ix=iloop+iiloop
                            jy=jloop+jjloop
                            kz=kloop+kkloop

                        ! Continue calculating the elements.  The next elements arise from the
                        ! different angular momentums portion of the GPT.

                            element00=PAx**(i-iloop) &
                            *PBx**(ii-iiloop) &
                            *PAy**(j-jloop) &
                            *PBy**(jj-jjloop) &
                            *PAz**(k-kloop) &
                            *PBz**(kk-kkloop) &
                            *xnumfact &
                            /(fact(iloop)*fact(iiloop)* &
                            fact(jloop)*fact(jjloop)* &
                            fact(kloop)*fact(kkloop)* &
                            fact(i-iloop)*fact(ii-iiloop)* &
                            fact(j-jloop)*fact(jj-jjloop)* &
                            fact(k-kloop)*fact(kk-kkloop))

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix))*(1+(-1)**(jy))*(1+(-1)**(kz))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 50
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix - jy - kz)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,ix/2
                                element = element * (dble(ix)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,jy/2
                                element = element * (dble(jy)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,kz/2
                                element = element * (dble(kz)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                            Y1=Y1+element
                            50 continue
!                            fmmonearray(0,0,1) = fmmonearray(0,0,1) + element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix))*(1+(-1)**(jy))*(1+(-1)**(kz+1))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 60
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix - jy - kz-1)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,ix/2
                                element = element * (dble(ix)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,jy/2
                                element = element * (dble(jy)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,(kz+1)/2
                                element = element * (dble(kz+1)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                            Y2=Y2+element
                            60 continue
!                            fmmonearray(1,0,1) = fmmonearray(1,0,1) + element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix+1))*(1+(-1)**(jy))*(1+(-1)**(kz))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 70
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix -1 - jy - kz)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,(ix+1)/2
                                element = element * (dble(ix+1)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,jy/2
                                element = element * (dble(jy)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,kz/2
                                element = element * (dble(kz)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                            Y3=Y3+element
                            70 continue
!                            fmmonearray(1,1,1) = fmmonearray(1,1,1) + element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix))*(1+(-1)**(jy+1))*(1+(-1)**(kz))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 80
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix - jy -1 - kz)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,ix/2
                                element = element * (dble(ix)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,(jy+1)/2
                                element = element * (dble(jy+1)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,kz/2
                                element = element * (dble(kz)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                            Y4=Y4+element
                            80 continue
!                            fmmonearray(1,1,2) = fmmonearray(1,1,2) + element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix))*(1+(-1)**(jy))*(1+(-1)**(kz+2))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 90
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix - jy - kz -2)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,ix/2
                                element = element * (dble(ix)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,jy/2
                                element = element * (dble(jy)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,(kz+2)/2
                                element = element * (dble(kz+2)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                            Y5=Y5+element
                            90 continue
!                            fmmonearray(2,0,1) = fmmonearray(2,0,1) + element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix+2))*(1+(-1)**(jy))*(1+(-1)**(kz))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 100
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix -2 - jy - kz)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,(ix+2)/2
                                element = element * (dble(ix+2)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,jy/2
                                element = element * (dble(jy)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,kz/2
                                element = element * (dble(kz)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                           Y5=Y5 - 0.5d0*element
                           Y8=Y8 + dsqrt(0.75d0)*element 
                           100 continue
!                            fmmonearray(2,0,1) = fmmonearray(2,0,1) - 0.5d0*element
!                            fmmonearray(2,2,1) = fmmonearray(2,2,1) + dsqrt(0.75d0)*element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix))*(1+(-1)**(jy+2))*(1+(-1)**(kz))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 110
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix - jy -2 - kz)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,ix/2
                                element = element * (dble(ix)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,(jy+2)/2
                                element = element * (dble(jy+2)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,kz/2
                                element = element * (dble(kz)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                            Y5=Y5 - 0.5d0*element
                            Y8=Y8 - dsqrt(0.75d0)*element

                            110 continue
!                            fmmonearray(2,0,1) = fmmonearray(2,0,1) - 0.5d0*element
!                            fmmonearray(2,2,1) = fmmonearray(2,2,1) - dsqrt(0.75d0)*element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix+1))*(1+(-1)**(jy))*(1+(-1)**(kz+1))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 120
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix - jy - kz -2)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,(ix+1)/2
                                element = element * (dble(ix+1)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,jy/2
                                element = element * (dble(jy)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,(kz+1)/2
                                element = element * (dble(kz+1)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                           Y6=Y6 + dsqrt(3.0d0)*element
                           120 continue
!                            fmmonearray(2,1,1) = fmmonearray(2,1,1) + dsqrt(3.0d0)*element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix))*(1+(-1)**(jy+1))*(1+(-1)**(kz+1))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 130
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix - jy - kz -2)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,ix/2
                                element = element * (dble(ix)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,(jy+1)/2
                                element = element * (dble(jy+1)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,(kz+1)/2
                                element = element * (dble(kz+1)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                           Y7=Y7+dsqrt(3.0d0)*element
                           130 continue
!                            fmmonearray(2,1,2) = fmmonearray(2,1,2) + dsqrt(3.0d0)*element

                        ! Check to see if this element is zero.

                            elementtemp=(1+(-1)**(ix+1))*(1+(-1)**(jy+1))*(1+(-1)**(kz))/8
                            if (elementtemp == 0)then
                               element=0.0d0
                               goto 140
                            endif

                        ! The next part arises from the integratation of a gaussian of arbitrary
                        ! angular momentum.

                            element=element00*g**(dble(-3 - ix - jy - kz -2)/2.d0)

                        ! Before the Gamma function code, a quick note. All gamma functions are
                        ! of the form:
                        ! 1
                        ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
                        ! 2
                        
                        ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
                        ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gamma
                        ! just requires a loop and multiplying by Pi^3/2

                            do iG=1,(ix+1)/2
                                element = element * (dble(ix+1)/2.d0-dble(iG) + .5d0)
                            enddo
                            do jG=1,(jy+1)/2
                                element = element * (dble(jy+1)/2.d0-dble(jG) + .5d0)
                            enddo
                            do kG=1,kz/2
                                element = element * (dble(kz)/2.d0-dble(kG) + .5d0)
                            enddo
                            element=element*pito3half

                        ! Now sum the whole thing into the overlap.

                           Y9=Y9 + dsqrt(3.0d0)*element

                           140 continue
!                            fmmonearray(2,2,2) = fmmonearray(2,2,2) + dsqrt(3.0d0)*element


                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

! The final step is multiplying in the K factor (from the gpt)

!    overlap = overlap* Exp(-((a*b*((Ax - Bx)**2.d0 + (Ay - By)**2.d0 &
!    + (Az - Bz)**2.d0))/(a + b)))

!    100 continue

    fmmonearray(0,0,1)=Y1
    fmmonearray(1,0,1)=Y2
    fmmonearray(1,1,1)=Y3
    fmmonearray(1,1,2)=Y4
    fmmonearray(2,0,1)=Y5
    fmmonearray(2,1,1)=Y6
    fmmonearray(2,1,2)=Y7
    fmmonearray(2,2,1)=Y8
    fmmonearray(2,2,2)=Y9

    return
    end 

    subroutine overlapzero(a,b,fmmtemparray)
    use quick_constants_module
    implicit double precision(a-h,o-z)

           real*8 Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Px,Py,Pz,g
           real*8 AA(3),BB(3),CC(3),PP(3)
           real*8 fmmtemparray(0:2,0:2,1:2)
!           common /xiaofmm/fmmonearray,AA,BB,CC,PP,g

!    pito3half=5.568327996831707845284817982118835702014
!    pito3half = pi**1.5
    g = a+b

    ssoverlap = pito3half*g**(-1.5d0)

!    ssoverlap = ssoverlap* Exp(-((a*b*((Ax - Bx)**2.d0 + &
!    (Ay - By)**2.d0 &
!    + (Az - Bz)**2.d0))/(a + b)))

    fmmtemparray(0,0,1)=ssoverlap

    end


