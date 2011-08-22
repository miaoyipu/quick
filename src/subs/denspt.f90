! Ed Brothers. January 17, 2002
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine denspt(gridx,gridy,gridz,densitya,densityb,gax,gay, &
    gaz,gbx,gby,gbz)
    use allmod
    implicit double precision(a-h,o-z)

! Given a point in space, this function calculates the densities and
! gradient  at that point.  The gradients are stored in the common block
! three element arrays ga and gb for alpha and beta electron gradients. Thus
! the x component of the alpha density is stored in ga(1).

    densitya=0.d0
!    densityb=0.d0

    gax=0.d0
    gay=0.d0
    gaz=0.d0
!    gbx=0.d0
!    gby=0.d0
!    gbz=0.d0

    do Ibas=1,nbasis
        DENSEIJ=quick_qm_struct%dense(Ibas,Ibas)
        If(DABS(DENSEIJ) < quick_method%DMCutoff)then
          continue
        else
        DENSEBIJ=quick_qm_struct%denseB(Ibas,Ibas)
!        call pteval(gridx,gridy,gridz,phi,dphidx,dphidy,dphidz,Ibas)
        phi=phixiao(Ibas)
        dphidx=dphidxxiao(Ibas)
        dphidy=dphidyxiao(Ibas)
        dphidz=dphidzxiao(Ibas)
                                quicktest = DABS(dphidx+dphidy+dphidz+ &
                                phi)
                                if (quicktest < quick_method%DMCutoff ) then
                                    continue
                                else

        densitya=densitya+DENSEIJ*phi*phi/2.0d0
!        densityb=densityb+DENSEBIJ*phi*phi
        gax=gax+DENSEIJ*phi*dphidx
        gay=gay+DENSEIJ*phi*dphidy
        gaz=gaz+DENSEIJ*phi*dphidz
!        gbx=gbx+DENSEBIJ*2.d0*phi*dphidx
!        gby=gby+DENSEBIJ*2.d0*phi*dphidy
!        gbz=gbz+DENSEBIJ*2.d0*phi*dphidz

        do Jbas=Ibas+1,nbasis
            DENSEIJ=quick_qm_struct%dense(Jbas,Ibas)
!            DENSEBIJ=quick_qm_struct%denseB(Jbas,Ibas)
!            call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy,dphi2dz,Jbas)
        phi2=phixiao(Jbas)
        dphi2dx=dphidxxiao(Jbas)
        dphi2dy=dphidyxiao(Jbas)
        dphi2dz=dphidzxiao(Jbas)
!                                quicktest = DABS( &
!                                phi2)
!                                if (quicktest < tol ) then
!                                    continue
!                                else

            densitya=densitya+DENSEIJ*phi*phi2
!            densityb=densityb+2.d0*DENSEBIJ*phi*phi2
            gax=gax+DENSEIJ*(phi*dphi2dx+phi2*dphidx)
            gay=gay+DENSEIJ*(phi*dphi2dy+phi2*dphidy)
            gaz=gaz+DENSEIJ*(phi*dphi2dz+phi2*dphidz)
!            gbx=gbx+DENSEBIJ*2.d0*(phi*dphi2dx+phi2*dphidx)
!            gby=gby+DENSEBIJ*2.d0*(phi*dphi2dy+phi2*dphidy)
!            gbz=gbz+DENSEBIJ*2.d0*(phi*dphi2dz+phi2*dphidz)
!                               endif
        enddo
                                       endif
      endif
    enddo

!    if( .NOT. quick_method%UNRST) then
!        densitya=densitya/2.d0
        densityb=densitya
!        gax =gax*.5d0
!        gay =gay*.5d0
!        gaz =gaz*.5d0
        gbx =gax
        gby =gay
        gbz =gaz
!    endif

    END subroutine denspt

