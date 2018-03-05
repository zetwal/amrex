
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_COORDSYS_F.H"
#include "AMReX_ArrayLim.H"

#define SDIM 4

! :: ----------------------------------------------------------
! :: SETVOL
! ::             Compute the volume of each cell
! ::
! :: INPUTS / OUTPUTS:
! ::  vol         <=  volume array
! ::  vlo,vhi      => index limits of vol array
! ::  offset       => shift to origin of computational domain
! ::  dx           => cell size
! ::  coord        => coordinate flag (0 = cartesian, 1 = RZ)
! :: ----------------------------------------------------------

  subroutine FORT_SETVOL(DIMS(reg),vol,DIMS(vol),offset,dx,coord)

    implicit none

    integer    DIMDEC(reg)
    integer    DIMDEC(vol)
    integer    coord
    REAL_T     dx(SDIM), offset(SDIM)
    REAL_T     vol(DIMV(vol))
       
    integer    i, j, k, m, n, p
    REAL_T     v
       
    if (coord .eq. 0) then

       ! cartesian

       v = dx(1)*dx(2)*dx(3)*dx(4)*dx(5)*dx(6)

       do p = ARG_L6(reg), ARG_H6(reg)
       do n = ARG_L5(reg), ARG_H5(reg)
       do m = ARG_L4(reg), ARG_H4(reg)
       do k = ARG_L3(reg), ARG_H3(reg)
          do j = ARG_L2(reg), ARG_H2(reg)
             do i = ARG_L1(reg), ARG_H1(reg)
                vol(i,j,k,m,n,p) = v
             end do
          end do
       end do
       end do
       end do
       end do

    else

       write(6,*) "FORT_SETVOLUME not define for coord = ",coord
       call bl_abort(" ")

    end if
       
  end subroutine FORT_SETVOL

! :: ----------------------------------------------------------
! :: SETAREA
! ::             Compute the area of given cell face
! ::
! :: INPUTS / OUTPUTS:
! ::  area        <=  area array
! ::  alo,ahi      => index limits of area array
! ::  offset       => shift to origin of computational domain
! ::  dx           => cell size
! ::  coord        => coordinate flag (0 =cartesian, 1 = RZ)
! :: ----------------------------------------------------------

  subroutine FORT_SETAREA(DIMS(reg),area,DIMS(area),offset,dx,dir,coord)

    implicit none

    integer    DIMDEC(reg)
    integer    DIMDEC(area)
    integer    coord, dir
    REAL_T     dx(SDIM), offset(SDIM)
    REAL_T     area(DIMV(area))
       
    integer    i, j, k, m, n p
    REAL_T     fa

    fa = 1.d0
       
    if (coord .eq. 0) then

       if (dir .eq. 0) then
          fa = dx(2)*dx(3)*dx(4)*dx(5)*dx(6)
       else if (dir .eq. 1) then
          fa = dx(1)*dx(3)*dx(4)*dx(5)*dx(6)
       else if (dir .eq. 2) then
          fa = dx(1)*dx(2)*dx(4)*dx(5)*dx(6)
       else if (dir .eq. 3) then
          fa = dx(1)*dx(2)*dx(3)*dx(5)*dx(6)
       else if (dir .eq. 4) then
          fa = dx(1)*dx(2)*dx(3)*dx(4)*dx(6)
       else if (dir .eq. 5) then
          fa = dx(1)*dx(2)*dx(3)*dx(4)*dx(5)
       else
          write(6,*) "FORT_SETAREA: invalid dir = ",dir
          call bl_abort(" ")
       end if

       do p = ARG_L6(reg), ARG_H6(reg)
       do n = ARG_L5(reg), ARG_H5(reg)
       do m = ARG_L4(reg), ARG_H4(reg)
       do k = ARG_L3(reg), ARG_H3(reg)
          do j = ARG_L2(reg), ARG_H2(reg)
             do i = ARG_L1(reg), ARG_H1(reg)
                area(i,j,k,m,n,p) = fa
             end do
          end do
       end do
       end do
       end do
       end do

    else

       write(6,*) "FORT_SETAREA not define for coord = ",coord
       call bl_abort(" ")

    end if
       
  end subroutine FORT_SETAREA

! :: SETDLOGA
! ::             Compute  d(log(A))/dr in each cell
! ::
! :: INPUTS / OUTPUTS:
! ::  dloga        <=  dloga array
! ::  dlo,dhi      => index limits of dloga array
! ::  offset       => shift to origin of computational domain
! ::  dx           => cell size
! ::  coord        => coordinate flag (0 = cartesian, 1 = RZ)
! :: ----------------------------------------------------------


  subroutine FORT_SETDLOGA(dloga,DIMS(dloga),offset,dx,dir,coord)

    implicit none

    integer    DIMDEC(dloga)
    integer    coord
    REAL_T     dx(SDIM), offset(SDIM)
    REAL_T     dloga(DIMV(dloga))
    integer dir
       
    integer    i, j, k, m, n, p
       
    if (coord .eq. 0) then

       ! cartesian

       do p = ARG_L6(dloga), ARG_H6(dloga)
       do n = ARG_L5(dloga), ARG_H5(dloga)
       do m = ARG_L4(dloga), ARG_H4(dloga)
       do k = ARG_L3(dloga), ARG_H3(dloga)
          do j = ARG_L2(dloga), ARG_H2(dloga)
             do i = ARG_L1(dloga), ARG_H1(dloga)
                dloga(i,j,k,m,n,p) = zero
             end do
          end do
       enddo
       enddo
       enddo
       enddo

    else 

       write(6,*)' Non-cartesian not allowed in 6D '
       call bl_abort(" ")

    endif

  end subroutine FORT_SETDLOGA
