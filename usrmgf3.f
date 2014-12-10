************************************************************************
*                                                                      *
      subroutine usrmgf3(xx,yy,zz,bxy,btot,
     &                   bbx,bby,bbz)
*                                                                      *
*        sample subroutine for user defined magnetic field.            *
*                                                                      *
*        Strength is given by mgf[T/m^2]                               *
*        Strength of additional field is given by gap[T]               *
*                                                                      *
*        input :                                                       *
*           xx, yy, zz    : position [cm]                              *
*           cmg           : strength of magnetic field [T/cm^2]        *
*           dmg           : additinal magnetic field [cm^2]            *
*                           dmg[cm^2] = gap[T] / cmg[T/cm^2]           *
*                                                                      *
*        output :                                                      *
*           bbx, bby, bbz : magnetic field [cm^2]                      *
*                           M[T] = cmg[T/cm^2] * bbx[cm^2]             *
*		 modified by Ye,YANG										   *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      double precision xcc, ycc, zcc, bxc, byc, bzc, btotc
      common /rzforc/ xcc, ycc, zcc, bxc, byc, bzc, btotc

*-----------------------------------------------------------------------
*		 dipole magnets
*		 output the magnetic field vector
*-----------------------------------------------------------------------
			
			xcc = xx * 10
			ycc = yy * 10
			zcc = zz * 10
			
			call readfieldfrommap()
			
			btot = sqrt( bxc**2 + byc**2 + bzc**2 )
			bxy = sqrt( bxc**2 + byc**2)
			bbx = bxc
			bby = byc
			bbz = bzc
			
!			if (btot .ne. 0) then
!				bbx = bxc
!				bby = byc
!				bbz = bzc
!			else
!				bbx = 0.0d0
!				bby = 0.0d0
!				bbz = 1.0d0
!			end if


      return
      end

