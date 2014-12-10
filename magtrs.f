*********************************************************************************
*                                                                      *
      subroutine elmtrs(elct,elcx,elcy,elcz,
     &                  bmgt,bmgx,bmgy,bmgz,chgp,
     &                  delt,dpr,udir,vdir,wdir)
*                                                                      *
*                                                                      *
*       particle transfer in void under electro magnetic field         *
*       uniform electric field and dipole magnet                       *
*       last modified by K.Niita on 2011/01/10                         *
*                                                                      *
*     input  :                                                         *
*       delt   : distance                                              *
*       elct   : electric field [MeV/nsec]                             *
*       elcx,elcy,elcz  : unit vector of electric field                *
*       bmgt   : magneric field [MeV/nsec]                             *
*       elcx,elcy,elcz  : unit vector of magnetic field                *
*       chgp   : charge state of the particle                          *
*                                                                      *
*     in common                                                        *
*                                                                      *
*       e(no)    : initial energy                                      *
*       x(no),y(no),z(no)                                              *
*                : initial position                                    *
*       u(no),v(no),w(no)                                              *
*                : initial direction                                   *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       ec(no)   : final energy                                        *
*       xc(no),yc(no),zc(no)                                           *
*                : final position                                      *
*                                                                      *
*       u(no),v(no),w(no)                                              *
*                : initial direction                                   *
*                                                                      *
*       udir,vdir,wdir                                                 *
*                : final direction                                     *
*                                                                      *
*       dpr      : final distance                                      *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA 2013/11/06
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
c      include 'mmbank.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      parameter ( pi   = 3.1415926535898d0 )
      parameter ( rlit = 29.97925d0 )

*-----------------------------------------------------------------------
*     zero electro magnetic field
*-----------------------------------------------------------------------

            if( elct .eq. 0.0d0 .and. bmgt .eq. 0.0d0 ) then

               xc(ibkxc+no) = x(ibkx+no) + delt * u(ibku+no)
               yc(ibkyc+no) = y(ibky+no) + delt * v(ibkv+no)
               zc(ibkzc+no) = z(ibkz+no) + delt * w(ibkw+no)

               udir = u(ibku+no)
               vdir = v(ibkv+no)
               wdir = w(ibkw+no)

               dpr = delt

               return

            end if

*-----------------------------------------------------------------------
*     initial coordinate and direction
*-----------------------------------------------------------------------

               xx0 = x(ibkx+no)
               yy0 = y(ibky+no)
               zz0 = z(ibkz+no)

               xx = xx0
               yy = yy0
               zz = zz0

               xp = u(ibku+no)
               yp = v(ibkv+no)
               zp = w(ibkw+no)

               ecc = e(ibke+no)
               rm0 = rtyp

*-----------------------------------------------------------------------
*      electro magnetic field
*-----------------------------------------------------------------------

            call elmgt1(xx,yy,zz,xp,yp,zp,ecc,rm0,delt,dpr,
     &                  elct,elcx,elcy,elcz,
     &                  bmgt,bmgx,bmgy,bmgz)

*-----------------------------------------------------------------------

               e(ibke+no)   = ecc
               ec(ibkec+no) = ecc

               xc(ibkxc+no) = xx
               yc(ibkyc+no) = yy
               zc(ibkzc+no) = zz

               u(ibku+no) = ( xx - xx0 ) / dpr
               v(ibkv+no) = ( yy - yy0 ) / dpr
               w(ibkw+no) = ( zz - zz0 ) / dpr

               udir = xp
               vdir = yp
               wdir = zp

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine elmgt1(xx,yy,zz,xp,yp,zp,ecc,rm0,delt,dpr,
     &                  elct,elcx,elcy,elcz,
     &                  bmgt,bmgx,bmgy,bmgz)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( rlit = 29.97925d0 )

*-----------------------------------------------------------------------

               pa0 = sqrt( ecc * ( ecc + 2.d0 * rm0 ) )

               rx0 = xx
               ry0 = yy
               rz0 = zz

               px0 = xp * pa0
               py0 = yp * pa0
               pz0 = zp * pa0

               et0 = sqrt( px0**2 + py0**2 + pz0**2 + rm0**2 )

               bx0 = px0 / et0
               by0 = py0 / et0
               bz0 = pz0 / et0

*-----------------------------------------------------------------------

               nt = 10

               bt0 = pa0 / et0
               deltt = delt / bt0 / rlit
               dt = deltt / dble( nt )

*-----------------------------------------------------------------------

  100    continue

*-----------------------------------------------------------------------
*           RKG Second order
*-----------------------------------------------------------------------

               rx1 = rx0 + 0.5 * dt * bx0 * rlit
               ry1 = ry0 + 0.5 * dt * by0 * rlit
               rz1 = rz0 + 0.5 * dt * bz0 * rlit

               fx0 = elct * elcx
     &             + bmgt * ( + bmgz * by0 - bmgy * bz0 )
               fy0 = elct * elcy
     &             + bmgt * ( - bmgz * bx0 + bmgx * bz0 )
               fz0 = elct * elcz
     &             + bmgt * ( + bmgy * bx0 - bmgx * by0 )

               px1 = px0 + 0.5 * dt * fx0
               py1 = py0 + 0.5 * dt * fy0
               pz1 = pz0 + 0.5 * dt * fz0

               et1 = sqrt( px1**2 + py1**2 + pz1**2 + rm0**2 )
               bx1 = px1 / et1
               by1 = py1 / et1
               bz1 = pz1 / et1

               rxf = rx0 + dt * bx1 * rlit
               ryf = ry0 + dt * by1 * rlit
               rzf = rz0 + dt * bz1 * rlit

               fx1 = elct * elcx
     &             + bmgt * ( + bmgz * by1 - bmgy * bz1 )
               fy1 = elct * elcy
     &             + bmgt * ( - bmgz * bx1 + bmgx * bz1 )
               fz1 = elct * elcz
     &             + bmgt * ( + bmgy * bx1 - bmgx * by1 )

               pxf = px0 + dt * fx1
               pyf = py0 + dt * fy1
               pzf = pz0 + dt * fz1

*-----------------------------------------------------------------------

               dpr = dsqrt( ( rxf - xx )**2
     &                    + ( ryf - yy )**2
     &                    + ( rzf - zz )**2 )

            if( abs( dpr - delt ) / delt .lt. 0.001 ) goto 200

            if( dpr .gt. delt ) then

               dt = dt / 2.0

            else if( dpr .lt. delt ) then

               rx0 = rxf
               ry0 = ryf
               rz0 = rzf

               px0 = pxf
               py0 = pyf
               pz0 = pzf

               et0 = sqrt( px0**2 + py0**2 + pz0**2 + rm0**2 )

               bx0 = px0 / et0
               by0 = py0 / et0
               bz0 = pz0 / et0

            end if

               goto 100

*-----------------------------------------------------------------------

  200    continue

               xx = rxf
               yy = ryf
               zz = rzf

               pf2 = pxf**2 + pyf**2 + pzf**2
               pf0 = sqrt( pf2 )
               et0 = sqrt( pf2 + rm0**2 )

               ecc = et0 - rm0

               xp = pxf / pf0
               yp = pyf / pf0
               zp = pzf / pf0

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine magtrs(a_mag,b_mag,s_mag,p_mag,t_mag,
     &                  delt,dpr,udir,vdir,wdir)
*                                                                      *
*                                                                      *
*       particle transfer in void under magnetic field                 *
*       created  by S.Meigo on 2000/07/01                              *
*       last modified by K.Niita on 2010/12/23                         *
*                                                                      *
*     input  :                                                         *
*       delt   : distance                                              *
*       a_mag  : magnet gap(mm)                                        *
*       b_mag  : magnet field at pole tip  [kG]                        *
*       s_mag  : speicies of magnet dypole:2, quad:4, sext:6, oct:8    *
*       p_mag  : phase of magnet                                       *
*       t_mag  : transform id                                          *
*                                                                      *
*     in common                                                        *
*                                                                      *
*       e(no)    : initial energy                                      *
*       x(no),y(no),z(no)                                              *
*                : initial position                                    *
*       u(no),v(no),w(no)                                              *
*                : initial direction                                   *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       ec(no)   : final energy                                        *
*       xc(no),yc(no),zc(no)                                           *
*                : final position                                      *
*                                                                      *
*       u(no),v(no),w(no)                                              *
*                : initial direction                                   *
*                                                                      *
*       udir,vdir,wdir                                                 *
*                : final direction                                     *
*                                                                      *
*       dpr      : final distance                                      *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
cFURUTA      include 'mmbank.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /gravit/ grav(3), igrav

      parameter ( alph = 5.7688252d0 )
      parameter ( gamn = 1.8324712d+8 )
      parameter ( grvt = 980.665d0 )
      parameter ( cvel = 2.997925d+10 )
      parameter ( pi   = 3.1415926535898d0 )
      parameter ( rlit = 29.97925d0 )
      parameter ( grvc = 980.665d-18 )

*-----------------------------------------------------------------------
*     transform
*-----------------------------------------------------------------------

            itrs = nint( t_mag )

            call trnsxx(x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                  xx0,yy0,zz0,itrs)

            call trnsuu(u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                  uu0,vv0,ww0,itrs)

*-----------------------------------------------------------------------
*     gravity
*-----------------------------------------------------------------------

               jgrav = 0

         if( igrav .ne. 0 .and. ityp .eq. 2 .and.
     &       e(ibke+no) .gt. 0.0 .and. e(ibke+no) .le. 1.d-6 ) then

               call trnsuu(grav(1),grav(2),grav(3),
     &                     gravx,gravy,gravz,itrs)

               jgrav = 1

         end if

*-----------------------------------------------------------------------
*     zero magnet field
*-----------------------------------------------------------------------

            if( b_mag .eq. 0.0d0  ) goto 1000

*-----------------------------------------------------------------------
*     which magnetic field and method
*-----------------------------------------------------------------------

               isp = nint( s_mag )

            if( isp .eq. 60 ) then

               isp = 6
               is6 = 0

            else if( isp .eq. 61 ) then

               isp = 6
               is6 = 1

            else if( isp .eq. 62 ) then

               isp = 100
               is6 = 6

            else if( isp .eq. 101 ) then

               isp = 100
               is6 = 1

            else if( isp .eq. 102 ) then

               isp = 100
               is6 = 2

            else if( isp .eq. 103 ) then

               isp = 100
               is6 = 3

            else if( isp .eq. 104 ) then

               isp = 100
               is6 = 4

            else if( isp .eq. 106 ) then

               isp = 100
               is6 = 6

            end if

*-----------------------------------------------------------------------
*     drift
*-----------------------------------------------------------------------

            if( ( isp .eq. 4 .or. isp .eq. 6 ) .and.
     &            abs(ww0) .lt. 1.0d-5 ) goto 1000

            if( ( isp .eq. 2 .or. isp .eq. 3 ) .and.
     &          uu0**2 + ww0**2 .lt. 1.0d-5 ) goto 1000

!	COMET magnetic field map reading option: 12

*-----------------------------------------------------------------------
*      magnet field
*-----------------------------------------------------------------------
            if( isp .eq. 2 ) then

                  b = b_mag
                  p = dsqrt(e(ibke+no)**2+2.*rtyp*e(ibke+no))/1000.		! momentum [GeV/c]

                  jchg = nzst(ibkzst+no)
                  cg = dble(jchg)      									! charge
cKN 2011/12/05    cg = dble(jtyp)
                  r0 = 33.356 * p / abs( cg * b ) * 100.d0				! radius [cm]

                  aw0 = dsqrt( uu0**2 + ww0**2 )						! direction of particle momentum

                  icm = -1
                  if( cg * b .lt. 0.0d0 ) icm = 1

               if( delt .gt. 2.0 * r0 / aw0 ) then

                  delt = 2.0 * r0 / aw0

               end if

               if( abs(ww0) .gt. 0.0d0 ) then

                  the2 = atan( uu0 / ww0 )								! angle between x and z

                  if( ww0 .lt. 0.0 ) the2 = the2 + pi

               else

                  the2 = pi / 2.0d0

                  if( uu0 .lt. 0.0 ) the2 = the2 + pi

               end if
  
                  cst = cos( the2 )
                  snt = sin( the2 )

*-----------------------------------------------------------------------

            else if( isp .eq. 3 ) then

				  call usrmgf3(x(ibkx+no),y(ibky+no),z(ibkz+no),bxy,btot,
     &                   bbx,bby,bbz)

				  angx = (pi/2-acos(bbz/btot))*180/pi
				  angy = 0.0
				  angz = (-pi/2+atan2(bby,bbx))*180/pi

				  rotx = pi/2-acos(bbz/btot)									! angle of around X (degree)
				  roty = 0.0													! angle of around Y (degree)
				  rotz = pi/2-atan2(bby,bbx)									! angle of around Z (degree)
				  
!				  write (*,*) "rotx, roty, rotz", angx, angy, angz

!		transform the coordinate

      			  cmgm11 = cos(rotz)*cos(roty)
			      cmgm12 = sin(rotz)*cos(rotx)+cos(rotz)*sin(roty)*sin(rotx)
      			  cmgm13 = sin(rotz)*sin(rotx)-cos(rotz)*sin(roty)*cos(rotx)
     			  cmgm21 = -sin(rotz)*cos(roty)
				  cmgm22 = cos(rotz)*cos(rotx)-sin(rotz)*sin(roty)*sin(rotx)
     			  cmgm23 = cos(rotz)*sin(rotx)+sin(rotz)*sin(roty)*cos(rotx)
   				  cmgm31 = sin(roty)
     			  cmgm32 = -cos(roty)*sin(rotx)
     			  cmgm33 = cos(roty)*cos(rotx)
				  
				  xx0 = cmgm11*x(ibkx+no)+cmgm21*y(ibky+no)+cmgm31*z(ibkz+no)
				  yy0 = cmgm12*x(ibkx+no)+cmgm22*y(ibky+no)+cmgm32*z(ibkz+no)
				  zz0 = cmgm13*x(ibkx+no)+cmgm23*y(ibky+no)+cmgm33*z(ibkz+no)
				  
				  uu0 = cmgm11*u(ibku+no)+cmgm21*v(ibkv+no)+cmgm31*w(ibkw+no)
				  vv0 = cmgm12*u(ibku+no)+cmgm22*v(ibkv+no)+cmgm32*w(ibkw+no)
				  ww0 = cmgm13*u(ibku+no)+cmgm23*v(ibkv+no)+cmgm33*w(ibkw+no)
				  
				  b = btot * 10											! [T] -> [kG]
                  p = dsqrt(e(ibke+no)**2+2.*rtyp*e(ibke+no))/1000.		! momentum [GeV/c]

!				  write (*,*) "position : ", x(ibkx+no),y(ibky+no),z(ibkz+no)
!				  write (*,*) "magnetic field : ", bbx, bby, bbz, btot

                  jchg = nzst(ibkzst+no)
                  cg = dble(jchg)      									! charge
cKN 2011/12/05    cg = dble(jtyp)
                  r0 = 33.356 * p / abs( cg * b ) * 100.d0				! radius [cm]

                  aw0 = dsqrt( uu0**2 + ww0**2 )						! direction of particle momentum

                  icm = -1
                  if( cg * b .lt. 0.0d0 ) icm = 1

               if( delt .gt. 2.0 * r0 / aw0 ) then

                  delt = 2.0 * r0 / aw0

               end if

               if( abs(ww0) .gt. 0.0d0 ) then

                  the2 = atan( uu0 / ww0 )								! angle between x and z

                  if( ww0 .lt. 0.0 ) the2 = the2 + pi

               else

                  the2 = pi / 2.0d0

                  if( uu0 .lt. 0.0 ) the2 = the2 + pi

               end if
				  
                  cst = cos( the2 )
                  snt = sin( the2 )

*-----------------------------------------------------------------------

            else if( isp .eq. 4 ) then

                  a = a_mag
                  b = b_mag
                  p = dsqrt(e(ibke+no)**2+2.*rtyp*e(ibke+no))/1000.

                  jchg = nzst(ibkzst+no)
                  cg = dble(jchg)
cKN 2011/12/05    cg = dble(jtyp)

                  uu = uu0 / ww0
                  vv = vv0 / ww0

                  xp0 = atan(uu)
                  yp0 = atan(vv)

*-----------------------------------------------------------------------

            else if( isp .eq. 6 ) then

                  omg = sqrt( alph * abs( b_mag ) )
                  sfc = gamn * omg / alph / 10000.0
                  grv = grvt / omg**2
                  vel = sqrt( 2. * e(ibke+no) / rtyp ) * cvel

                  xp0 = vel / omg * uu0
                  yp0 = vel / omg * vv0
                  zp0 = vel / omg * ww0

*-----------------------------------------------------------------------

            else if( isp .eq. 100 ) then

                  omg = sqrt( alph * abs( b_mag ) )
                  sfc = gamn * omg / alph / 10000.0
                  grv = grvt / omg**2
                  vel = sqrt( 2. * e(ibke+no) / rtyp ) * cvel

                  cmg = abs( b_mag ) / 10000.d0
                  dmg = a_mag / cmg

                  xp0 = vel / omg * uu0
                  yp0 = vel / omg * vv0
                  zp0 = vel / omg * ww0

            end if

*-----------------------------------------------------------------------

         deltv = delt

         icc = 0
 5000    icc = icc + 1

*-----------------------------------------------------------------------

            if( isp .eq. 2 ) then

                  dll = deltv * aw0
                  al1 = dll / r0

                  dpr = sqrt( ( 2.0 * r0 * sin(al1/2.0) )**2
     &                      + ( vv0 / aw0 * dll )**2 )

*-----------------------------------------------------------------------
!	COMET magnetic field option:

            else if( isp .eq. 3 ) then

                  dll = deltv * aw0
                  al1 = dll / r0

                  dpr = sqrt( ( 2.0 * r0 * sin(al1/2.0) )**2
     &                      + ( vv0 / aw0 * dll )**2 )

*-----------------------------------------------------------------------

            else if( isp .eq. 4 ) then

                  xx = xx0
                  yy = yy0
                  xp = xp0
                  yp = yp0
                  dz = deltv * ww0

                  call quad(p,b,a,dz,xx,xp,yy,yp,cg)

                  xxc = xx
                  yyc = yy
                  zzc = zz0 + dz

                  dpr = dsqrt( ( xx0 - xxc )**2
     &                       + ( yy0 - yyc )**2
     &                       + ( zz0 - zzc )**2 )

*-----------------------------------------------------------------------

            else if( isp .eq. 6 ) then

                     xx = xx0
                     xp = xp0
                     yy = yy0
                     yp = yp0
                     zz = zz0
                     zp = zp0

                     dz = deltv * ww0
                     thet = dz / zp0

*-----------------------------------------------------------------------
*                 initial spin
*-----------------------------------------------------------------------

                     ssr = spx(ibkspx+no)**2
     &                   + spy(ibkspy+no)**2
     &                   + spz(ibkspz+no)**2

                  if( ssr .lt. 1.d-8 ) then

                        pol = p_mag
                        sn = 0.d0

                     if( pol .lt. -1.d0 ) then
                        pol  =  0.d0
                     else if( pol .gt. 1.d0 ) then
                        pol  =  100.d0
                        ipal = 0
                     end if

                     if( pol .ge. -1.d0 .and. pol .le. 1.d0 ) then

                        pol = ( pol + 1.d0 ) / 2.d0
                        ipal = 1
                        if( unirn(dummy) .gt. pol ) ipal = -1

                        sn = ( xx**2 + yy**2 ) / 2.d0

                     end if

                     if( sn .gt. 0.0d0 .and. ipal .ne. 0 ) then

                        sx = ipal * ( yy**2 - xx**2 ) / 2.d0 / sn
                        sy = ipal *   xx * yy  / sn
                        sz = 0.0d0

                     else

                        th = 2.0d0 * pi * unirn(dummy)
                        cs = 2.d0 *  unirn(dummy) - 1.d0
                        sn = sqrt( 1.d0 - cs**2 )
                        sx = sn * cos( th )
                        sy = sn * sin( th )
                        sz = cs

                     end if

                  else

                       call trnsuu(spx(ibkspx+no),
     &                             spy(ibkspy+no),
     &                             spz(ibkspz+no),sx,sy,sz,itrs)

                  end if

*-----------------------------------------------------------------------

               if( is6 .eq. 0 ) then

                     call sexts0(dz,xx,xp,yy,yp,sx,sy,sz,thet)

                     xxc = xx
                     yyc = yy
                     zzc = zz + dz

                     dpr = dsqrt( ( xx0 - xxc )**2
     &                          + ( yy0 - yyc )**2
     &                          + ( zz0 - zzc )**2 )

               else if( is6 .eq. 1 ) then

                  call sexts1(xx,xp,yy,yp,zz,zp,sx,sy,sz,
     &                        thet,delt,dpr,sfc,
     &                        grv,gravx,gravy,gravz)

                     xxc = xx
                     yyc = yy
                     zzc = zz

               end if

*-----------------------------------------------------------------------
*        special version
*-----------------------------------------------------------------------

            else if( isp .eq. 100 ) then

                     xx = xx0
                     xp = xp0
                     yy = yy0
                     yp = yp0
                     zz = zz0
                     zp = zp0

                     dz = deltv * ww0
                     thet = dz / zp0

*-----------------------------------------------------------------------
*                 initial spin
*-----------------------------------------------------------------------

                     ssr = spx(ibkspx+no)**2
     &                   + spy(ibkspy+no)**2
     &                   + spz(ibkspz+no)**2

                  if( ssr .lt. 1.d-8 ) then

                        pol = p_mag
                        bba = 0.d0

                     if( pol .lt. -1.d0 ) then
                        pol  =  0.d0
                     else if( pol .gt. 1.d0 ) then
                        pol  =  100.d0
                        ipal = 0
                     end if

                     if( pol .ge. -1.d0 .and. pol .le. 1.d0 ) then

                        pol = ( pol + 1.d0 ) / 2.d0
                        ipal = 1
                        if( unirn(dummy) .gt. pol ) ipal = -1

                        call magbdx(is6,xx,yy,zz,dmg,cmg,
     &                              bbx,bby,bbz,bba,
     &                              dxx,dyx,dzx,
     &                              dxy,dyy,dzy,
     &                              dxz,dyz,dzz)

                     end if

                     if( bba .gt. 0.0d0 .and. ipal .ne. 0 ) then

                        sx = ipal * bbx / bba
                        sy = ipal * bby / bba
                        sz = ipal * bbz / bba

                     else

                        th = 2.0d0 * pi * unirn(dummy)
                        cs = 2.d0 *  unirn(dummy) - 1.d0
                        sn = sqrt( 1.d0 - cs**2 )
                        sx = sn * cos( th )
                        sy = sn * sin( th )
                        sz = cs

                     end if

                  else

                        call trnsuu(spx(ibkspx+no),
     &                              spy(ibkspy+no),
     &                              spz(ibkspz+no),sx,sy,sz,itrs)

                  end if

*-----------------------------------------------------------------------

                  call sexts2(xx,xp,yy,yp,zz,zp,sx,sy,sz,
     &                        thet,delt,dpr,sfc,dmg,cmg,
     &                        grv,gravx,gravy,gravz,is6)

                     xxc = xx
                     yyc = yy
                     zzc = zz

            end if

*-----------------------------------------------------------------------
*           adjust the distanc to delt
*-----------------------------------------------------------------------

            if( icc .le. 100 .and.
     &          abs( dpr - delt ) / delt .gt. 0.001 ) then

               factv = ( delt - dpr ) / delt
               factv = sign( min( 1.0d0, abs( factv ) ), factv )

               deltv = deltv + deltv * factv * 0.9

               goto 5000

            end if

*-----------------------------------------------------------------------

            if( isp .eq. 2 ) then

               cs = cos( al1 )
               sn = sin( al1 )
               zz = r0 * sn
               xx = ( r0 - r0 * cs ) * icm

               zzc = zz0 + zz * cst - xx * snt			! position
               xxc = xx0 + zz * snt + xx * cst
               yyc = yy0 + dll * vv0 / aw0

               zz = cs
               xx = sn * icm

               vz = ( zz * cst - xx * snt ) * aw0
               vx = ( zz * snt + xx * cst ) * aw0
               vy = vv0

               av = dsqrt( vx**2 + vy**2 + vz**2 )

               uur = vx / av		! direction
               vvr = vy / av
               wwr = vz / av
               
           else if( isp .eq. 3 ) then

               cs = cos( al1 )
               sn = sin( al1 )
               zz = r0 * sn
               xx = ( r0 - r0 * cs ) * icm

               zzc = zz0 + zz * cst - xx * snt			! final position
               xxc = xx0 + zz * snt + xx * cst
               yyc = yy0 + dll * vv0 / aw0

               zz = cs
               xx = sn * icm

               vz = ( zz * cst - xx * snt ) * aw0
               vx = ( zz * snt + xx * cst ) * aw0
               vy = vv0

               av = dsqrt( vx**2 + vy**2 + vz**2 )

               uur = vx / av		! direction
               vvr = vy / av
               wwr = vz / av

            else if( isp .eq. 4 ) then

               vx = ww0 * tan( xp )
               vy = ww0 * tan( yp )
               vz = ww0 * 1.0d0

               av = dsqrt( vx**2 + vy**2 + vz**2 )

               uur = vx / av
               vvr = vy / av
               wwr = vz / av

            else if( isp .eq. 6 ) then

               av = dsqrt( xp**2 + yp**2 + zp**2 )
               vl = av * omg

               ec(ibkec+no) = 0.5 * rtyp * ( vl / cvel )**2

               uur = xp / av
               vvr = yp / av
               wwr = zp / av

            else if( isp .eq. 100 ) then

               av = dsqrt( xp**2 + yp**2 + zp**2 )
               vl = av * omg

               ec(ibkec+no) = 0.5 * rtyp * ( vl / cvel )**2

               uur = xp / av
               vvr = yp / av
               wwr = zp / av

               sv = dsqrt( sx**2 + sy**2 + sz**2 )

               sx = sx / sv
               sy = sy / sv
               sz = sz / sv

            end if

               uuc = ( xxc - xx0 ) / dpr
               vvc = ( yyc - yy0 ) / dpr
               wwc = ( zzc - zz0 ) / dpr

*-----------------------------------------------------------------------
*        transform
*-----------------------------------------------------------------------

			if ( isp .ne. 3 ) then

            call trnsxv(xxc,yyc,zzc,
     &                  xc(ibkxc+no),yc(ibkyc+no),zc(ibkzc+no),itrs)

            call trnsuv(uuc,vvc,wwc,
     &                  u(ibku+no),v(ibkv+no),w(ibkw+no),itrs)

            call trnsuv(uur,vvr,wwr,
     &                  udir,vdir,wdir,itrs)

            call trnsuv(sx,sy,sz,
     &           spx(ibkspx+no),spy(ibkspy+no),spz(ibkspz+no),itrs)
     		
     		end if
     
     		! return to the original coordinate
     		if( isp .eq. 3 ) then
     		
     		xc(ibkxc+no) = cmgm11*xxc+cmgm12*yyc+cmgm13*zzc
     		yc(ibkyc+no) = cmgm21*xxc+cmgm22*yyc+cmgm23*zzc
     		zc(ibkzc+no) = cmgm31*xxc+cmgm32*yyc+cmgm33*zzc
     		
     		u(ibku+no) = cmgm11*uuc+cmgm12*vvc+cmgm13*wwc
     		v(ibkv+no) = cmgm21*uuc+cmgm22*vvc+cmgm23*wwc
     		w(ibkw+no) = cmgm31*uuc+cmgm32*vvc+cmgm33*wwc
     		
     		udir = cmgm11*uur+cmgm12*vvr+cmgm13*wwr
     		vdir = cmgm21*uur+cmgm22*vvr+cmgm23*wwr
     		wdir = cmgm31*uur+cmgm32*vvr+cmgm33*wwr
     		
     		spx(ibkspx+no) = cmgm11*sx+cmgm12*sy+cmgm13*sz
     		spy(ibkspy+no) = cmgm21*sx+cmgm22*sy+cmgm23*sz
     		spz(ibkspz+no) = cmgm31*sx+cmgm32*sy+cmgm33*sz
     		
     		! debug
!     		xc(ibkxc+no) = xxc
!     		yc(ibkyc+no) = yyc
!     		zc(ibkzc+no) = zzc
     		
!     		u(ibku+no) = uuc
!     		v(ibkv+no) = vvc
!     		w(ibkv+no) = wwc
     		
!     		udir = uur
!     		vdir = vvr
!     		wdir = wwr
     		
!     		spx(ibkspx+no) = sx
!     		spy(ibkspy+no) = sy
!     		spz(ibkspz+no) = sz
     		
     		end if

            return

*-----------------------------------------------------------------------

 1000 continue

               xc(ibkxc+no) = x(ibkx+no) + delt * u(ibku+no)
               yc(ibkyc+no) = y(ibky+no) + delt * v(ibkv+no)
               zc(ibkzc+no) = z(ibkz+no) + delt * w(ibkw+no)

            if( jgrav .ne. 0 ) then

               ekin = e(ibke+no)
               vel  = sqrt( 2.0 * ekin / rtyp ) * rlit
               timd = delt / vel

               vxx = u(ibku+no) * vel - grvc * grav(1) * timd
               vyy = v(ibkv+no) * vel - grvc * grav(2) * timd
               vzz = w(ibkw+no) * vel - grvc * grav(3) * timd

               vab = sqrt( vxx**2 + vyy**2 + vzz**2 )

               ec(ibkec+no) = 0.5 * rtyp * ( vab / rlit )**2

               xc(ibkxc+no) = xc(ibkxc+no)
     &                      - 0.5 * grvc * grav(1) * timd**2
               yc(ibkyc+no) = yc(ibkyc+no)
     &                      - 0.5 * grvc * grav(2) * timd**2
               zc(ibkzc+no) = zc(ibkzc+no)
     &                      - 0.5 * grvc * grav(3) * timd**2

               u(ibku+no)   = vxx / vab
               v(ibkv+no)   = vyy / vab
               w(ibkw+no)   = vzz / vab

               delt = sqrt( ( xc(ibkxc+no) - x(ibkx+no) )**2
     &                    + ( yc(ibkyc+no) - y(ibky+no) )**2
     &                    + ( zc(ibkzc+no) - z(ibkz+no) )**2 )

            end if

*-----------------------------------------------------------------------

               udir = u(ibku+no)
               vdir = v(ibkv+no)
               wdir = w(ibkw+no)

               dpr = delt

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine magbdx(is6,xx,yy,zz,dmg,cmg,
     &                  bbx,bby,bbz,bba,
     &                  dxx,dyx,dzx,
     &                  dxy,dyy,dzy,
     &                  dxz,dyz,dzz)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
*     diple only
*-----------------------------------------------------------------------

      if( is6 .eq. 2 ) then

            bbx =  1.d0
            bby =  0.d0
            bbz =  0.d0

            bba = sqrt( bbx**2 + bby**2 + bbz**2 )

            dxx =  0.d0
            dyx =  0.d0
            dzx =  0.d0

            dxy =  0.d0
            dyy =  0.d0
            dzy =  0.d0

            dxz =  0.d0
            dyz =  0.d0
            dzz =  0.d0
            

*-----------------------------------------------------------------------
*     quadrapole and dipole
*-----------------------------------------------------------------------

      else if( is6 .eq. 4 ) then

            bbx =  xx
            bby = -yy
            bbz = dmg

            bba = sqrt( bbx**2 + bby**2 + bbz**2 )

            dxx =  1.d0
            dyx =  0.d0
            dzx =  0.d0

            dxy =  0.d0
            dyy = -1.d0
            dzy =  0.d0

            dxz =  0.d0
            dyz =  0.d0
            dzz =  0.d0
            

*-----------------------------------------------------------------------
*     sextapole and dipole
*-----------------------------------------------------------------------

      else if( is6 .eq. 6 ) then

            bbx = ( yy**2 - xx**2 ) / 2.d0
            bby = xx * yy
            bbz = dmg

            bba = sqrt( bbx**2 + bby**2 + bbz**2 )

            dxx = - xx
            dyx =   yy
            dzx = 0.d0

            dxy =   yy
            dyy =   xx
            dzy = 0.d0

            dxz = 0.d0
            dyz = 0.d0
            dzz = 0.d0

*-----------------------------------------------------------------------
*     from user subroutine 1
*-----------------------------------------------------------------------

      else if( is6 .eq. 1 ) then

!			write (*,*) "position:", xx, yy, zz

            call usrmgf1(xx,yy,zz,dmg,cmg,
     &                   bbx,bby,bbz,
     &                   dxx,dyx,dzx,
     &                   dxy,dyy,dzy,
     &                   dxz,dyz,dzz)

            bba = sqrt( bbx**2 + bby**2 + bbz**2 )
!            write (*,*) "position:", xx, yy, zz
!            write (*,*) "magnetic field:",bbx, bby, bbz, bba

*-----------------------------------------------------------------------
*     from user subroutine 3
*-----------------------------------------------------------------------

      else if( is6 .eq. 3 ) then
      
!      		write (*,*) "x, y, z: ", xx, yy, zz

            call usrmgf3(xx,yy,zz,dmg,cmg,
     &                   bbx,bby,bbz)

            bba = sqrt( bbx**2 + bby**2 + bbz**2 )
!            write (*,*) "total magnetic field:", bba
            
*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sexts2(xx,xp,yy,yp,zz,zp,sx,sy,sz,
     &                  thet,delt,dpr,sfc,dmg,cmg,
     &                  grv,gravx,gravy,gravz,is6)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi   = 3.1415926535898d0 )

*-----------------------------------------------------------------------

            xx0 = xx
            vx0 = xp
            yy0 = yy
            vy0 = yp
            zz0 = zz
            vz0 = zp

            sx0 = sx
            sy0 = sy
            sz0 = sz

            sxf = sx
            syf = sy
            szf = sz

*-----------------------------------------------------------------------

            nt  = 10
            dt0 = thet / dble( nt )
            ispm = 0
            dt = dt0
            idt = 0

*-----------------------------------------------------------------------

  100    continue

               call magbdx(is6,xx0,yy0,zz0,dmg,cmg,
     &                     bbx0,bby0,bbz0,bba0,
     &                     dxx0,dyx0,dzx0,
     &                     dxy0,dyy0,dzy0,
     &                     dxz0,dyz0,dzz0)

*-----------------------------------------------------------------------

               sn = bba0
               dts = dt0
               if( idt .eq. 1 ) dts = dt

*-----------------------------------------------------------------------

            if( sn * sfc .lt. 1000.0 ) then

               ispm = 1

               dt = dts
               if( sn * sfc .gt. 0.d0 )
     &         dt = min( dts, 0.05 / ( sfc * sn ) )

            else

               ispm = 0
               dt = dts

            end if

*-----------------------------------------------------------------------
*           RKG Second order
*-----------------------------------------------------------------------

         if( ispm .eq. 0 ) then

*-----------------------------------------------------------------------
*           spin rotation
*-----------------------------------------------------------------------

               sba = ( sx0 * bbx0 + sy0 * bby0 + sz0 * bbz0 ) / bba0
               sba = max( -1.d0, sba )
               sba = min(  1.d0, sba )

               costh = bbz0 / bba0
               rt2 = ( bbx0 / bba0 )**2 + ( bby0 / bba0 )**2

            if( rt2 .eq. 0.0d0 ) then

               sinth = 0.0d0
               cosphi= 1.0d0
               sinphi= 0.0d0

            else

               rt = sqrt(rt2)
               sinth  = rt
               cosphi =   ( bbx0 / bba0 ) / rt
               sinphi = - ( bby0 / bba0 ) / rt

            end if

               adc = sx0
               bdc = sy0
               gdc = sz0

               t1  = costh  * adc + sinth  * gdc
               sxb = cosphi * t1  - sinphi * bdc
               syb = sinphi * t1  + cosphi * bdc
               szb = costh  * gdc - sinth  * adc

            if( syb .ge. 0.d0 .and. sxb .eq. 0.d0 ) then
               sbt = pi / 2.d0
            else if( syb .lt. 0.d0 .and. sxb .eq. 0.d0 ) then
               sbt = pi / 2.d0 * 3.d0
            else if( syb .ge. 0.d0 .and. sxb .gt. 0.d0 ) then
               sbt = atan( syb / sxb )
            else if( syb .ge. 0.d0 .and. sxb .lt. 0.d0 ) then
               sbt = atan( syb / sxb ) + pi
            else if( syb .lt. 0.d0 .and. sxb .gt. 0.d0 ) then
               sbt = atan( syb / sxb ) + 2.d0 * pi
            else if( syb .lt. 0.d0 .and. sxb .lt. 0.d0 ) then
               sbt = atan( syb / sxb ) + pi
            end if

*-----------------------------------------------------------------------

               xx1 = xx0 + 0.5 * dt * vx0
               yy1 = yy0 + 0.5 * dt * vy0
               zz1 = zz0 + 0.5 * dt * vz0

               fx0 = - (  bbx0 / bba0 * dxx0
     &                 +  bby0 / bba0 * dyx0
     &                 +  bbz0 / bba0 * dzx0 ) * sba
     &               - grv * gravx

               fy0 = - (  bbx0 / bba0 * dxy0
     &                 +  bby0 / bba0 * dyy0
     &                 +  bbz0 / bba0 * dzy0 ) * sba
     &               - grv * gravy

               fz0 = - (  bbx0 / bba0 * dxz0
     &                 +  bby0 / bba0 * dyz0
     &                 +  bbz0 / bba0 * dzz0 ) * sba
     &               - grv * gravz

               vx1 = vx0 + 0.5 * dt * fx0
               vy1 = vy0 + 0.5 * dt * fy0
               vz1 = vz0 + 0.5 * dt * fz0

               xxf = xx0 + dt * vx1
               yyf = yy0 + dt * vy1
               zzf = zz0 + dt * vz1

               call magbdx(is6,xx1,yy1,zz1,dmg,cmg,
     &                     bbx1,bby1,bbz1,bba1,
     &                     dxx1,dyx1,dzx1,
     &                     dxy1,dyy1,dzy1,
     &                     dxz1,dyz1,dzz1)

               fx1 = - (  bbx1 / bba1 * dxx1
     &                 +  bby1 / bba1 * dyx1
     &                 +  bbz1 / bba1 * dzx1 ) * sba
     &               - grv * gravx

               fy1 = - (  bbx1 / bba1 * dxy1
     &                 +  bby1 / bba1 * dyy1
     &                 +  bbz1 / bba1 * dzy1 ) * sba
     &               - grv * gravy

               fz1 = - (  bbx1 / bba1 * dxz1
     &                 +  bby1 / bba1 * dyz1
     &                 +  bbz1 / bba1 * dzz1 ) * sba
     &               - grv * gravz

               vxf = vx0 + dt * fx1
               vyf = vy0 + dt * fy1
               vzf = vz0 + dt * fz1

               call magbdx(is6,xxf,yyf,zzf,dmg,cmg,
     &                     bbxf,bbyf,bbzf,bbaf,
     &                     dxxf,dyxf,dzxf,
     &                     dxyf,dyyf,dzyf,
     &                     dxzf,dyzf,dzzf)

*-----------------------------------------------------------------------
*           spin rotation
*-----------------------------------------------------------------------

               costh = bbzf / bbaf
               rt2 = ( bbxf / bbaf )**2 + ( bbyf / bbaf )**2

            if( rt2 .eq. 0.0d0 ) then

               sinth = 0.0d0
               cosphi= 1.0d0
               sinphi= 0.0d0

            else

               rt = sqrt(rt2)
               sinth  = rt
               cosphi = ( bbxf / bbaf ) / rt
               sinphi = ( bbyf / bbaf ) / rt

            end if

               cosh  = sba
               sinh  = sqrt( 1.d0 - cosh**2 )

               thet = sbt + bba1 * sfc * dt

               cost = cos( thet )
               sint = sin( thet )

               adc = sinh * cost
               bdc = sinh * sint
               gdc = cosh

               t1  = costh  * adc + sinth  * gdc
               sxf = cosphi * t1  - sinphi * bdc
               syf = sinphi * t1  + cosphi * bdc
               szf = costh  * gdc - sinth  * adc

*-----------------------------------------------------------------------

         else

*-----------------------------------------------------------------------

               xx1 = xx0 + 0.5 * dt * vx0
               yy1 = yy0 + 0.5 * dt * vy0
               zz1 = zz0 + 0.5 * dt * vz0

               fx0 = - (  sx0 * dxx0
     &                 +  sy0 * dyx0
     &                 +  sz0 * dzx0 )
     &               - grv * gravx

               fy0 = - (  sx0 * dxy0
     &                 +  sy0 * dyy0
     &                 +  sz0 * dzy0 )
     &               - grv * gravy

               fz0 = - (  sx0 * dxz0
     &                 +  sy0 * dyz0
     &                 +  sz0 * dzz0 )
     &               - grv * gravz

               vx1 = vx0 + 0.5 * dt * fx0
               vy1 = vy0 + 0.5 * dt * fy0
               vz1 = vz0 + 0.5 * dt * fz0

               gx0 = sfc * ( sy0 * bbz0 - sz0 * bby0 )
               gy0 = sfc * ( sz0 * bbx0 - sx0 * bbz0 )
               gz0 = sfc * ( sx0 * bby0 - sy0 * bbx0 )

               sx1 = sx0 + 0.5 * dt * gx0
               sy1 = sy0 + 0.5 * dt * gy0
               sz1 = sz0 + 0.5 * dt * gz0

               xxf = xx0 + dt * vx1
               yyf = yy0 + dt * vy1
               zzf = zz0 + dt * vz1

               call magbdx(is6,xx1,yy1,zz1,dmg,cmg,
     &                     bbx1,bby1,bbz1,bba1,
     &                     dxx1,dyx1,dzx1,
     &                     dxy1,dyy1,dzy1,
     &                     dxz1,dyz1,dzz1)

               fx1 = - (  sx1 * dxx1
     &                 +  sy1 * dyx1
     &                 +  sz1 * dzx1 )
     &               - grv * gravx

               fy1 = - (  sx1 * dxy1
     &                 +  sy1 * dyy1
     &                 +  sz1 * dzy1 )
     &               - grv * gravy

               fz1 = - (  sx1 * dxz1
     &                 +  sy1 * dyz1
     &                 +  sz1 * dzz1 )
     &               - grv * gravz

               vxf = vx0 + dt * fx1
               vyf = vy0 + dt * fy1
               vzf = vz0 + dt * fz1

               gx1 = sfc * ( sy1 * bbz1 - sz1 * bby1 )
               gy1 = sfc * ( sz1 * bbx1 - sx1 * bbz1 )
               gz1 = sfc * ( sx1 * bby1 - sy1 * bbx1 )

               sxf = sx0 + dt * gx1
               syf = sy0 + dt * gy1
               szf = sz0 + dt * gz1

         end if

*-----------------------------------------------------------------------

            dpr = dsqrt( ( xxf - xx )**2
     &                 + ( yyf - yy )**2
     &                 + ( zzf - zz )**2 )

         if( abs( dpr - delt ) / delt .lt. 0.001 ) goto 200

         if( dpr .gt. delt ) then

            dt = dt / 2.0
            idt = 1

         else if( dpr .lt. delt ) then

            xx0 = xxf
            yy0 = yyf
            zz0 = zzf

            vx0 = vxf
            vy0 = vyf
            vz0 = vzf

            sx0 = sxf
            sy0 = syf
            sz0 = szf

         end if

             ssr = sqrt( sx0**2 + sy0**2 + sz0**2 )

           if( ssr .gt. 1.05 .or. ssr .lt. 0.95 )  then
             write(*,*) 'ERROR: integration is diverged !!', ssr
             write(*,*) '     ispm =', ispm
             write(*,*) '   sn*sfc =', sn * sfc
             write(*,*) '       dt =', dt
             stop 777
           end if

            goto 100

*-----------------------------------------------------------------------

  200    continue

            xx = xxf
            yy = yyf
            zz = zzf

            xp = vxf
            yp = vyf
            zp = vzf

            sx = sxf
            sy = syf
            sz = szf

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sexts1(xx,xp,yy,yp,zz,zp,sx,sy,sz,
     &                  thet,delt,dpr,sfc,
     &                  grv,gravx,gravy,gravz)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi   = 3.1415926535898d0 )

*-----------------------------------------------------------------------

            xx0 = xx
            vx0 = xp
            yy0 = yy
            vy0 = yp
            zz0 = zz
            vz0 = zp

            sx0 = sx
            sy0 = sy
            sz0 = sz

            sxf = sx
            syf = sy
            szf = sz

*-----------------------------------------------------------------------

            nt  = 10
            dt0 = thet / dble( nt )
            ispm = 0
            dt = dt0
            idt = 0

*-----------------------------------------------------------------------

  100    continue

               sn = ( xx0**2 + yy0**2 ) / 2.d0
               dts = dt0
               if( idt .eq. 1 ) dts = dt

*-----------------------------------------------------------------------

            if( sn * sfc .lt. 1000.0 ) then

               ispm = 1
               dt = dts
               if( sn * sfc .gt. 0.d0 )
     &         dt = min( dts, 0.05 / ( sfc * sn ) )

            else

               ispm = 0
               dt = dts

            end if

*-----------------------------------------------------------------------
*           RKG Second order
*-----------------------------------------------------------------------

         if( ispm .eq. 0 ) then

*-----------------------------------------------------------------------
*           spin rotation
*-----------------------------------------------------------------------

               sba = ( ( yy0**2 - xx0**2 ) / 2.d0 * sx0
     &                 + xx0 * yy0 * sy0 ) / sn
               sba = max( -1.d0, sba )
               sba = min(  1.d0, sba )

               costh  = 0.d0
               sinth  = 1.d0
               cosphi =   ( yy0**2 - xx0**2 ) / 2.d0 / sn
               sinphi = - xx0 * yy0 / sn

               adc = sx0
               bdc = sy0
               gdc = sz0

               t1  = costh  * adc + sinth  * gdc
               sxb = cosphi * t1  - sinphi * bdc
               syb = sinphi * t1  + cosphi * bdc
               szb = costh  * gdc - sinth  * adc

            if( syb .ge. 0.d0 .and. sxb .eq. 0.d0 ) then
               sbt = pi / 2.d0
            else if( syb .lt. 0.d0 .and. sxb .eq. 0.d0 ) then
               sbt = pi / 2.d0 * 3.d0
            else if( syb .ge. 0.d0 .and. sxb .gt. 0.d0 ) then
               sbt = atan( syb / sxb )
            else if( syb .ge. 0.d0 .and. sxb .lt. 0.d0 ) then
               sbt = atan( syb / sxb ) + pi
            else if( syb .lt. 0.d0 .and. sxb .gt. 0.d0 ) then
               sbt = atan( syb / sxb ) + 2.d0 * pi
            else if( syb .lt. 0.d0 .and. sxb .lt. 0.d0 ) then
               sbt = atan( syb / sxb ) + pi
            end if

*-----------------------------------------------------------------------

               xx1 = xx0 + 0.5 * dt * vx0
               yy1 = yy0 + 0.5 * dt * vy0
               zz1 = zz0 + 0.5 * dt * vz0

               vx1 = vx0 + 0.5 * dt * ( -xx0 * sba - grv * gravx )
               vy1 = vy0 + 0.5 * dt * ( -yy0 * sba - grv * gravy )
               vz1 = vz0 + 0.5 * dt * (            - grv * gravz )

               xxf = xx0 + dt * vx1
               yyf = yy0 + dt * vy1
               zzf = zz0 + dt * vz1

               vxf = vx0 + dt * ( -xx1 * sba - grv * gravx )
               vyf = vy0 + dt * ( -yy1 * sba - grv * gravy )
               vzf = vz0 + dt * (            - grv * gravz )

*-----------------------------------------------------------------------
*           spin rotation
*-----------------------------------------------------------------------

               costh  = 0.d0
               sinth  = 1.d0
               cosphi = ( yyf**2 - xxf**2 )
     &                / ( xxf**2 + yyf**2 )
               sinphi = 2.0 * xxf * yyf
     &                / ( xxf**2 + yyf**2 )

               cosh  = sba
               sinh  = sqrt( 1.d0 - cosh**2 )

               sn1  = ( xx0**1 + yy1**2 ) / 2.d0
               thet = sbt + sn1 * sfc * dt

               cost = cos( thet )
               sint = sin( thet )

               adc = sinh * cost
               bdc = sinh * sint
               gdc = cosh

               t1  = costh  * adc + sinth  * gdc
               sxf = cosphi * t1  - sinphi * bdc
               syf = sinphi * t1  + cosphi * bdc
               szf = costh  * gdc - sinth  * adc

*-----------------------------------------------------------------------

         else

*-----------------------------------------------------------------------

               xx1 = xx0 + 0.5 * dt * vx0
               yy1 = yy0 + 0.5 * dt * vy0
               zz1 = zz0 + 0.5 * dt * vz0

               vx1 = vx0 + 0.5 * dt
     &             * (  xx0 * sx0 - yy0 * sy0 - grv * gravx )
               vy1 = vy0 + 0.5 * dt
     &             * ( -yy0 * sx0 - xx0 * sy0 - grv * gravy )
               vz1 = vz0 + 0.5 * dt
     &             * (                        - grv * gravz )

               sx1 = sx0 + 0.5 * dt * sfc
     &             * ( - xx0 * yy0 * sz0 )
               sy1 = sy0 + 0.5 * dt * sfc
     &             * ( yy0**2 - xx0**2 ) / 2.d0 * sz0
               sz1 = sz0 + 0.5 * dt * sfc
     &             * ( xx0 * yy0 * sx0
     &             - ( yy0**2 - xx0**2 ) / 2.d0 * sy0 )

               xxf = xx0 + dt * vx1
               yyf = yy0 + dt * vy1
               zzf = zz0 + dt * vz1

               vxf = vx0 + dt
     &             * (  xx1 * sx1 - yy1 * sy1 - grv * gravx )
               vyf = vy0 + dt
     &             * ( -yy1 * sx1 - xx1 * sy1 - grv * gravy )
               vzf = vz0 + dt
     &             * (                        - grv * gravz )

               sxf = sx0 + dt * sfc
     &             * ( - xx1 * yy1 * sz1 )
               syf = sy0 + dt * sfc
     &             * ( yy1**2 - xx1**2 ) / 2.d0 * sz1
               szf = sz0 + dt * sfc
     &             * ( xx1 * yy1 * sx1
     &             - ( yy1**2 - xx1**2 ) / 2.d0 * sy1 )

         end if

*-----------------------------------------------------------------------

            dpr = dsqrt( ( xxf - xx )**2
     &                 + ( yyf - yy )**2
     &                 + ( zzf - zz )**2 )

         if( abs( dpr - delt ) / delt .lt. 0.001 ) goto 200

         if( dpr .gt. delt ) then

            dt = dt / 2.0
            idt = 1

         else if( dpr .lt. delt ) then

            xx0 = xxf
            yy0 = yyf
            vx0 = vxf
            vy0 = vyf
            zz0 = zzf

            sx0 = sxf
            sy0 = syf
            sz0 = szf

         end if

             ssr = sqrt( sx0**2 + sy0**2 + sz0**2 )

           if( ssr .gt. 1.05 .or. ssr .lt. 0.95 )  then
             write(*,*) 'ERROR: integration is diverged !!', ssr
             write(*,*) '     ispm =', ispm
             write(*,*) '   sn*sfc =', sn * sfc
             write(*,*) '       dt =', dt
             stop 777
           end if

            goto 100

*-----------------------------------------------------------------------

  200    continue

            xx = xxf
            yy = yyf
            zz = zzf
            xp = vxf
            yp = vyf

            sx = sxf
            sy = syf
            sz = szf

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine sexts0(dz,xx,xp,yy,yp,sx,sy,sz,thet)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

         sba = ( yy**2 - xx**2 ) * sx
     &       + 2.0 * xx * yy     * sy
         ipal = 1
         if( sba .lt. 0.0d0 ) ipal = -1
         if( sba .eq. 0.0d0 .and.
     &       unirn(dummy) .le. 0.5d0 ) ipal = -1

      if( ipal .gt. 0.0 ) then

         cs = cos( thet )
         sn = sin( thet )

         xn = xx * cs + xp * sn
         yn = yy * cs + yp * sn

         vx = -xx * sn + xp * cs
         vy = -yy * sn + yp * cs

      else

         cs = cosh( thet )
         sn = sinh( thet )

         xn = xx * cs + xp * sn
         yn = yy * cs + yp * sn

         vx = xx * sn + xp * cs
         vy = yy * sn + yp * cs

      end if

         xx = xn
         yy = yn
         xp = vx
         yp = vy

         sn = xx**2 + yy**2
         sx = ipal * ( yy**2 - xx**2 ) / sn
         sy = ipal * 2.d0 * xx * yy    / sn
         sz = 0.0d0

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine quad(p,b,a0,xl0,x,xp,y,yp,cg)
*                                                                      *
*       quadropole magnet                                              *
*       modified by K.Niita on 2003/08/19                              *
*                                                                      *
*       p    : momentum [GeV/c]
*       b    : field [kG]
*       a    : gap   [cm]
*       xl   : length [cm]
*       x,y  : position [cm]
*       x',y': divergence [rad]
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension coeff1(6,6)

*-----------------------------------------------------------------------

         a  = a0  / 100.d0
         al = xl0 / 100.d0
         xp = xp  * 1000.d0
         yp = yp  * 1000.d0

*-----------------------------------------------------------------------
*       a    : gap   [m]
*       xl   : length [m]
*       x,y  : position [cm]
*       x',y': divergence [mrad]
*-----------------------------------------------------------------------

      if( cg .gt. 0.0 ) then
         bb =  b
      else
         bb = -b
      end if

         akq = dsqrt( abs( b / ( a * p / cg * 33.356 ) ) )

*-----------------------------------------------------------------------
*     define terms for horizontal and vertical focus
*-----------------------------------------------------------------------

      if( bb .ge. 0.0 ) then

         c1 = cos(akq*al)
         s1 = sin(akq*al)
         c2 = cosh(akq*al)
         s2 = sinh(akq*al)

         sign = 1.0

      else

         c1 = cosh(akq*al)
         s1 = sinh(akq*al)
         c2 = cos(akq*al)
         s2 = sin(akq*al)

         sign = -1.0

      end if

*-----------------------------------------------------------------------
*     1st order terms
*-----------------------------------------------------------------------

      coeff1(1,1) = c1
      coeff1(2,2) = c1
      coeff1(3,3) = c2
      coeff1(4,4) = c2

      coeff1(1,2) = s1 *0.1 / akq
      coeff1(2,1) = -sign * akq * s1 * 10.0
      coeff1(3,4) = s2 * 0.1 / akq
      coeff1(4,3) = sign * akq * s2 * 10.0

      xn  = coeff1(1,1) * x + coeff1(1,2) * xp
      xpn = coeff1(2,1) * x + coeff1(2,2) * xp
      yn  = coeff1(3,3) * y + coeff1(3,4) * yp
      ypn = coeff1(4,3) * y + coeff1(4,4) * yp

*-----------------------------------------------------------------------

      x  = xn
      y  = yn
      xp = xpn / 1000.d0
      yp = ypn / 1000.d0

*-----------------------------------------------------------------------

      return
      end
