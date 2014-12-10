************************************************************************
*                                                                      *
      subroutine partrs(mark,markp,nbeta,icge,itmak)
*                                                                      *
*       particle transport                                             *
*       and region check                                               *
*       modified by K.Niita on 2003/10/12                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       icge   : error code                                            *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
cFURUTA      include 'mmbank.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      common /eparm/  esmax, esmin, emin(20)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /tcntl/  icntl, inucr
      common /paraj/  mstz(200), parz(200)

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /inout/  in,io
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /cgerr/  nlost, ilost, igerr, icger, ncger
      common /cgstr/  novp, nrovp(3,1000)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /gravit/ grav(3), igrav

*-----------------------------------------------------------------------

      common /elmgrg/ nelctf, nereg, inerc, inert, kelcs

      common /insprd/ initsp
!$OMP THREADPRIVATE(/insprd/)
      common /magreg/ nmreg, ingrc, ingrt, kmags, magtin, imgusr
      common /srsph/  pmphs
!$OMP THREADPRIVATE(/srsph/)
      common /argcns/ kcar, kczs, kcze

cFURUTA      dimension arg(1),zsmcs(1),zemcs(1)
cFURUTA      equivalence ( das, arg, zsmcs, zemcs )
      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      data initsp /0/

*-----------------------------------------------------------------------
*     initialization of spread
*-----------------------------------------------------------------------

         if( initsp .eq. 0 .and. nspred .ne. 0 ) then

            initsp = initsp + 1

            call sprdint

         end if

*-----------------------------------------------------------------------
*     initialization of time frag
*-----------------------------------------------------------------------

            itmak = 0

*-----------------------------------------------------------------------
*     initial jcoll = 0
*-----------------------------------------------------------------------

            jcoll = 0

*-----------------------------------------------------------------------
*        off region particles ( for source particles )
*-----------------------------------------------------------------------

            if( mat .le. -1 ) then

               mark  = -1
               nbeta = 2
               icge  = 0

               ec(ibkec+no) = e(ibke+no)
               xc(ibkxc+no) = x(ibkx+no)
               yc(ibkyc+no) = y(ibky+no)
               zc(ibkzc+no) = z(ibkz+no)

               return

            end if

*-----------------------------------------------------------------------
*        energy cutoff particles
*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

               emint = emin(ityp)

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------

            if( ( e(ibke+no) .le. emint .and. mat .gt. 0 ) .or.
     &            e(ibke+no) .eq. 0.0d0 ) then

               mark  = 1
               nbeta = 3
               icge  = 0

               ec(ibkec+no) = e(ibke+no)
               xc(ibkxc+no) = x(ibkx+no)
               yc(ibkyc+no) = y(ibky+no)
               zc(ibkzc+no) = z(ibkz+no)

               return

            end if

*-----------------------------------------------------------------------
*        check initial cell
*-----------------------------------------------------------------------

               ici = 1

            call gomsor(x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                  u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                  nmed(ibknmd+no),iblz(ibkblz+no),
     &                  mark,markp,ici)

               if( mark .le. -2  ) goto 100

               markp = 1

*-----------------------------------------------------------------------
*     getflt : [ fpl ] sampling of flight path length
*-----------------------------------------------------------------------

            if( icntl .ne. 5 ) then

               call getflt(fpl)

            else

               fpl = 1.0d+10

               goto 30

            end if

*-----------------------------------------------------------------------
*        electron transfer
*-----------------------------------------------------------------------

         if( ( ityp .eq. 12 .or. ityp .eq. 13 )  .and.
     &         mat .ne. 0 .and. fpl .lt. 1.0d+10 ) then

               call eletr(mark,markp,fpl,nbeta,itmak)

               goto 100

         end if

*-----------------------------------------------------------------------
*        charge particle or neutron transfer under magnetic field
*           nmreg    : number of region of magnetic field
*-----------------------------------------------------------------------

         if( nwsprd .gt. 0 .and. nmreg .gt. 0 ) then

                  idsm = ingrt
                  jdsm = 0

                  kdsm = kmags
                  ldsm = 0

            do i = 1, nmreg

                  jj  = 0

                  jdsm = jdsm + 1
                  ntrn = idas(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas(idsm+jdsm)

                  ldsm  = ldsm + 1
                  a_mag = das(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  b_mag = das(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  s_mag = das(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  p_mag = das(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  t_mag = das(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  ldsm  = ldsm + 1
                  u_mag = das(kdsm+ldsm)

               do j = 1, ntrn

                  call tregck(iblz1,ilev1,ilat1,
     &                        mtrn,idas(idsm+jdsm+1),jj,icc)

                  if( icc .ne. 0 ) goto 20

               end do

                  jdsm = jdsm + mtrn

            end do

               goto 25

   20       continue

*-----------------------------------------------------------------------

               isp = nint( s_mag )

            if( isp .ne.  2  .and. isp .ne. 3  .and. isp .ne.   4 .and.
     &          isp .ne. 60  .and. isp .ne.  61 .and.
     &          isp .ne. 62  .and.
     &          isp .ne. 100 .and. isp .ne. 101 .and.
     &          isp .ne. 102 .and. isp .ne. 103 .and.
     &          isp .ne. 104 .and. isp .ne. 106 ) goto 25
            if( ( isp .eq. 2 .or. isp .eq. 3 .or. isp .eq. 4 ) .and.
     &            jtyp .eq. 0 ) goto 25
            if( ( isp .eq.  60 .or. isp .eq.  61 .or.
     &            isp .eq.  62 .or.
     &            isp .eq. 100 .or. isp .eq. 101 .or.
     &            isp .eq. 102 .or. isp .eq. 103 .or.
     &            isp .eq. 104 .or. isp .eq. 106 ) .and.
     &            ityp .ne. 2 ) goto 25

*-----------------------------------------------------------------------
*           Wobbler magnet for charged particles for IHI
*-----------------------------------------------------------------------

            if( p_mag .gt. -1000.0d0 ) then

               if( jtyp .ne. 0 ) then

                  b_mag = sin( pmphs + p_mag ) * b_mag

               end if

            end if

*-----------------------------------------------------------------------

               call magfld(mark,markp,fpl,nbeta,itmak,
     &                     a_mag,b_mag,s_mag,p_mag,t_mag,u_mag)

               goto 100

         end if

   25    continue

*-----------------------------------------------------------------------
*        charge particle transfer under electro magnetic field
*           nereg    : number of region of electro magnetic field
*-----------------------------------------------------------------------

         if( nelctf .gt. 0 .and. nereg .gt. 0 .and. 
     &       jtyp .ne. 0 ) then

                  idsm = inert
                  jdsm = 0

                  kdsm = kelcs
                  ldsm = 0

            do i = 1, nereg

                  jj  = 0

                  jdsm = jdsm + 1
                  ntrn = idas(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas(idsm+jdsm)

                  ldsm  = ldsm + 1
                  s_elf = das(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  s_mgf = das(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  t_elf = das(kdsm+ldsm)
                  ldsm  = ldsm + 2
                  t_mgf = das(kdsm+ldsm)
                  ldsm  = ldsm + 1

               do j = 1, ntrn

                  call tregck(iblz1,ilev1,ilat1,
     &                        mtrn,idas(idsm+jdsm+1),jj,icc)

                  if( icc .ne. 0 ) goto 21

               end do

                  jdsm = jdsm + mtrn

            end do

               goto 26

   21       continue

*-----------------------------------------------------------------------

               call elmgfd(mark,markp,fpl,nbeta,itmak,
     &                     s_elf,s_mgf,t_elf,t_mgf,e_chs)

               goto 100

         end if

   26    continue

*-----------------------------------------------------------------------
*        proton beam transfer by fpl in medium
*        beam spread
*           only for nspred > 0, name(no) = 1, proton, no void
*           nfcs(no) =0,1 not 2 ; no spread for forced collisions
*-----------------------------------------------------------------------

         if( jtyp .ne. 0     .and.
     &       nspred .gt. 0   .and. mat .ne. 0 .and.
     &       name(ibknam+no) .eq. 1 .and.
     &       nfcs(ibknfc+no) .ne. 2 .and.
     &       arg(kcar+mat) .gt. 0.0     ) then

               call sprd(mark,markp,fpl,nbeta,itmak)

               goto 100

         end if

*-----------------------------------------------------------------------
*        particle transfer by fpl in free space or neutron gravity
*-----------------------------------------------------------------------

   30 continue

            if( ityp .eq. 2 .and. igrav .ne. 0 .and.
     &          e(ibke+no) .le. 1.d-6 ) then

               call ngravt(mark,markp,fpl,nbeta,itmak)

            else

               call parfre(mark,markp,fpl,nbeta,itmak)

            end if

               goto 100

*-----------------------------------------------------------------------

  100 continue

*-----------------------------------------------------------------------
*        surface cross: reset the forced collision flag
*-----------------------------------------------------------------------

         if( mark .eq. 0 .or. mark .eq. 2 ) then

            nfcs(ibknfc+no) = 0

         end if

*-----------------------------------------------------------------------
*        surface cross: keep importance of the previous non-void cell
*-----------------------------------------------------------------------

         if( ( mark .eq. 0 .or. mark .eq. 2 ) .and. mat .gt. 0 ) then

            wtnz(ibkwnz+no) = aimp(ityp,iblz1,ilev1,ilat1,ii1)

         end if

*-----------------------------------------------------------------------
*        error in CG/GG ( mark = -2, icge = -2 ) : lost particles
*-----------------------------------------------------------------------

         if( mark .eq. -2 ) then

               icge = -1

               ilost = ilost + 1

                  write(6,'(/''*** Lost particle in CG/GG *** lost ='',
     &                      i3)') ilost

               if( npe .gt. 0 ) 
     &            write(6,'( '' my ip       = '',i3)') me

                  write(6,'( '' nbch ncs no = '',3i10)')
     &                  nobch, nocas, no
                  write(6,'( '' ityp, ktyp, e(no) = '',i4,2x,i9,e17.8)')
     &                  ityp, ktyp, e(ibke+no)
                  write(6,'( ''        mark ='',3i6)')
     &                  mark
                  write(6,'( '' reg. ini fin     ='',5i6)')
     &                  iblz1, iblz2
                  write(6,'( '' mat. ini fin     ='',5i6)')
     &                  idmn(mat), idmn(nmed(ibknmd+no))
                  write(6,'( '' x, y, z :'',3e17.8)')
     &                  x(ibkx+no),y(ibky+no),z(ibkz+no)
                  write(6,'( '' u, v, w :'',3e17.8)')
     &                  u(ibku+no),v(ibkv+no),w(ibkw+no)
            if( ilost .ge. nlost ) then

               write(6,'(/''*** Number of lost error exceeds nlost'',
     &               '' : nlost ='',i3)') nlost
               call parastop( 802 )

            end if

               return

         end if

*-----------------------------------------------------------------------
*        error in CG/GG ( mark = -3, -4, -5 ): region error
*        goto 81, start again after gomsor
*-----------------------------------------------------------------------

         if( mark .le. -3 ) then

                  icge  = icge + 1
                  icger = icger + 1

*-----------------------------------------------------------------------
*              booking the error and write information
*-----------------------------------------------------------------------

                  nb1 = iblz1
                  nb2 = iblz2

                  if( nb1 .lt. nb2 ) then

                     nb3 = nb1
                     nb1 = nb2
                     nb2 = nb3

                  end if
!$OMP CRITICAL (nrovp_crit)
                  do i = 1, novp

                     if( nrovp(1,i) .eq. nb1 .and.
     &                   nrovp(2,i) .eq. nb2 ) then

                        nrovp(3,i) = nrovp(3,i) + 1

                        goto 200

                     end if

                  end do

                  novp = novp + 1

                  if( novp .le. 1000 ) then

                     nrovp(1,novp) = nb1
                     nrovp(2,novp) = nb2
                     nrovp(3,novp) = 1

                  end if

  200             continue
!$OMP END CRITICAL (nrovp_crit)
*-----------------------------------------------------------------------
*                 wirte error information
*-----------------------------------------------------------------------

            if( icger .le. nlost ) then

                  write(6,'(/''*** warning in region check *** no ='',
     &                  i6)') icger

               if( npe .gt. 0 ) 
     &            write(6,'( '' my ip       = '',i3)') me

                  write(6,'( '' nbch ncs no = '',3i10)')
     &                  nobch, nocas, no
                  write(6,'( '' ityp, e(no) = '',i4,1x,e17.8)')
     &                  ityp, e(ibke+no)
                  write(6,'( ''        mark ='',3i6)')
     &                  mark
                  write(6,'( '' reg. ini fin     ='',5i6)')
     &                  iblz1, iblz2
                  write(6,'( '' mat. ini fin     ='',5i6)')
     &                  idmn(mat), idmn(nmed(ibknmd+no))
                  write(6,'( '' x, y, z :'',3e17.8)')
     &                  x(ibkx+no),y(ibky+no),z(ibkz+no)
                  write(6,'( '' xc,yc,zc:'',3e17.8)')
     &                  xc(ibkxc+no),yc(ibkyc+no),zc(ibkzc+no)
                  write(6,'( '' u, v, w :'',3e17.8)')
     &                  u(ibku+no),v(ibkv+no),w(ibkw+no)

            end if

*-----------------------------------------------------------------------
*           icge gt igerr
*-----------------------------------------------------------------------

            if( icge .gt. igerr ) then

                  if( icger .le. nlost ) then

                   write(6,'(''*** failed to recovering ***'')')
        write(6,'(''*** LOST particle! Geometry error occurrs ***'')')

                  endif

                  ncger = ncger + 1

                  icge  = -2

                  return

*-----------------------------------------------------------------------
*           icge > 1
*-----------------------------------------------------------------------

            else if( icge .ge. 1 ) then

               if( mark .ne. -4 ) then

                  x(ibkx+no) = x(ibkx+no) + u(ibku+no) * parz(28)
                  y(ibky+no) = y(ibky+no) + v(ibkv+no) * parz(28)
                  z(ibkz+no) = z(ibkz+no) + w(ibkw+no) * parz(28)

               end if

                  ici   = 0
                  mark  = 1
                  markp = 0

                  call gomsor(x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                        u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                        nmed(ibknmd+no),iblz(ibkblz+no),
     &                        mark,markp,ici)

                     if( mark .le. -2 ) then

                        icge  = -2
                        return

                     end if

                  xc(ibkxc+no) = x(ibkx+no)
                  yc(ibkyc+no) = y(ibky+no)
                  zc(ibkzc+no) = z(ibkz+no)
                  ec(ibkec+no) = e(ibke+no)

                  return

            end if

*-----------------------------------------------------------------------
*        recovered errors
*-----------------------------------------------------------------------

         else if( icge .gt. 0 ) then

                  if( icger .le. nlost )
     &            write(6,'(''*** succeeded in recovering *** icge ='',
     &                      i5)') icge

                  icge = -3

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine parfre(mark,markp,dist,nbeta,itmak)
*                                                                      *
*                                                                      *
*       particle transport by d                                        *
*       and region check                                               *
*       modified by K.Niita on 2005/02/02                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distance                                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
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
      common /eparm/  esmax, esmin, emin(20)
      common /paraj/  mstz(200), parz(200)
      common /tcntl/  icntl, inucr
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx

*-----------------------------------------------------------------------
*        dv = fpl : flight length
*-----------------------------------------------------------------------

                  iff = 0

                  dv  = dist

                  esav = e(ibke+no)
                  tsav = t(ibkt+no)
                  xsav = x(ibkx+no)
                  ysav = y(ibky+no)
                  zsav = z(ibkz+no)

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

               emint = emin(ityp)

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        delt1 : maximum flight step
*        delt  : flight step
*-----------------------------------------------------------------------

                  delt1 = min( parz(27), dist )

               if( jtyp .ne. 0 .and. mat .gt. 0 .and.
     &             nedisp .ne. 0 )
     &             delt1 = min( parz(163), delt1 )

                  delt  = delt1

*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt, dv=dv-delt
*        nbeta = 2 : dv < delt,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

            call gomdis(dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                  u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                  mark,markp,nmed(ibknmd+no),iblz(ibkblz+no))

               if( mark .le. -2 ) goto 600

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

                  ee1 = e(ibke+no)
                  ecc = ee1
                  delt2 = delt

            if( jtyp .ne. 0 .and. mat .ne. 0 .and.
     &          icntl .ne. 5 ) then

                  call rainge(ee1,rng1,mat,ityp,ktyp,jtyp,rtyp)
                  call ecol(ecc,delt2,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

            end if

                  ec(ibkec+no) = ecc

*-----------------------------------------------------------------------
*        next step or collision or stop
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

               if( delt2 .lt. delt ) then

                  nbeta = 3
                  delt = delt2

               end if

               if( ecc .le. emint )  nbeta = 3

*-----------------------------------------------------------------------

               if( nbeta .eq. 1 ) then

                  goto 800

               else

                  goto 700

               end if

*-----------------------------------------------------------------------
*        cross surface or stop
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. mat .ne. 0 .and.
     &          icntl .ne. 5 ) then

                  dpr1 = dpr

                  call ecol(ecc,dpr1,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

                  ec(ibkec+no) = ecc

               if( ecc .le. emint .and. idcyc(ktyp) .eq. 1 ) then

                  call rainge(emint,rngs,mat,ityp,ktyp,jtyp,rtyp)

                  nbeta = 3
                  delt  = rng1 - rngs
                  ec(ibkec+no) = emint

                  goto 700

               else if( dpr1 .lt. dpr .or. ecc .le. emint ) then

                  nbeta = 3
                  delt  = dpr1

                  goto 700

               end if

            end if

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

               call gomnew(1,dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                     xc(ibkxc+no),yc(ibkyc+no),zc(ibkzc+no),
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     ec(ibkec+no),
     &                     nmed(ibknmd+no),iblz(ibkblz+no),
     &                     mark,markp)

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------
*              call tally if itstep = 1,  ncol = 15
*-----------------------------------------------------------------------

               if( itstep .ne. 0 .and. nedisp .ne. 0 ) then

                  call timtrs(itmak,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no)
                  tsav = tc(ibktc+no)
                  xsav = xc(ibkxc+no)
                  ysav = yc(ibkyc+no)
                  zsav = zc(ibkzc+no)

                  ncol = 15

                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------

                  e(ibke+no) = ec(ibkec+no)
                  t(ibkt+no) = tc(ibktc+no)
                  x(ibkx+no) = xc(ibkxc+no)
                  y(ibky+no) = yc(ibkyc+no)
                  z(ibkz+no) = zc(ibkzc+no)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------

  600 continue

                  e(ibke+no) = esav
                  t(ibkt+no) = tsav
                  x(ibkx+no) = xsav
                  y(ibky+no) = ysav
                  z(ibkz+no) = zsav

                  call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sprd(mark,markp,dist,nbeta,itmak)
*                                                                      *
*                                                                      *
*       proton beam spread by Coulomb                                  *
*       and region check                                               *
*       modified by K.Niita on 2010/12/23                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distance                                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
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
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /eparm/  esmax, esmin, emin(20)
      common /paraj/  mstz(200), parz(200)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx

*-----------------------------------------------------------------------
*        dv = dist : flight length
*-----------------------------------------------------------------------

               iff = 0

               dv  = dist

               esav = e(ibke+no)
               tsav = t(ibkt+no)
               xsav = x(ibkx+no)
               ysav = y(ibky+no)
               zsav = z(ibkz+no)

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

               emint = emin(ityp)

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        delt0 : minimum flight step
*        delt1 : maximum flight step
*-----------------------------------------------------------------------

                  delt0 = parz(26)

                  delt1 = min( parz(27), dist )

               if( jtyp .ne. 0 .and. mat .gt. 0 .and.
     &             nedisp .ne. 0 )
     &             delt1 = min( parz(163), delt1 )

                  delt  = delt1

*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt,     dv = dv - delt
*        nbeta = 2 : delt = dv,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------

               if( delt .le. delt0 ) then

                  idelt = 1

               else

                  idelt = 0

               end if

*-----------------------------------------------------------------------
*        transport particle by delt through sprdtrs(delt)
*-----------------------------------------------------------------------

                  uprv = u(ibku+no)
                  vprv = v(ibkv+no)
                  wprv = w(ibkw+no)

               call sprdtrs(delt)

               call gomupp(mark,markp,
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no))

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

               call gomdis(dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     mark,markp,nmed(ibknmd+no),iblz(ibkblz+no))

                  if( mark .le. -2 ) goto 600

*-----------------------------------------------------------------------
*        back to the initial by reducing the delta
*-----------------------------------------------------------------------

         if( delt .ge. dpr .and. idelt .eq. 0 ) then

                  if( nbeta .eq. 1 ) dv = dv + delt

                  u(ibku+no) = uprv
                  v(ibkv+no) = vprv
                  w(ibkw+no) = wprv

                  call gomupp(mark,markp,
     &                        u(ibku+no),v(ibkv+no),w(ibkw+no))

                  delt = max( delt0, dpr * 0.9 )

                  nmed(ibknmd+no) = mat
                  iblz(ibkblz+no) = iblz1

                  goto 500

         end if

*-----------------------------------------------------------------------
*        energy range and final energy
*-----------------------------------------------------------------------

                  ee1 = e(ibke+no)
                  delt2 = delt

                  call rainge(ee1,rng1,mat,ityp,ktyp,jtyp,rtyp)
                  call ecol(ecc,delt2,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

*-----------------------------------------------------------------------
*        next step or collision or stop
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

               if( delt2 .lt. delt ) then

                  nbeta = 3
                  delt = delt2

               end if

               if( ecc .le. emint )  nbeta = 3

                  ec(ibkec+no) = ecc

*-----------------------------------------------------------------------

               if( nbeta .eq. 1 ) then

                  goto 800

               else

                  goto 700

               end if

*-----------------------------------------------------------------------
*        cross surface or stop
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

                  dpr1 = dpr

                  call ecol(ecc,dpr1,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

                  ec(ibkec+no) = ecc

               if( ecc .le. emint .and. idcyc(ktyp) .eq. 1 ) then

                  call rainge(emint,rngs,mat,ityp,ktyp,jtyp,rtyp)

                  nbeta = 3
                  delt  = rng1 - rngs
                  ec(ibkec+no) = emint

                  goto 700

               else if( dpr1 .lt. dpr .or. ecc .le. emint ) then

                  nbeta = 3
                  delt  = dpr1

                  goto 700

               end if

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

               call gomnew(1,dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                     xc(ibkxc+no),yc(ibkyc+no),zc(ibkzc+no),
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     ec(ibkec+no),
     &                     nmed(ibknmd+no),iblz(ibkblz+no),
     &                     mark,markp)

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------
*              call tally if itstep = 1, ncol = 15
*-----------------------------------------------------------------------

               if( nedisp .ne. 0 .and. itstep .ne. 0 ) then

                  call timtrs(itmak,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no)
                  tsav = tc(ibktc+no)
                  xsav = xc(ibkxc+no)
                  ysav = yc(ibkyc+no)
                  zsav = zc(ibkzc+no)

                  ncol = 15

                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------

                  e(ibke+no) = ec(ibkec+no)
                  t(ibkt+no) = tc(ibktc+no)
                  x(ibkx+no) = xc(ibkxc+no)
                  y(ibky+no) = yc(ibkyc+no)
                  z(ibkz+no) = zc(ibkzc+no)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------

  600 continue

                  e(ibke+no) = esav
                  t(ibkt+no) = tsav
                  x(ibkx+no) = xsav
                  y(ibky+no) = ysav
                  z(ibkz+no) = zsav

                  call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ngravt(mark,markp,dist,nbeta,itmak)
*                                                                      *
*                                                                      *
*       neutron with gravity                                           *
*       and region check                                               *
*       modified by K.Niita on 2010/12/23                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distance                                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
cFURUTA      include 'mmbank.inc'

      parameter ( grvc = 980.665d-18 )
      parameter ( rlit = 29.97925d0 )

*-----------------------------------------------------------------------


      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /eparm/  esmax, esmin, emin(20)
      common /paraj/  mstz(200), parz(200)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /gravit/ grav(3), igrav

*-----------------------------------------------------------------------
*        dv = dist : flight length
*-----------------------------------------------------------------------

               iff = 0

               dv  = dist

               esav = e(ibke+no)
               tsav = t(ibkt+no)
               xsav = x(ibkx+no)
               ysav = y(ibky+no)
               zsav = z(ibkz+no)

*-----------------------------------------------------------------------
*        delt0 : minimum flight step
*        delt1 : maximum flight step
*-----------------------------------------------------------------------

                  delt0 = parz(26)

                  delt1 = min( parz(27), dist )
                  delt1 = min( parz(163), delt1 )

                  delt  = delt1

*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt,     dv = dv - delt
*        nbeta = 2 : delt = dv,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------

               if( delt .le. delt0 ) then

                  idelt = 1

               else

                  idelt = 0

               end if

*-----------------------------------------------------------------------
*        transport particle by delt
*-----------------------------------------------------------------------

                  delts = delt

                  uprv = u(ibku+no)
                  vprv = v(ibkv+no)
                  wprv = w(ibkw+no)

*-----------------------------------------------------------------------
*           gravity
*-----------------------------------------------------------------------

               xc(ibkxc+no) = x(ibkx+no) + delt * u(ibku+no)
               yc(ibkyc+no) = y(ibky+no) + delt * v(ibkv+no)
               zc(ibkzc+no) = z(ibkz+no) + delt * w(ibkw+no)

               vel = sqrt( 2.0 * e(ibke+no) / rtyp ) * rlit
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

*-----------------------------------------------------------------------

               call gomupp(mark,markp,
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no))

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

               call gomdis(dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     mark,markp,nmed(ibknmd+no),iblz(ibkblz+no))

                  if( mark .le. -2 ) goto 600

*-----------------------------------------------------------------------
*        back to the initial by reducing the delta
*-----------------------------------------------------------------------

         if( delt .ge. dpr .and. idelt .eq. 0 ) then

                  if( nbeta .eq. 1 ) dv = dv + delts

                  u(ibku+no) = uprv
                  v(ibkv+no) = vprv
                  w(ibkw+no) = wprv

               call gomupp(mark,markp,
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no))

                  delt = max( delt0, dpr * 0.9 )

                  nmed(ibknmd+no) = mat
                  iblz(ibkblz+no) = iblz1

                  goto 500

         end if

*-----------------------------------------------------------------------
*        next step or collisions
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

            if( nbeta .eq. 1 ) then

               goto 800

            else

               goto 700

            end if

*-----------------------------------------------------------------------
*        cross surface and search new cell
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

               call gomnew(1,dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                     xc(ibkxc+no),yc(ibkyc+no),zc(ibkzc+no),
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     ec(ibkec+no),
     &                     nmed(ibknmd+no),iblz(ibkblz+no),
     &                     mark,markp)

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------
*              call tally if itstep = 1, ncol = 15
*-----------------------------------------------------------------------

               if( itstep .ne. 0 ) then

                  call timtrs(itmak,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no)
                  tsav = tc(ibktc+no)
                  xsav = xc(ibkxc+no)
                  ysav = yc(ibkyc+no)
                  zsav = zc(ibkzc+no)

                  ncol = 15

                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------

                  e(ibke+no) = ec(ibkec+no)
                  t(ibkt+no) = tc(ibktc+no)
                  x(ibkx+no) = xc(ibkxc+no)
                  y(ibky+no) = yc(ibkyc+no)
                  z(ibkz+no) = zc(ibkzc+no)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------

  600 continue

                  e(ibke+no) = esav
                  t(ibkt+no) = tsav
                  x(ibkx+no) = xsav
                  y(ibky+no) = ysav
                  z(ibkz+no) = zsav

                  call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine magfld(mark,markp,dist,nbeta,itmak,
     &                  a_mag,b_mag,s_mag,p_mag,t_mag,u_mag)
*                                                                      *
*                                                                      *
*       charge particle or neutron transfer under magnetic field       *
*       and region check                                               *
*       modified by K.Niita on 2005/03/30                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distace                                               *
*       a_mag  : magnet gap(mm)                                        *
*       b_mag  : magnet field at pole tip  [kG]                        *
*       s_mag  : speicies of magnet dypole:2, quad:4, sext:6, oct:8    *
*       p_mag  : phase of magnetic field for charge particle           *
*                polarization of neuteron                              *
*       t_mag  : transform id                                          *
*       u_mag  : critical time of time dependent magnetic field        *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
cFURUTA      include 'mmbank.inc'

*-----------------------------------------------------------------------

      parameter ( rlit = 29.97925d0 )

      common /inout/  in,io
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common /paraj/  mstz(200), parz(200)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /eparm/  esmax, esmin, emin(20)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /magreg/ nmreg, ingrc, ingrt, kmags, magtin, imgusr

*-----------------------------------------------------------------------
*        dv = dist : distance
*-----------------------------------------------------------------------

               iff = 0

               dv  = dist

               esav = e(ibke+no)
               tsav = t(ibkt+no)
               xsav = x(ibkx+no)
               ysav = y(ibky+no)
               zsav = z(ibkz+no)

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

               emint = emin(ityp)

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        delt0 : minimum flight step
*        delt1 : maximum flight step ( parz(109) for deltg )
*                                    ( parz(27)  for normal )
*                                    ( parz(163) for deltc )
*                                    ( parz(171) for deltt msec)
*-----------------------------------------------------------------------

                  delt0 = parz(26)

                  delt1 = min( parz(27), dist )
                  delt1 = min( parz(109), delt1 )

               if( jtyp .ne. 0 .and. mat .gt. 0 .and.
     &             nedisp .ne. 0 )
     &             delt1 = min( parz(163), delt1 )

               if( u_mag .gt. -1.0d+9 .and. esav .lt. 1.e-6 ) then

                     timd = delt1 * sqrt( rtyp / 2.0 / esav )
     &                    / rlit * 1.e-6

                  if( timd .gt. parz(171) ) then

                     delt1 = parz(171) / sqrt( rtyp / 2.0 / esav )
     &                     * rlit / 1.e-6

                  end if

               end if


                  delt  = delt1

*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
*           Time depedent magnet
*-----------------------------------------------------------------------

            if( u_mag .gt. -1.0d+9 ) then

                  ptime = abs(t(ibkt+no)) * 1.e-6

               if( imgusr .eq. 1 ) then

                  call usrmgt1(b_mag,u_mag,ptime,o_mag)

               else if( imgusr .eq. 2 ) then

                  call usrmgt2(b_mag,u_mag,ptime,o_mag)

               end if

            else

                  o_mag = b_mag

            end if

*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt,     dv = dv - delt
*        nbeta = 2 : delt = dv,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------

               if( delt .le. delt0 ) then

                  idelt = 1

               else

                  idelt = 0

               end if

*-----------------------------------------------------------------------
*        transport particle by delt(distance) under magnetic field
*-----------------------------------------------------------------------

                  delts = delt

                  uprv = u(ibku+no)
                  vprv = v(ibkv+no)
                  wprv = w(ibkw+no)

                  sxsv = spx(ibkspx+no)
                  sysv = spy(ibkspy+no)
                  szsv = spz(ibkspz+no)

               call magtrs(a_mag,o_mag,s_mag,p_mag,t_mag,
     &                     delt,dpr,udir,vdir,wdir)

               call gomupp(mark,markp,
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no))

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

               call gomdis(dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     mark,markp,nmed(ibknmd+no),iblz(ibkblz+no))

                  if( mark .le. -2 ) goto 600

*-----------------------------------------------------------------------
*        back to the initial by reducing the delta
*-----------------------------------------------------------------------

         if( delt .ge. dpr .and. idelt .eq. 0 ) then

                  if( nbeta .eq. 1 ) dv = dv + delts

                  u(ibku+no) = uprv
                  v(ibkv+no) = vprv
                  w(ibkw+no) = wprv

                  spx(ibkspx+no) = sxsv
                  spy(ibkspy+no) = sysv
                  spz(ibkspz+no) = szsv

               call gomupp(mark,markp,
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no))

                  delt = max( delt0, dpr * 0.9 )

                  nmed(ibknmd+no) = mat
                  iblz(ibkblz+no) = iblz1

                  goto 500

         end if

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

                  ee1 = e(ibke+no)
                  ecc = ee1
                  delt2 = delt

            if( jtyp .ne. 0 .and. mat .ne. 0 ) then

                  call rainge(ee1,rng1,mat,ityp,ktyp,jtyp,rtyp)
                  call ecol(ecc,delt2,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

            end if

                  ec(ibkec+no) = ecc

*-----------------------------------------------------------------------
*        next step or collision or stop
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

               if( delt2 .lt. delt ) then

                  nbeta = 3
                  delt = delt2

               end if

               if( ecc .le. emint )  nbeta = 3

*-----------------------------------------------------------------------
*              for dipole, delt is sometimes changed ( small value )
*-----------------------------------------------------------------------

               if( nbeta .eq. 1 ) then

                  dv = dv + ( delts - delt )

                  goto 800

               else if( nbeta .eq. 2 ) then

                  dv = dv - delt

                  if( dv .le. 0.0d0 ) goto 700

                  nbeta = 1

                  goto 800

               else

                  goto 700

               end if

*-----------------------------------------------------------------------
*        cross surface or stop
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. mat .ne. 0 ) then

                  dpr1 = dpr

                  call ecol(ecc,dpr1,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

                  ec(ibkec+no) = ecc

               if( ecc .le. emint .and. idcyc(ktyp) .eq. 1 ) then

                  call rainge(emint,rngs,mat,ityp,ktyp,jtyp,rtyp)

                  nbeta = 3
                  delt  = rng1 - rngs
                  ec(ibkec+no) = emint

                  goto 700

               else if( dpr1 .lt. dpr .or. ecc .le. emint ) then

                  nbeta = 3
                  delt  = dpr1

                  goto 700

               end if

            end if

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

               call gomnew(1,dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                     xc(ibkxc+no),yc(ibkyc+no),zc(ibkzc+no),
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     ec(ibkec+no),
     &                     nmed(ibknmd+no),iblz(ibkblz+no),
     &                     mark,markp)

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

                  u(ibku+no) = udir
                  v(ibkv+no) = vdir
                  w(ibkw+no) = wdir

               call gomupp(mark,markp,
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no))

*-----------------------------------------------------------------------
*              call tally if itstep = 1, ncol = 15
*-----------------------------------------------------------------------

               if( itstep .ne. 0 ) then

                  call timtrs(itmak,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no)
                  tsav = tc(ibktc+no)
                  xsav = xc(ibkxc+no)
                  ysav = yc(ibkyc+no)
                  zsav = zc(ibkzc+no)

                  ncol = 15

                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------

                  e(ibke+no) = ec(ibkec+no)
                  t(ibkt+no) = tc(ibktc+no)
                  x(ibkx+no) = xc(ibkxc+no)
                  y(ibky+no) = yc(ibkyc+no)
                  z(ibkz+no) = zc(ibkzc+no)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------

  600 continue

                  e(ibke+no) = esav
                  t(ibkt+no) = tsav
                  x(ibkx+no) = xsav
                  y(ibky+no) = ysav
                  z(ibkz+no) = zsav

                  call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine elmgfd(mark,markp,dist,nbeta,itmak,
     &                  s_elf,s_mgf,t_elf,t_mgf,e_chs)
*                                                                      *
*                                                                      *
*       charge particle under electro magnetic field                   *
*       and region check                                               *
*       modified by K.Niita on 2011/01/10                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distace                                               *
*       s_elf  : electric field (kV/cm)                                *
*       s_mgf  : dipole magnet field  [kG]                             *
*       t_elf  : transform for electirc field                          *
*       t_mgf  : transform for magnetic field                          *
*       e_chs  : charge state for this region                          *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************
cKN
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
cKN   include 'mmbank.inc'

*-----------------------------------------------------------------------

      parameter ( rlit = 29.97925d0 )

      common /inout/  in,io
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common /paraj/  mstz(200), parz(200)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /eparm/  esmax, esmin, emin(20)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx

      common /elmgrg/ nelctf, nereg, inerc, inert, kelcs

*-----------------------------------------------------------------------
*        transform ( electric field = x,  magnetic field = y )
*-----------------------------------------------------------------------
cKN 2013/07/29
               itre = nint( t_elf )
               call trnsuv(1.d0,0.d0,0.d0,elcx,elcy,elcz,itre)

               itrm = nint( t_mgf )
               call trnsuv(0.d0,1.d0,0.d0,bmgx,bmgy,bmgz,itrm)

*-----------------------------------------------------------------------
*        charge state
*-----------------------------------------------------------------------

               jchg = nzst(ibkzst+no)
               chgp = dble( jchg )

*-----------------------------------------------------------------------
*        constant
*-----------------------------------------------------------------------

               elct = chgp * s_elf * 1.d-3 * rlit
               bmgt = chgp * s_mgf / 3.3356 * rlit

*-----------------------------------------------------------------------
*        dv = dist : distance
*-----------------------------------------------------------------------

               iff = 0

               dv  = dist

               esav = e(ibke+no)
               tsav = t(ibkt+no)
               xsav = x(ibkx+no)
               ysav = y(ibky+no)
               zsav = z(ibkz+no)

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

               emint = emin(ityp)

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        delt0 : minimum flight step
*        delt1 : maximum flight step ( parz(109) for deltg )
*                                    ( parz(27)  for normal )
*                                    ( parz(163) for deltc )
*-----------------------------------------------------------------------

                  delt0 = parz(26)

                  delt1 = min( parz(27), dist )
                  delt1 = min( parz(109), delt1 )

               if( jtyp .ne. 0 .and. mat .gt. 0 .and.
     &             nedisp .ne. 0 )
     &            delt1 = min( parz(163), delt1 )

                  delt  = delt1

*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt,     dv = dv - delt
*        nbeta = 2 : delt = dv,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------

               if( delt .le. delt0 ) then

                  idelt = 1

               else

                  idelt = 0

               end if

*-----------------------------------------------------------------------
*        transport particle by delt(distance) under magnetic field
*-----------------------------------------------------------------------

                  delts = delt

                  uprv = u(ibku+no)
                  vprv = v(ibkv+no)
                  wprv = w(ibkw+no)

*-----------------------------------------------------------------------

               call elmtrs(elct,elcx,elcy,elcz,
     &                     bmgt,bmgx,bmgy,bmgz,chgp,
     &                     delt,dpr,udir,vdir,wdir)

*-----------------------------------------------------------------------

               call gomupp(mark,markp,
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no))

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

               call gomdis(dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     mark,markp,nmed(ibknmd+no),iblz(ibkblz+no))

                  if( mark .le. -2 ) goto 600

*-----------------------------------------------------------------------
*        back to the initial by reducing the delta
*-----------------------------------------------------------------------

         if( delt .ge. dpr .and. idelt .eq. 0 ) then

                  if( nbeta .eq. 1 ) dv = dv + delts

                  u(ibku+no) = uprv
                  v(ibkv+no) = vprv
                  w(ibkw+no) = wprv

               call gomupp(mark,markp,
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no))

                  delt = max( delt0, dpr * 0.9 )

                  nmed(ibknmd+no) = mat
                  iblz(ibkblz+no) = iblz1

                  goto 500

         end if

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

                  ee1 = e(ibke+no)
                  ecc = ee1
                  delt2 = delt

            if( jtyp .ne. 0 .and. mat .ne. 0 ) then

                  call rainge(ee1,rng1,mat,ityp,ktyp,jtyp,rtyp)
                  call ecol(ecc,delt2,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

            end if

                  ec(ibkec+no) = ecc

*-----------------------------------------------------------------------
*        next step or collision or stop
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

               if( delt2 .lt. delt ) then

                  nbeta = 3
                  delt = delt2

               end if

               if( ecc .le. emint )  nbeta = 3

*-----------------------------------------------------------------------
*              for dipole, delt is sometimes changed ( small value )
*-----------------------------------------------------------------------

               if( nbeta .eq. 1 ) then

                  dv = dv + ( delts - delt )

                  goto 800

               else if( nbeta .eq. 2 ) then

                  dv = dv - delt

                  if( dv .le. 0.0d0 ) goto 700

                  nbeta = 1

                  goto 800

               else

                  goto 700

               end if

*-----------------------------------------------------------------------
*        cross surface or stop
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. mat .ne. 0 ) then

                  dpr1 = dpr

                  call ecol(ecc,dpr1,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

                  ec(ibkec+no) = ecc

               if( ecc .le. emint .and. idcyc(ktyp) .eq. 1 ) then

                  call rainge(emint,rngs,mat,ityp,ktyp,jtyp,rtyp)

                  nbeta = 3
                  delt  = rng1 - rngs
                  ec(ibkec+no) = emint

                  goto 700

               else if( dpr1 .lt. dpr .or. ecc .le. emint ) then

                  nbeta = 3
                  delt  = dpr1

                  goto 700

               end if

            end if

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

               call gomnew(1,dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                     xc(ibkxc+no),yc(ibkyc+no),zc(ibkzc+no),
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     ec(ibkec+no),
     &                     nmed(ibknmd+no),iblz(ibkblz+no),
     &                     mark,markp)

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

                  u(ibku+no) = udir
                  v(ibkv+no) = vdir
                  w(ibkw+no) = wdir

               call gomupp(mark,markp,
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no))

*-----------------------------------------------------------------------
*              call tally if itstep = 1, ncol = 15
*-----------------------------------------------------------------------

               if( itstep .ne. 0 ) then

                  call timtrs(itmak,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no)
                  tsav = tc(ibktc+no)
                  xsav = xc(ibkxc+no)
                  ysav = yc(ibkyc+no)
                  zsav = zc(ibkzc+no)

                  ncol = 15

                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------

                  e(ibke+no) = ec(ibkec+no)
                  t(ibkt+no) = tc(ibktc+no)
                  x(ibkx+no) = xc(ibkxc+no)
                  y(ibky+no) = yc(ibkyc+no)
                  z(ibkz+no) = zc(ibkzc+no)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------

  600 continue

                  e(ibke+no) = esav
                  t(ibkt+no) = tsav
                  x(ibkx+no) = xsav
                  y(ibky+no) = ysav
                  z(ibkz+no) = zsav

                  call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine timtrs(itmak,mark)
*                                                                      *
*       time transport of the particle                                 *
*       and check the time limit                                       *
*       modified by K.Niita on 2004/01/07                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       itmak  : =0, less than tmax                                    *
*                =1, greater than tmax                                 *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
cFURUTA      include 'mmbank.inc'

      parameter ( rlit = 29.97925d0 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common /tparm/  tmax(20)

*-----------------------------------------------------------------------
*     t(ibkt+no) < 0 : timer stopped the clock.
*-----------------------------------------------------------------------

               itmak = 0

           if( t(ibkt+no) .lt. 0.d0 ) return

*-----------------------------------------------------------------------

               ekin = ( ec(ibkec+no) + e(ibke+no) ) / 2.0
               dist = sqrt( ( xc(ibkxc+no) - x(ibkx+no) )**2
     &                    + ( yc(ibkyc+no) - y(ibky+no) )**2
     &                    + ( zc(ibkzc+no) - z(ibkz+no) )**2 )

                  timd = 0.0

               if( ekin .gt. 0.1 .or. rtyp .eq. 0.0d0 ) then

                  timd = dist * ( ekin + rtyp )
     &                 / sqrt( ekin * ( ekin + 2.0 * rtyp ) )
     &                 / rlit

               else if( ekin .gt. 0.0 ) then

                  timd = dist * sqrt( rtyp / 2.0 / ekin ) / rlit

               end if

                  tc(ibktc+no) = t(ibkt+no) + timd

*-----------------------------------------------------------------------
*           time over
*-----------------------------------------------------------------------

            if( tc(ibktc+no) .gt. tmax(ityp) ) then

                  timd = ( tc(ibktc+no) - tmax(ityp) )

                  dd = dist

               if( ekin .gt. 0.1 .or. rtyp .eq. 0.0d0 ) then

                  dd = dist - timd / ( ekin + rtyp )
     &                      * sqrt( ekin * ( ekin + 2.0 * rtyp ) )
     &                      * rlit

               else if( ekin .gt. 0.0 ) then

                  dd = dist - timd / sqrt( rtyp / 2.0 / ekin ) * rlit

               end if

                  xc(ibkxc+no) = x(ibkx+no) + u(ibku+no) * dd
                  yc(ibkyc+no) = y(ibky+no) + v(ibkv+no) * dd
                  zc(ibkzc+no) = z(ibkz+no) + w(ibkw+no) * dd

                  itmak = 1
                  mark  = 1

            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine eletr(mark,markp,dist,nbeta,itmak)
*                                                                      *
*                                                                      *
*       electron transfer                                              *
*       and region check                                               *
*       modified by K.Niita on 2003/10/12                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distance                                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
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
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /eparm/  esmax, esmin, emin(20)

*-----------------------------------------------------------------------

cKN 2012/10/31
      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, noz, mtel
!$OMP THREADPRIVATE(/elect/)

c      common /qscom/qsexpect
c!$OMP THREADPRIVATE(/qscom/)

*-----------------------------------------------------------------------

      data delt0 /1.d-5/
      data delt1 /1.d-6/

*-----------------------------------------------------------------------

               nbeta = 2

*-----------------------------------------------------------------------
*        delt = dist : range of electron in ns substeop
*        delt0 : small distance before boundary
*-----------------------------------------------------------------------

               delt  = dist

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

               call gomdis(dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     mark,markp,nmed(ibknmd+no),iblz(ibkblz+no))

                  if( mark .le. -2 ) return

*-----------------------------------------------------------------------
*        cross the boundary
*-----------------------------------------------------------------------
*             mark: description
*              -1 : outgoing to the void region
*               0 : pass the forward surface
*               2 : reflect surface
*-----------------------------------------------------------------------

         if( delt .ge. 1.0d+10 .or.
     &     ( delt .ge. dpr .and. dpr .lt. delt0 ) ) then

               call gomnew(1,dpr,x(ibkx+no),y(ibky+no),z(ibkz+no),
     &                     xc(ibkxc+no),yc(ibkyc+no),zc(ibkzc+no),
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     ec(ibkec+no),
     &                     nmed(ibknmd+no),iblz(ibkblz+no),
     &                     mark,markp)

               ec(ibkec+no) = e(ibke+no)

               ns = 0

               return

*-----------------------------------------------------------------------
*        collision at just before the boundary
*-----------------------------------------------------------------------

         else if( delt .ge. dpr .and. dpr .ge. delt1 ) then

               call gomupr(mark,markp,dpr-delt1)

               d = dpr - delt1

*-----------------------------------------------------------------------
*        collision at delt
*-----------------------------------------------------------------------

         else if( delt .lt. dpr ) then

               call gomupr(mark,markp,delt)

               d = delt

         end if

*-----------------------------------------------------------------------
*        change u and energy
*-----------------------------------------------------------------------

               uint(1) = u(ibku+no)
               uint(2) = v(ibkv+no)
               uint(3) = w(ibkw+no)

               am = 1.
               call defelc(am,d,1.d0,1.d0,1,
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no),
     &                     uint(1),uint(2),uint(3))

               call gomupp(mark,markp,
     &                     u(ibku+no),v(ibkv+no),w(ibkw+no))

               delc = d
               eint = e(ibke+no)

               ec(ibkec+no) = max( 0.0d0, e(ibke+no) - d * qs )
               noz = n1
               qo  = qs

*-----------------------------------------------------------------------
*        time evolution
*-----------------------------------------------------------------------

               call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sprdtrs(delt)
*                                                                      *
*       proton beam spread by Coulomb                                  *
*       modified by K.Niita on 2003/10/07                              *
*                                                                      *
*       transport particle by delt(distance)                           *
*                                                                      *
*         initial values                                               *
*                                                                      *
*              e(no), x(no), y(no), z(no), u(no), v(no), w(no)         *
*                                                                      *
*           final values                                               *
*                                                                      *
*              xc(no), yc(no), zc(no), u(no), v(no), w(no)             *
*                                                                      *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
cFURUTA      include 'mmbank.inc'

      parameter ( pi  = 3.1415926535898d0 )

*-----------------------------------------------------------------------

      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common /argcns/ kcar, kczs, kcze

cFURUTA      dimension arg(1),zsmcs(1),zemcs(1)
cFURUTA      equivalence ( das, arg, zsmcs, zemcs )
      common /kmat1g/ kmat(kvlmax)
      dimension dnel(1), denh(1), zz(1), a(1), den(1)
      equivalence ( das, dnel, denh, zz, a, den )

*-----------------------------------------------------------------------

         if( delt .le. 0.0d0 ) then

               xc(ibkxc+no) = x(ibkx+no)
               yc(ibkyc+no) = y(ibky+no)
               zc(ibkzc+no) = z(ibkz+no)

               return

         end if

*-----------------------------------------------------------------------
*        nspred = 1 : NMTC original
*               = 2 : First order Moliere model ( by Meigo and Harada )
*   not comlete = 3 : Third order Moliere model ( by Meigo )
*-----------------------------------------------------------------------

         if( nspred .eq. 1 .or. nspred .eq. 2 .or.
     &     ( nspred .eq. 10 .and.
     &       ityp .ne. 1 .and. ityp .lt. 15 ) ) then

               sqxox = arg(kcar+mat) * dsqrt(delt)

            if( nspred .eq. 1 ) then

               sigth = 0.369 * abs(jtyp)
     &               / (e(ibke+no)*(1.+rtyp/(e(ibke+no)+rtyp)))
     &               * sqxox / 0.866

               sigx = delt * sigth

            else if( nspred .eq. 2 .or. nspred .eq. 10 ) then

               sqxox1= sqxox / sqrt(431.273)
*                                   431.273 = 716.4 * 6.02e23 /1e24

               sigth = 13.6 * abs(jtyp)
     &               / (e(ibke+no)*(1.+rtyp/(e(ibke+no)+rtyp)))
     &               * sqxox1 * ( 1.0 + 0.038 * log( sqxox1**2 ) )

               sigx = delt * sigth

            end if

*-----------------------------------------------------------------------
cKN
               sigx = min( delt, sigx )
cKN
*-----------------------------------------------------------------------
cKN 2005/02/07

  560          r1  = gaurn(x1)
               r2  = gaurn(x2)

               rm = sigx * sqrt( r1**2 + r2**2 )

               aa  = delt**2 - rm**2

               if( aa .lt. 0.0d0 ) goto 560

               thet = unirn(dummy) * 2.0 * pi
               cs   = cos(thet)
               sn   = sin(thet)

               xpr = rm * cs
               ypr = rm * sn

cKN 2005/02/07
*-----------------------------------------------------------------------

         else if( nspred .eq. 3 ) then

               sqxox = arg(kcar+mat) * dsqrt(delt)

               p   = sqrt( e(ibke+no)**2 + 2. * rtyp * e(ibke+no) )
               bet = p / rtyp / sqrt( 1. + ( p / rtyp )**2 )

  570          continue

               sud = denh(kmat(mat)+2)
               sum = sud * 2.0
     &             * ( log( 1. + 3.34
     &             * ( 1.0 / 137.036 / bet )**2 ) )

               nel = nint( dnel(kmat(mat)+1) )

            do i = 1, nel

               zzik =  zz(kmat(mat)+(i-1)*3+3)
                aik =   a(kmat(mat)+(i-1)*3+4)
               dnik = den(kmat(mat)+(i-1)*3+5)

               sum = sum + dnik * zzik * ( zzik + 1.0 )
     &             * ( log( 1. + 3.34
     &             * ( zzik / 137.036 / bet )**2 ) )

               sud = sud + dnik

            end do

               zxmcs = sum / sud

cKN
               chiccmcs = 0.26057 * zsmcs(kczs+mat) * sud

c              chiccmcs = 0.156915 * facfms * zsmcs
c    &                        * mparm(3,krad) / mparm(2,krad)
c              densat = avagn * mparm(3,krad) / mparm(2,krad)
c                       avagn/6.0221367e23/
cKN

               chick2 = chiccmcs / ( bet * p )**2 * delt

               chia2 = 2.016e-5 / p**2
     &               * exp( ( zxmcs - zemcs(kcze+mat) )
     &               / zsmcs(kczs+mat) )

               omega = chick2 / ( 1.167 * chia2 )

            if( omega .lt. 2.8 ) then

               call ruthscat(omega,chia2,xpr,ypr)

            else

               call moliere(omega,chick2,xpr,ypr)

            end if

               xpr = delt * xpr
               ypr = delt * ypr

               aa  = delt**2 - xpr**2 - ypr**2

            if( aa .lt. 0.0d0 ) goto 570

*-----------------------------------------------------------------------
*        angle straggling by ATIMA
*-----------------------------------------------------------------------

         else if( nspred .eq. 10 ) then

               ene = e(ibke+no)

            if( ityp .ge. 15 ) then

               ap  = dble( ktyp - ktyp / 1000000 * 1000000 )
               zp  = dble( ktyp / 1000000 )

            else if( ityp .eq. 1 ) then

               ap  = 1.d0
               zp  = 1.d0

            end if

               iway = 4

               call atima(ap,zp,ene,rtyp,mat,rng,delt,rm,iway)

               aa  = delt**2 - rm**2

               thet = unirn(dummy) * 2.0 * pi
               cs   = cos(thet)
               sn   = sin(thet)

               xpr = rm * cs
               ypr = rm * sn

         end if

*-----------------------------------------------------------------------

               zpr = dsqrt(aa)

               csth = w(ibkw+no)
               snth = u(ibku+no)**2 + v(ibkv+no)**2

            if( snth .le. 0.0d0 ) then

               csphi = 1.0d0
               snphi = 0.0d0

            else

               snth  = dsqrt(snth)
               csphi = u(ibku+no) / snth
               snphi = v(ibkv+no) / snth

            end if

               cord = csth * xpr + snth * zpr

               xc(ibkxc+no) = x(ibkx+no) + csphi * cord - snphi * ypr
               yc(ibkyc+no) = y(ibky+no) + snphi * cord + csphi * ypr
               zc(ibkzc+no) = z(ibkz+no) - snth  * xpr  + csth  * zpr

               u(ibku+no) = ( csphi * cord - snphi * ypr ) / delt
               v(ibkv+no) = ( snphi * cord + csphi * ypr ) / delt
               w(ibkw+no) = ( -snth * xpr  + csth  * zpr ) / delt

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine sprdint
*                                                                      *
*                                                                      *
*       initialization of proton beam spread by Coulomb                *
*       modified by K.Niita on 10/02/2000                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       arg(kcar+i): out put in common                                 *
*                                                                      *
************************************************************************
      use MEMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx
      common /cparm/ maxbch,maxcas
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /argcns/ kcar, kczs, kcze

cFURUTA      dimension arg(1),zsmcs(1),zemcs(1)
cFURUTA      equivalence ( das, arg, zsmcs, zemcs )
      common /kmat1g/ kmat(kvlmax)
      dimension dnel(1), denh(1), zz(1), a(1), den(1)
      equivalence ( das, dnel, denh, zz, a, den )

*-----------------------------------------------------------------------

      if( nspred .eq. 0 ) return

*-----------------------------------------------------------------------
*     initialization of Coulomb spreading
*        nspred = 1 : NMTC original
*               = 2 : First order Moliere model ( by Meigo and Harada )
*   not comlete = 3 : Third order Moliere model ( by Meigo )
*-----------------------------------------------------------------------

      do 80 k = 1, mxmat

         if( nspred .eq. 1 ) then

            sumarg = denh(kmat(k)+2) * 21.132d+0 * 0.498d+0

         else if( nspred .eq. 2 .or. nspred .eq. 10 ) then

            sumarg = denh(kmat(k)+2) * 11.319d+0

         else if( nspred .eq. 3 ) then

            sumarg = denh(kmat(k)+2) * 11.319d+0

            sumzzs= denh(kmat(k)+2)* 2.
            sumzze= 0.0
            sumden= denh(kmat(k)+2)

         end if

            iii = nint( dnel(kmat(k)+1) )

         do 75 i = 1, iii

               zzik = zz(kmat(k)+(i-1)*3+3)
                aik =  a(kmat(k)+(i-1)*3+4)
               dnik = den(kmat(k)+(i-1)*3+5)

            if( nspred .eq. 1 ) then

               sumarg = sumarg +
     &                ( dnik * zzik * (zzik+1.d+0)
     &                * ( 10.566d+0 - 0.333d+0 * log( zzik * aik ) )
     &                * 0.498d+0 )

            else if( nspred .eq. 2 .or. nspred .eq. 10 ) then

               sumarg = sumarg +
     &                  dnik * zzik * (zzik+1.d+0)
     &                * log( 287.d0 / sqrt( zzik ) )

            else if( nspred .eq. 3 ) then

               sumarg = sumarg +
     &                  dnik * zzik * (zzik+1.d+0)
     &                * log( 287.d0 / sqrt( zzik ) )

               sumzzs = sumzzs + dnik * zzik * (zzik+1.d+0)
               sumzze = sumzze + dnik * zzik * (zzik+1.d+0)
     &                * log(zzik)*(-2./3.)
               sumden = sumden + dnik

            end if

   75    continue

               arg(kcar+k) = sqrt(sumarg)

            if( nspred .eq. 3 ) then

               if( sumden .gt. 0.0 ) then

                  zsmcs(kczs+k) = sumzzs / sumden
                  zemcs(kcze+k) = sumzze / sumden

               else

                  zsmcs(kczs+k) = 0.0
                  zemcs(kcze+k) = 0.0

               end if

            end if

   80 continue

*-----------------------------------------------------------------------

      return
      end


C add by meigo
c  *********************************************************
      subroutine moliere(omega,chic2,dthxz,dthyz)
c
c  determine scattering angle from moliere distribution
c
c  derived from geant program  gmolie
c
c     omega : mean number of scatters
c     chic2 : characteristic scattering angle (squared)
c.    *                                                                *
c.    * thred(na)=reduced angles of moliere theory                     *
c.    *                                                                *
c.    * f0i(na),f1i(na),f2i(na)= integrale of moliere functions        *
c.    *                                                                *
c     include 'icommon.inc'
!
      implicit real*8 (a-h,o-z)
      
      parameter (eneper = 2.7182818d0) !,pi=3.141592654d0)
      dimension tint(40),arg(4),val(4),thred(40),f0i(40),f1i(40),f2i(40)
      data thred/
     +   0.00, 0.10, 0.20, 0.30
     +,  0.40, 0.50, 0.60, 0.70
     +,  0.80, 0.90, 1.00, 1.10
     +,  1.20, 1.30, 1.40, 1.50
     +,  1.60, 1.70, 1.80, 1.90
     +,  2.00, 2.20, 2.40, 2.60
     +,  2.80, 3.00, 3.20, 3.40
     +,  3.60, 3.80, 4.00, 5.00
     +,  6.00, 7.00, 8.00, 9.00
     +, 10.00,11.00,12.00,13.00/
      data f0i/
     +  0.000000e+00 ,0.995016e-02 ,0.392106e-01 ,0.860688e-01
     + ,0.147856e+00 ,0.221199e+00 ,0.302324e+00 ,0.387374e+00
     + ,0.472708e+00 ,0.555142e+00 ,0.632121e+00 ,0.701803e+00
     + ,0.763072e+00 ,0.815480e+00 ,0.859142e+00 ,0.894601e+00
     + ,0.922695e+00 ,0.944424e+00 ,0.960836e+00 ,0.972948e+00
     + ,0.981684e+00 ,0.992093e+00 ,0.996849e+00 ,0.998841e+00
     + ,0.999606e+00 ,0.999877e+00 ,0.999964e+00 ,0.999990e+00
     + ,0.999998e+00 ,0.999999e+00 ,0.100000e+01 ,0.100000e+01
     + ,0.100000e+01 ,0.100000e+01 ,0.100000e+01 ,0.100000e+01
     + ,1.,1.,1.,1./
      data f1i/
     +  0.000000e+00,0.414985e-02,0.154894e-01,0.310312e-01
     + ,0.464438e-01,0.569008e-01,0.580763e-01,0.468264e-01
     + ,0.217924e-01,-0.163419e-01,-0.651205e-01,-0.120503e+00
     + ,-0.178272e+00,-0.233580e+00,-0.282442e+00,-0.321901e+00
     + ,-0.350115e+00,-0.366534e+00,-0.371831e+00,-0.367378e+00
     + ,-0.354994e+00,-0.314803e+00,-0.266539e+00,-0.220551e+00
     + ,-0.181546e+00,-0.150427e+00,-0.126404e+00,-0.107830e+00
     + ,-0.933106e-01,-0.817375e-01,-0.723389e-01,-0.436650e-01
     + ,-0.294700e-01,-0.212940e-01,-0.161406e-01,-0.126604e-01
     + ,-0.102042e-01,-0.840465e-02,-0.704261e-02,-0.598886e-02/
      data f2i/
     +  0.0,0.121500e-01,0.454999e-01,0.913000e-01
     + ,0.137300e+00,0.171400e+00,0.183900e+00,0.170300e+00
     + ,0.132200e+00,0.763000e-01,0.126500e-01,-0.473500e-01
     + ,-0.936000e-01,-0.119750e+00,-0.123450e+00,-0.106300e+00
     + ,-0.732800e-01,-0.312400e-01,0.128450e-01,0.528800e-01
     + ,0.844100e-01,0.114710e+00,0.106200e+00,0.765830e-01
     + ,0.435800e-01,0.173950e-01,0.695001e-03,-0.809500e-02
     + ,-0.117355e-01,-0.125449e-01,-0.120280e-01,-0.686530e-02
     + ,-0.385275e-02,-0.231115e-02,-0.147056e-02,-0.982480e-03
     + ,-0.682440e-03,-0.489715e-03,-0.361190e-03,-0.272582e-03/
*
*     ------------------------------------------------------------------
*
*
      pi=4.*atan(1.0)
      twopi=2.*pi
      
* *** compute theta angle from moliere distribution
      chic  = sqrt(chic2)
      costh=1.
      sinth=0.
      th   =0.
      if(omega.le.eneper)go to 90
      cnst=log(omega)
      b=5.
c
      do 10 l=1,10
         if(abs(b).lt.1.e-10)then
            b=1.e-10
         endif
         db=-(b-log(abs(b))-cnst)/(1.-1./b)
         b=b+db
         if(abs(db).le.0.0001)go to 20
   10 continue
c
      go to 90
c
   20 continue
      if(b.le.0.)go to 90
      binv = 1./b
      tint(1) = 0.
c
      do 30 ja=2,4
         tint(ja)=f0i(ja)+(f1i(ja)+f2i(ja)*binv)*binv
   30 continue
c
      nmax = 4
   40 continue
c      xint = ran1(idum)
      xint = unirn(dummy)
c
      do 50 na=3,40
         if(na.gt.nmax) then
            tint(na)=f0i(na)+(f1i(na)+f2i(na)*binv)*binv
            nmax=na
         endif
         if(xint.le.tint(na-1)) go to 60
   50 continue
c
      if(xint.le.tint(40)) then
         na=40
         goto 60
      else
         tmp=1.-(1.-b*(1.-xint))**5
         if(tmp.le.0.)go to 40
         thri=5./tmp
         go to 80
      endif
c
   60 continue
      na = max(na-1,3)
      na3 = na-3
c
      do 70 m=1,4
         na3m=na3+m
         arg(m)=tint(na3m)
         val(m)=thred(na3m)**2
   70 continue
c

      call polint(arg,val,4,xint,thri,err)    ! num rec
c
   80 continue

      th = chic * sqrt( abs( b * thri ) )

cc      write(6,'(4e13.5)') chic,b,thri,th
cc      pause 'check!'

      if(th.gt.pi)go to 40

      sinth = sin(th)

c     test=th*(ran1(idum))**2
      test=th*(unirn(dummy))**2
      if(test.gt.sinth)go to 40

      goto 100

   90 continue
*
* *** calculate sine and cosine of a random angle between 0 and 360 deg
*
  100 phi = unirn(dummy) * twopi

c     tth = tan(th)
c     dthxz = atan( tth * cos(phi) )
c     dthyz = atan( tth * sin(phi) )

      dthxz = sinth * cos(phi) 
      dthyz = sinth * sin(phi) 

      return
      end

c  *********************************************************
      subroutine ruthscat(omega,chia2,dthxz,dthyz)
c
c  determine scattering angle from rutherford distribution
c
c  derived form geant program  gmcoul
c
c     omega : number of scatters
c     chia2 : coulomb screening angle (squared)
c
c     include 'icommon.inc'
      
      implicit real*8 (a-h,o-z)
c
c  get actual # of scatters from poisson distribution
      omega0 = 1.167 * omega
cc    write(6,*) 'call poidev '
      nsc = poidev(omega0,idum)    ! num rec
      dthxz = 0.
      dthyz = 0.
      if( nsc .le. 0 ) go to 900
c
      do i=1,nsc
c
90     continue
       rn = unirn(dummy)
c      if( rn .lt. 1e-10) go to 90
       thet = sqrt( chia2 * ((1./rn) - 1.) )

       tth = tan(thet)

       twopi=2.*4.*atan(1.0)
       phi = twopi* unirn(dummy)

c       dthxz = dthxz + atan( tth * cos(phi) )
c       dthyz = dthyz + atan( tth * sin(phi) )

        dthxz = dthxz + sin(tth) * cos(phi)
        dthyz = dthyz + sin(tth) * sin(phi)

      end do
c
900   continue
      return
      end

c  *********************************************************
      function poidev(xm,idum)

      implicit real*8 (a-h,o-z)

      integer idum
c     real poidev,xm,pi
      parameter (pi=3.141592654)
cu    uses gammln,ran1
c     real alxm,em,g,oldm,sq,t,y,gammln,ran1
      save alxm,g,oldm,sq
!$OMP THREADPRIVATE(alxm,g,oldm,sq)
      data oldm /-1./
      if (xm.lt.12.)then
        if (xm.ne.oldm) then
          oldm=xm
          g=exp(-xm)
        endif
        em=-1
        t=1.
2       em=em+1.
        t=t*unirn(dummy)
        if (t.gt.g) goto 2
      else
        if (xm.ne.oldm) then
          oldm=xm
          sq=sqrt(2.*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1.)
        endif
1       y=tan(pi*unirn(dummy))
        em=sq*y+xm
        if (em.lt.0.) goto 1
        em=int(em)
        t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
        if (unirn(dummy).gt.t) goto 1
      endif
      poidev=em
      return
      end
c  ********************************************************
      subroutine polint(xa,ya,n,x,y,dy)

      implicit real*8 (a-h,o-z)

      integer n,nmax
      real*8 dy,x,y,xa(n),ya(n)
      parameter (nmax=10)
      integer i,m,ns
      real*8 den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
c  T.Sato 2012/8/10, delete 'pause' command
          if(den.eq.0.) write(*,*) 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
c  ********************************************************
      FUNCTION gammln(xx)

      implicit real*8 (a-h,o-z)

c     REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
!$OMP THREADPRIVATE(cof,stp)
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
