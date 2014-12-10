************************************************************************
*                                                                      *
      subroutine usrtally(ncol)
*                                                                      *
*        sample subroutine for user defined tally.                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*       ncol = 1 : start of calculation                                *
*              2 : end of calculation                                  *
*              3 : end of a batch                                      *
*              4 : source                                              *
*              5 : detection of geometry error                         *
*              6 : recovery of geometry error                          *
*              7 : termination by geometry error                       *
*              8 : termination by weight cut-off                       *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : (n,x) reaction                                      *
*             14 : (n,n'x) reaction                                    *
*             15 : sequential transport only for tally                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        In the distributed-memory parallel computing,                 *
*          npe : total number of used Processor Elements               *
*          me : ID number of each processor                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        In the shared memory parallel computing,                      *
*          ipomp : ID number of each core                              *
*          npomp : total number of used cores                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        nocas : current event number in this batch                    :
*        nobch : current batch number                                  *
*        rcasc : real number of NOCAS+maxcas*(NOBCH-1)                 *
*        rsouin : sum of the weight of source particle                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        no : cascade id in this event                                 *
*        idmn(mat) : material id                                       *
*        ityp : particle type                                          *
*        ktyp : particle kf-code                                       *
*        jtyp : charge number of the particle                          *
*        mtyp : baryon number of the particle                          *
*        rtyp : rest mass of the particle (MeV)                        *
*        oldwt : weight of the particle at (x,y,z)                     *
*        qs : dE/dx of electron at (x,y,z)                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        iblz1 : cell id at (x,y,z)                                    *
*        iblz2 : cell id after crossing                                *
*        ilev1 : level structure id of the cell at (x,y,z)             *
*          ilat1(i,j)                                                  *
*        ilev2 : level structure id of the cell after crossing         *
*          ilat2(i,j)                                                  *
*        costha : cosine of theta on surface crossing                  *
*        uang(1) : x of position(?) at surface crossing                *
*        uang(2) : y of position(?) at surface crossing                *
*        uang(3) : z of position(?) at surface crossing                *
*        nsurf : surface number                                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        name(no) : collision number of the particle                   *
*        ncnt(1,no) : values of counter 1                              *
*        ncnt(2,no) : values of counter 2                              *
*        ncnt(3,no) : values of counter 3                              *
*        wt(no) : weight of the particle at (xc,yc,zc)                 *
*        u(no) : x, y, z-components of unit vector of                  *
*        v(no) :momentum of the particle                               *
*        w(no) :                                                       *
*        e(no) : energy of the particle at (x,y,z) (MeV)               *
*        t(no) : time of the particle at (x,y,z) (nsec)                *
*        x(no) : x, y, z-position coordinates of                       *
*        y(no) :the preceding event point (cm)                         *
*        z(no) :                                                       *
*        ec(no) : energy of the particle at (xc,yc,zc) (MeV)           *
*        tc(no) : time of the particle at (xc,yc,zc) (nsec)            *
*        xc(no) : x, y, z-position coordinates of                      *
*        yc(no) :the particle (cm)                                     *
*        zc(no) :                                                      *
*        spx(no) : x, y, z-components of unit vector of                *
*        spy(no) :spin direction of the particle                       *
*        spz(no) :                                                     *
*        nzst(no)                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        nclsts : the number of produced particle and nucleus          *
*        mathz : Z number of the mother nucleus                        *
*        mathn : N number of the mother nucleus                        *
*        jcoll : reaction type id1                                     *
*        kcoll : reaction type id2                                     *
*                                                                      *
*----------------------------------------------------------------------*
*        jcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        jcoll : =  0, nothing happen                                  *
*                =  1, Hydrogen collisions                             *
*                =  2, Particle Decays                                 *
*                =  3, Elastic collisions                              *
*                =  4, High Energy Nuclear collisions                  *
*                =  5, Heavy Ion reactions                             *
*                =  6, Neutron reactions by data                       *
*                =  7, Photon reactions by data                        *
*                =  8, Electron reactions by data                      *
*                =  9, Proton reactions by data                        *
*                = 10, Neutron event mode                              *
*                = 11, Delta Ray production                            *
*                = 12, Muon interaction                                *
*                = 13, Photon by EGS5                                  *
*                = 14, Electron by EGS5                                *
*                                                                      *
*----------------------------------------------------------------------*
*        kcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        kcoll : =  0, normal                                          *
*                =  1, high energy fission                             *
*                =  2, high energy absorption                          *
*                =  3, low energy n elastic                            *
*                =  4, low energy n non-elastic                        *
*                =  5, low energy n fission                            *
*                =  6, low energy n absorption                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        nclsts   : total number of out going particles and nuclei     *
*                                                                      *
*        iclusts(nclsts)                                               *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4, status of the particle 0: real, <0 : dead      *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        qclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, px (GeV/c)                                     *
*                  = 2, py (GeV/c)                                     *
*                  = 3, pz (GeV/c)                                     *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight                                         *
*                  = 9, time                                           *
*                  = 10, x                                             *
*                  = 11, y                                             *
*                  = 12, z                                             *
*                                                                      *
*        jcount(i,nclsts)                                              *
*                                                                      *
************************************************************************
*                                                                      *
*        kf code table                                                 *
*                                                                      *
*           kf-code: ityp :  description                               *
*                                                                      *
*             2212 :   1  :  proton                                    *
*             2112 :   2  :  neutron                                   *
*              211 :   3  :  pion (+)                                  *
*              111 :   4  :  pion (0)                                  *
*             -211 :   5  :  pion (-)                                  *
*              -13 :   6  :  muon (+)                                  *
*               13 :   7  :  muon (-)                                  *
*              321 :   8  :  kaon (+)                                  *
*              311 :   9  :  kaon (0)                                  *
*             -321 :  10  :  kaon (-)                                  *
*                                                                      *
*               11 :  12  :  electron                                  *
*              -11 :  13  :  positron                                  *
*               22 :  14  :  photon                                    *
*          1000002 :  15  :  deuteron                                  *
*          1000003 :  16  :  triton                                    *
*          2000003 :  17  :  3he                                       *
*          2000004 :  18  :  alpha                                     *
*      Z*1000000+A :  19  :  nucleus                                   *
*                                                                      *
*           kf-code of the other transport particles (ityp=11)         *
*               12 :         nu_e                                      *
*               14 :         nu_mu                                     *
*              221 :         eta                                       *
*              331 :         eta'                                      *
*             -311 :         k0bar                                     *
*            -2112 :         nbar                                      *
*            -2212 :         pbar                                      *
*             3122 :         Lanbda0                                   *
*             3222 :         Sigma+                                    *
*             3212 :         Sigma0                                    *
*             3112 :         Sigma-                                    *
*             3322 :         Xi0                                       *
*             3312 :         Xi-                                       *
*             3334 :         Omega-                                    *
*                                                                      *
************************************************************************

      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me
      common /ipomp0/ipomp,npomp
!$OMP THREADPRIVATE(/ipomp0/)

*-----------------------------------------------------------------------

      common /cusrtally/ iusrtally, iudtf(50)
      common /cudtpara/ udtpara(0:9)

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /rcomon/ rcasc

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)
      common /celdb/  idsn(kvlmax), idtn(kvlmax)

cKN 2012/10/31
      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, noz, mtel
!$OMP THREADPRIVATE(/elect/)

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

      common /taliin/ rsouin, nzztin, nrgnin
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      
C      common /userdefined/ io, para1
      integer pdg, evt, bch, evtcas
	  double precision x1, y1, z1, px, py, pz, ee, wi
	  common /root/ pdg, evt, bch, evtcas, x1, y1, z1, px, py, pz, ee, 
     &              wi

      io = iudtf(1)

*-----------------------------------------------------------------------
*     modified by Ye, Yang
*-----------------------------------------------------------------------

        if( iblz1 .eq. udtpara(0)) then
            pdg = ktyp
            evt = no
            bch = nobch
            evtcas = nocas
            x1 = x(no)
            y1 = y(no)
            z1 = z(no)
            px = u(no)
            py = w(no)
            pz = v(no)
            ee = e(no)
            wi = wt(no)
C            print *, id
			call rootlink()
C				write(io, 1000) ktyp, x(no), y(no), z(no), u(no), w(no), v(no),
C     &				 e(no), wt(no), t(no), ncnt(1,no), ncnt(2,no), ncnt(3,no),
C     &				 nocas, nobch, no
C1000			format(i9, 9e12.4, 3i9, ' 999  999  999  999', 3i9)
        end if

*-----------------------------------------------------------------------
*        ncol = 1, 2, 3
*-----------------------------------------------------------------------

!         if( ncol .eq. 1 .or. ncol .eq. 2 .or. ncol .eq. 3 ) then

!               write(io,1000) ncol
! 1000          format('NCOL=',/i3)
!               if (ncol .ne. 2) write(io,*)
!               return

!         end if

*-----------------------------------------------------------------------
*        ncol = 4 - 15
*-----------------------------------------------------------------------

!         if( ncol .ge. 4 ) then

!               if( ncol .eq. 4 ) write(io,'(//)')
!               write(io,1000) ncol

*-----------------------------------------------------------------------

!            if( ncol .eq. 4 ) then

!               write(io,1001) nocas, nobch, rcasc, rsouin
! 1001          format('NOCAS,NOBCH,RCASC,RSOUIN=',/2i9,1p2d16.8)

!            end if

*-----------------------------------------------------------------------

 !              write(io,1002)no,idmn(mat),ityp,ktyp,jtyp,mtyp,rtyp,oldwt
 !1002          format('NO,MAT,ITYP,KTYP,JTYP,MTYP,RTYP,OLDWT='
 !    &              ,/6i9,1p2d16.8)

!            if( ityp .eq. 12 .or. ityp .eq. 13 ) then

!               write(io,1003) qs
! 1003          format('QS=',/1p1d16.8)

!            end if

*-----------------------------------------------------------------------

!               write(io,1004) iblz1,iblz2, ilev1,ilev2
! 1004          format('IBLZ1,IBLZ2,ILEV1,ILEV2=',/4i9)

!            if( ilev1 .gt. 0 ) then

!               write(io,1005) ( ( ilat1(i,j), i=1,5 ), j=1,ilev1 )
! 1005          format('ILAT1(i,j)=',5i9)

!            end if

!            if( ilev2 .gt. 0 ) then

!               write(io,1006) ( ( ilat2(i,j), i=1,5 ), j=1,ilev2 )
! 1006          format('ILAT2(i,j)=',5i9)

!            end if

*-----------------------------------------------------------------------

!            if( ncol .eq. 10 ) then

!               write(io,1007) costha,uang(1),uang(2),uang(3),idsn(nsurf)
! 1007          format('COS,UANG(1,2,3),NSURF=',/1p4d16.8,i9)

!            end if

*-----------------------------------------------------------------------

!               write(io,1008) name(no), (ncnt(i,no),i=1,3)
! 1008          format('NAME,NCNT(1,2,3)=',/4i9)

!               write(io,1009) wt(no),
!     &                     u(no),v(no),w(no)
! 1009          format('WT,U,V,W=',/1p4d16.8)

!               write(io,1010) e(no), t(no),
!     &                     x(no),y(no),z(no)
! 1010          format('E,T,X,Y,Z=',/1p5d16.8)

!               write(io,1011) ec(no), tc(no),
!     &                     xc(no),yc(no),zc(no)
! 1011          format('EC,TC,XC,YC,ZC=',/1p5d16.8)

!               write(io,1012)
!     &              spx(no),spy(no),spz(no)
! 1012          format('SPX,SPY,SPZ=',/1p3d16.8)

!               write(io,1013) nzst(no)
! 1013          format('NZST=',/i9)

*-----------------------------------------------------------------------

!            if( ncol .eq. 13 .or. ncol .eq. 14 ) then

!                  write(io,1014) nclsts
! 1014             format('NCLSTS=',/i9)

!               if( nclsts .gt. 0 ) then

!                     write(io,1015) mathz, mathn, jcoll, kcoll
! 1015                format('MATHZ,MATHN,JCOLL,KCOLL=',4i9)

!                  do i = 1, nclsts

!                     write(io,1016) i, iclusts(i)
! 1016                format('ICLUSTS(i=',i9,')=',/i9)

!                     write(io,1017) i, ( jclusts(j,i), j=0,8)
! 1017                format('JCLUST(j=0-8,i=',i9,')=',/9i9)

!                     write(io,1018) i, ( qclusts(j,i), j=0,12)
! 1018                format('QCLUST(j=0-12,i=',i9,')='
!     &                    ,/1p5d16.8,/1p5d16.8,/1p3d16.8)

!                     write(io,1019) i, ( jcount(j,i),  j=1,3)
! 1019                format('JCOUNT(j=1-3,i=',i9,')=',/3i9)

!                  end do

!               end if

!            end if

*-----------------------------------------------------------------------

!         end if
!         write(io,*)

*-----------------------------------------------------------------------

      return
      end
