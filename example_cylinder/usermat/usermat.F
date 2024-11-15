*deck,usermat      USERDISTRIB  parallel                                gal
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and 3D/1D beam.
c
c           A 3D material constitutive model can be used for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be used.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for 
c       a plasticity model, which is the same as TB, BISO,
c       for different stress states. 
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c       This routine calls four routines,
c       usermat3d.F, usermatps.F usermatbm.F and usermat1d.F, w.r.t.
c       the corresponding stress states.
c       Each routine can be also a usermat routine for the specific 
c       element.
c
c*************************************************************************
c Copyright ANSYS.  All Rights Reserved.
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,i)               Number of state variables
c      nProp     (int,sc,i)               Number of material constants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(ncomp),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variables
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      epsPl   (dp,ar(ncomp),io)          plastic strain
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,o)                loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),o)    material jacobian matrix
c      tsstif   (dp,ar(2),o)              transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresses and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
      EXTERNAL         usermat3d, usermatps, usermatbm, usermat1d

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8
c
c*************************************************************************
c
      IF(ncomp .GE. 4) THEN
c ***    3d, plane strain and axisymmetric example
         call usermat3d (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      ELSE IF(nDirect.eq. 2 .and. ncomp .EQ. 3) THEN
c ***    plane stress example
         call usermatps (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      ELSE IF(ncomp .EQ. 3) THEN
c ***    3d beam example
         call usermatbm (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      ELSE IF(ncomp .EQ. 1) THEN
c ***    1d beam example
         call usermat1d (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      END IF
      return
      end
*deck,usermat1d    USERDISTRIB  parallel                                gal
      subroutine usermat1d(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of 1D truss element (LINK180). 
c       The plasticity model is the same as TB, BISO.
c
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variables
c            ustatev(1)                     - equivalent plastic strain
c            ustatev(2) - ustatev(1+ncomp)   - plastic strain vector
c            ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          mcomp
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, SMALL
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 mcomp      = 1
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      sigDev   (dp,ar(6  ),l)            deviatoric stress tensor
c      dfds     (dp,ar(6  ),l)            derivative of the yield function 
c      JM       (dp,ar(6,6),l)            2D matrix for a 4 order tensor
c      pEl      (dp,sc     ,l)            hydrostatic pressure stress
c      qEl      (dp,sc     ,l)            von-mises stress
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c
c
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8

      DOUBLE PRECISION qEl,   pleq_t,  sigy_t , sigy,
     &                 dpleq, pleq,    signTens,
     &                 young, posn,    sigy0,   dsigdep, 
     &                 twoG,  fratio
c*************************************************************************
c
      keycut   = 0
      dsigdep  = ZERO 
      pleq_t   = ustatev(1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
      young    = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)
c *** calculate plastic slope
      dsigdep  = young*prop(4)/(young-prop(4))
      twoG     = young / (ONE+posn)
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = HALF * twoG
c
c *** calculate elastic stiffness matrix 
c
      dsdeEl(1,1)= young
c
c *** calculate the trial stress and 
c     copy elastic moduli dsdeEl to material Jacobian matrix
      sigElp(1)   = stress(1)
      dsdePl(1,1) = dsdeEl(1,1)
      sigElp(1)   = sigElp(1) + dsdeEl(1,1) * dStrain(1)
c *** sign of predicted stress
      signTens = sign (ONE, sigElp(1))
c *** compute von-mises equivalent stress
      qEl = abs(sigElp(1))
c *** compute current yield stress
      sigy    = sigy0 + dsigdep * pleq
c
      fratio = qEl / sigy - ONE
c *** check for yielding
      IF (sigy .LE. ZERO.or.fratio .LE. -SMALL) GO TO 500
c
      sigy_t   = sigy
c *** initial guess of incremental equivalent plastic strain   
      dpleq    = (qEl - sigy) / young
      pleq     = pleq_t + dpleq
      sigy     = sigy0 + dsigdep * pleq
c
c ***  update plastic strains, stresses
      epsPl(1) = epsPl(1) + dpleq * signTens
      stress(1) =  signTens * sigy
c
c ***  update plastic strains
      epseq  = pleq
c *** Update state variables
      ustatev(1) = pleq
      ustatev(2) = epsPl(1)
c *** Update plastic work
      sedPl = sedPl + HALF * (sigy_t + sigy) * dpleq
c
c *** Material Jcobian matrix
c
      dsdePl(1,1) = dsdeEl(1,1) * dsigdep /(dsdeEl(1,1) + dsigdep)
c *** Allow a small number for Jcobian matrix if it is an ideal plasticity
      if(dsdePl(1,1).LE.ZERO) dsdePl(1,1) = SMALL*dsdeEl(1,1)
c
      goto 600
  500 continue

c *** Update stress in case of elastic/unloading
      stress(1) = sigElp(1)

  600 continue
c *** elastic strain energy
      sedEl = HALF * stress(1) * (Strain(1)+dStrain(1)-epsPl(1))
c *** update state variables
      ustatev(nStatev) = sigy
c
      return
      end
*deck,usermatbm    USERDISTRIB  parallel                                gal
      subroutine usermatbm(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for 
c       a plasticity model in 3D beam(188, 189). The plasticity
c       model is the same as TB, BISO.
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nStatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nStatev),io)       statev
c           ustatev(1)                     - equivalent plastic strain
c           ustatev(2) - ustatev(1+ncomp)  - plastic strain vector
c           ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D
c      ncomp   4   for plane strain, axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3D beam (nShear = 2), beam188/189
c      ncomp   1   for 1D beam, link180
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for Plane strain and axisymmetry
c          11, 22, 12                for Plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain, axisymmetric
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam plasticity
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), sigi(ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          NEWTON, mcomp
      DOUBLE PRECISION HALF, ONE, TWO, SMALL, SQTWOTHIRD,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 NEWTON     = 20,
     &                 mcomp      = 3
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(3  ),l)            trial stress
c      dsdeEl   (dp,ar(3,3),l)            elastic moduli
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      gamma    (dp,sc     ,l)            variable for solving incremental equivalent plastic strain
c      dgamma   (dp,sc     ,l)            correction of gamma
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      funcf    (dp,sc     ,l)            nonlinear function to be solved for gamma
c      dFdep    (dp,sc     ,l)            derivative of nonlinear function over gamma
c
c --- temperary variables for solution purpose
c      i, j
c      c1, c2, c3, fratio
c      wk1(3), wk2(3), wk3(3), wk4(3) vector working arrays
c
      EXTERNAL         vmove, vzero, vapb1, vamb1,get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), 
     &                 wk1(3), wk2(3), wk3(3), wk4(3)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8

      INTEGER          i, j, k
      DOUBLE PRECISION pleq_t,  sigy_t , sigy,
     &                 cpleq, dpleq,   pleq,    twoG,    et,
     &                 young, posn,    sigy0,   dsigdep, 
     &                 gamma, dgamma,  dfdga,   dplga,   fratio,
     &                 funcFb,funcFb2, funcf,   dFdep,
     &                 c1, c2, c3, c4, c5
      DOUBLE PRECISION pv(3)
      data pv/TWOTHIRD, TWO, TWO/
c*************************************************************************
c
      keycut   = 0
c *** equivalent plastic strain at beginning of time step
      pleq_t   = ustatev(ncomp+1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and slope of stress-strain
      young    = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)
      et       = prop(4)
c *** calculate plastic slope
      dsigdep  = young * et/(young - et)
      twoG     = young / (ONE+posn)
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = HALF * twoG
c
c *** calculate elastic stiffness matrix
c
      call vzero(dsdeEl(1,1), ncomp * ncomp)
      c1 = twoG * HALF
      dsdeEl (1,1) = young
      dsdeEl (2,2) = c1
      dsdeEl (3,3) = c1
      DO i = 1, ncomp
         wk3(i) = dsdeEl(i,i)
      END DO
c *** calculate predicted strain 
      call vmove(Strain(1), wk1(1), ncomp)
      call vapb1(wk1(1), dStrain(1), ncomp)
      call vamb1(wk1(1), ustatev(1), ncomp)

c
c *** get initial stress
      call vzero(sigi(1),ncomp)
      i = ncomp
      call get_ElmData ('ISIG', elemId,kDomIntPt, i, sigi)

c
c *** calculate the trial stress and 
c     copy elastic moduli dsdeEl to material Jacobian matrix
      call vmove(dsdeEl(1,1), dsdePl(1,1), ncomp * ncomp)
      do i=1,ncomp
         sigElp(i) = wk3(i) * wk1(i) + sigi(i)
      end do
c
      funcFb2 = ZERO
      DO i = 1, ncomp
        funcFb2 = funcFb2 + pv(i) * sigElp(i) * sigElp(i)
      END DO
      funcFb = sqrt(funcFb2)

c *** compute current yield stress
      sigy    = sigy0 + dsigdep * pleq
c
      fratio = funcFb/sigy - SQTWOTHIRD
c *** check for yielding
      IF (fratio .LE. -SMALL) GO TO 500
      sigy_t  = sigy

      DO i = 1, ncomp
         wk3(i) = wk3(i) * pv(i)
      END DO

      gamma    = ZERO
      dplga    = ZERO
      dfdga    = ZERO
      dpleq    = ZERO
      pleq     = pleq_t 

c *** Local New-Raphson procedure for solving the gamma and 
c     thus the incremental equivalent plastic strain
      DO k=1,NEWTON
         funcFb2 = ZERO
         dfdga   = ZERO
         DO j = 1 , ncomp
            c1 = ONE + gamma * wk3(j)
            c1 = ONE / c1
            c2 = sigElp(j) * c1
            wk4(j) = c2
            funcFb2 = funcFb2 + pv(j) * c2 * c2
            c2 = c2 * c2 * c1 * wk3(j) * pv(j)
            dfdga   = dfdga - c2
         END DO
         funcFb   = sqrt(funcFb2)
c ***    derivative of funcFb w.r.t. gamma
         dfdga   = dfdga / funcFb

c ***    calculate the incremental equivalent plastic strain
         dpleq    = gamma * SQTWOTHIRD * funcFb
c ***    update the total equivalent plastic strain
         pleq     = pleq_t + dpleq
c ***    current yield stress
         sigy     = sigy0 + dsigdep * pleq
c ***    calculate the residual
         funcf    = funcFb - SQTWOTHIRD * sigy
c ***    derivative of incremental equivalent plastic strain w.r.t. gamma
         dplga    = SQTWOTHIRD * (gamma * dfdga + funcFb)
c ***    derivative of residual function w.r.t. gamma
         dFdep    = dfdga - SQTWOTHIRD * dsigdep * dplga
c ***    correction of gamma
         dgamma   = -funcf / dFdep 
         gamma    = gamma   + dgamma
c ***    check for negative gamma
         gamma    = max (gamma, sqTiny)
         fratio   = funcf/ sigy
c
c ***    Check for convergence of local New-Raphson iteration
         IF (((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (abs(dgamma) .LT. ONEDM02*gamma)) .OR.
     &       ((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (dgamma      .LE. sqTiny  ) .AND.
     &        ( gamma      .LE. sqTiny  )))  GO TO 100

      END DO
c
c *** Uncovergence, set keycut to 1 for bisection/cutback the load increment
      keycut   = 1
      GO TO 990
 100  CONTINUE
c
c *** update stresses
      call vmove(wk4(1), stress(1), ncomp)

c *** calculate incremental plastic strain 
      DO j = 1, ncomp
         wk2(j) = gamma * pv(j) * wk4(j)
      END DO
c *** update plastic strains
      call vapb1(epsPl(1),wk2(1),ncomp)

c *** Update state variables
      ustatev(ncomp+1) = pleq
      do i=1,ncomp
         ustatev(i) = epsPl(i)
      end do

c *** update plastic work
      sedPl     = sedPl + HALF * (sigy_t+sigy) * dpleq

      c1     = TWOTHIRD * dsigdep
      c3     = c1 * funcFb2 / (ONE - c1 * gamma)
      DO j = 1 , ncomp
         c1 = ONE / (ONE + gamma * wk3(j))
         wk3(j) = wk3(j) * c1 / pv(j)
      END DO
      DO j = 1 , ncomp
         wk4(j) = wk4(j) * pv(j)
      END DO
      DO j = 1 , ncomp
         c3 = c3 + wk4(j) * wk4(j) * wk3(j)
      END DO
      DO j = 1 , ncomp
         wk4(j) = wk4(j) * wk3(j)
      END DO

      c3     = ONE / c3
      DO i=1,ncomp
         dsdePl(i,i) = wk3(i)
      END DO

c *** Calculate the plastic Jacobian

      DO i=1,ncomp
         DO j=1,ncomp
            dsdePl(i,j) =    dsdePl(i,j) - c3 * wk4(i) * wk4(j)
         END DO
      END DO

      goto 600

  500 continue

c *** Update stress in case of elastic/unloading
      call vmove(sigElp(1),stress(1),ncomp)

  600 continue
c *** elastic strain energy
      sedEl = ZERO
      DO i = 1 , ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i)-epsPl(i))
      END DO
      sedEl = sedEl * HALF
      ustatev(nStatev) = funcFb / SQTWOTHIRD

 990  CONTINUE
c
      return
      end
*deck,usermatps    USERDISTRIB  parallel                                gal
      subroutine usermatps(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, 
     &                   var1, var2, var3, var4, var5, 
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of plane stress state (such as PLANE182,
c       PLANE183 or SHELL181). The plasticity model is the same 
c       as ANSYS TB, BISO.
c
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nStatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nStatev),io)      user state variables
c            ustatev(1)                     - equivalent plastic strain
c            ustatev(2) - ustatev(1+ncomp)   - plastic strain vector
c            ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress, 
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), sigi(ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          NEWTON, mcomp
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, 
     &                 SQTWOTHIRD, SQTWO1,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 SQTWO1     = 0.707106769084930420d0,
     &                 NEWTON     = 20,
     &                 mcomp      = 6
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      gamma    (dp,sc     ,l)            variable for solving incremental equivalent plastic strain
c      dgamma   (dp,sc     ,l)            correction of gamma
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      funcf    (dp,sc     ,l)            nonlinear function to be solved for dpleq
c      dFdep    (dp,sc     ,l)            derivative of nonlinear function over dpleq
c
c --- temperary variables for solution purpose
c      i, j
c      con1eOv2qEl, oneOv3G, qElOv3G, con1, con2, fratio
c
      EXTERNAL         vmove, vzero, vapb1, get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), 
     &                 wk1(3), wk2(3), wk3(3), wk4(3)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8

      INTEGER          i, j, k
      DOUBLE PRECISION pleq_t,  sigy_t , sigy,
     &                 dpleq,   pleq,    twoG,    et,
     &                 young, posn,    sigy0,   dsigdep, tEo1pm,
     &                 gamma, dgamma,  dfdga,   dplga, 
     &                 funcFb,funcFb2, funcf,   dFdep,   fratio,
     &                 con1,  con2,    con3,  con4,
     &                 con2p1, ocon2p1,
     &                 ocon2p2, con4p1, ocon4p1, ocon4p2,
     &                 c1, c2, c3,c4, c5,dperr(3)
c*************************************************************************
c
      keycut   = 0
c *** equivalent plastic strain at beginning of time step
      pleq_t   = ustatev(1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and slope of stress-strain
      young    = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)
      et       = prop(4)
c *** calculate plastic slope
      dsigdep  = young * et/(young - et)
      twoG     = young / (ONE+posn)
      tEo1pm   = THIRD * young /(ONE - posn)
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = HALF * twoG
c
c *** calculate elastic stiffness matrix (3d)
c
      c1 = ONE - posn * posn
      c2 = young / c1
      c3 = posn * c2
      dsdeEl (1,1) = c2
      dsdeEl (1,2) = c3
      dsdeEl (1,3) = ZERO
      dsdeEl (2,2) = c2
      dsdeEl (2,3) = ZERO
      dsdeEl (3,3) = HALF * twoG
      do i=1,ncomp-1
        do j=i+1,ncomp
          dsdeEl(j,i)=dsdeEl(i,j)
        end do
      end do
      call vmove(ustatev(2), epsPl(1), ncomp)
c *** calculate elastic strain
      do i=1,ncomp
         wk1(i) = Strain(i) - epsPl(i) + dStrain(i)
      end do

c
c *** get initial stress
      call vzero(sigi(1),ncomp)
      i = ncomp
      call get_ElmData ('ISIG', elemId,kDomIntPt, i, sigi)

c
c *** calculate the trial stress and 
c     copy elastic moduli dsdeEl to material Jacobian matrix
      do i=1,ncomp
         sigElp(i) = ZERO
         do j=1,ncomp
            dsdePl(j,i) = dsdeEl(j,i)
            sigElp(i) = sigElp(i) + dsdeEl(j,i) * wk1(j)
         end do
         sigElp(i) = sigElp(i) + sigi(i)
      end do
c
      wk1(1)   = SQTWO1 * ( sigElp(1) +  sigElp(2) )
      wk1(2)   = SQTWO1 * (-sigElp(1) +  sigElp(2) )
      wk1(3)   = sigElp(3)
c
      funcFb2 =  THIRD * wk1(1)  * wk1(1) +
     &   wk1(2)  * wk1(2) + TWO *  wk1(3)  * wk1(3)

      funcFb  = sqrt(funcFb2)

c *** compute current yield stress
      sigy    = sigy0 + dsigdep * pleq
c
      fratio = funcFb/sigy - SQTWOTHIRD
c *** check for yielding
      IF (fratio .LE. -SMALL) GO TO 500
      sigy_t  = sigy

      gamma    = ZERO
      dplga    = ZERO
      dfdga    = ZERO
      con1     = THIRD * wk1(1) * wk1(1)
      con2     = wk1(2) * wk1(2) + TWO * wk1(3) * wk1(3)
      con3     = - con1 * tEo1pm
      con4     = - con2 * twoG
      funcf    = funcFb - SQTWOTHIRD * sigy
      dfdga    =  (con3 + con4) / funcFb
      dFdep    = dfdga - TWOTHIRD * dsigdep * funcFb
      dgamma   = -funcf / dFdep
      gamma    = gamma   + dgamma
      gamma    = max (gamma, sqTiny)

      DO k=1,NEWTON

         con2p1   = ONE + tEo1pm * gamma
         con4p1   = ONE + twoG * gamma
         ocon2p1  = ONE / con2p1
         ocon2p2  = ocon2p1 * ocon2p1
         ocon4p1  = ONE / con4p1
         ocon4p2  = ocon4p1 * ocon4p1
         funcFb2  =  con1 * ocon2p2 + con2 * ocon4p2
         funcFb   = sqrt(funcFb2)
c
c ***    calculate the incremental equivalent plastic strain
         dpleq    = gamma * SQTWOTHIRD * funcFb
c ***    update the total equivalent plastic strain
         pleq     = pleq_t + dpleq
c ***    current yield stress
         sigy     = sigy0 + dsigdep * pleq

c ***    calculate the residual
         funcf    = funcFb - SQTWOTHIRD * sigy
c ***    derivative of funcFb w.r.t. gamma
         dfdga    =  ( con3 * ocon2p2 * ocon2p1
     &               + con4 * ocon4p2 * ocon4p1) / funcFb
c ***    derivative of incremental equivalent plastic strain w.r.t. gamma
         dplga    = SQTWOTHIRD * (gamma * dfdga + funcFb)
c ***    derivative of residual function w.r.t. gamma
         dFdep  = dfdga - SQTWOTHIRD * dsigdep * dplga
         fratio   = funcf/ sigy

         dgamma   = -funcf / dFdep
         gamma    = gamma   + dgamma
         gamma    = max (gamma, sqTiny)
c
         IF (((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (abs(dgamma) .LT. ONEDM02*gamma)) .OR.
     &       ((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (dgamma      .LE. sqTiny  ) .AND.
     &        ( gamma      .LE. sqTiny  )))  GO TO 100

      END DO
c
c *** Uncovergence, set keycut to 1 for bisect/cut
      keycut   = 1
      GO TO 990
 100  CONTINUE
c
c *** update stresses
      c1 = ONE / ( ONE + gamma * tEo1pm)
      c2 = ONE / ( ONE + gamma * twoG  )
      stress(1) = SQTWO1 * (c1 * wk1(1) - c2 * wk1(2))
      stress(2) = SQTWO1 * (c1 * wk1(1) + c2 * wk1(2))
      stress(3) = c2 * wk1(3)
c
      wk1(1)    = SQTWO1 * ( stress(1) +  stress(2))
      wk1(2)    = SQTWO1 * (-stress(1) +  stress(2))
      wk1(3)    = stress(3)
      con1  = SQTWO1 * gamma
c *** calculate incremental plastic strain, wk2(i)
      wk2(1)    = con1 * ( THIRD * wk1(1) -  wk1(2))
      wk2(2)    = con1 * ( THIRD * wk1(1) +  wk1(2))
      wk2(3)    = TWO * gamma * wk1 (3)

c ***  update plastic strains
      call vapb1(epsPl(1),wk2(1),ncomp)
      epseq     = pleq
      ustatev(1) = pleq
c *** calculate plastic work
      sedPl     = sedPl + HALF * (sigy_t+sigy) * dpleq

c *** consistent tangent stiffness matrix
      con3      = TWOTHIRD*dsigdep*funcFb2/(ONE-TWOTHIRD*dsigdep*gamma)
      call vmove(wk1(1),wk4(1),3)
      c1 = wk4(1)
      c2 = wk4(2)
      c3 = wk4(3)
      c4 = tEo1pm / (ONE + tEo1pm * gamma)
      c5 = twoG / (ONE + twoG * gamma)
      wk4(1) = SQTWO1 * ( c4 * c1 - c5 * c2)
      wk4(2) = SQTWO1 * ( c4 * c1 + c5 * c2)
      wk4(3) = c5 * c3

      call vmove(wk4(1),wk3(1),3)
      c1 = wk3(1)
      c2 = wk3(2)
      wk3(1) = SQTWO1 * ( c1+c2)
      wk3(2) = SQTWO1 * (-c1+c2)

      con3     = con3 + THIRD*wk1(1)*wk3(1)+
     &           wk1(2)*wk3(2)+TWO*wk1(3)*wk3(3)
      c1  = 3.0d0 * tEo1pm / (ONE + tEo1pm * gamma)
      c2  = twoG / (ONE + twoG * gamma)
      call vzero(dsdePl(1,1),9)
      dsdePl(1,1) = HALF * (c1 + c2)
      dsdePl(2,1) = HALF * (c1 - c2)
      dsdePl(1,2) = dsdePl(2,1)
      dsdePl(2,2) = dsdePl(1,1)
      dsdePl(3,3) = HALF * c2
      con3     = ONE / con3
      DO i=1,3
         DO j=1,3
            dsdePl(i,j) = dsdePl(i,j) - con3 * wk4(i) * wk4(j)
         END DO
      END DO

      goto 600

  500 continue

c *** Update stress in case of elastic/unloading
      call vmove(sigElp(1),stress(1),ncomp)

  600 continue
c *** elastic strain and put into epsZZ
      epsZZ    = -posn/young * (stress(1) + stress(2))
c *** add plastic strain to total strain epsZZ
      epsZZ    = epsZZ - (epsPl(1) + epsPl(2))

c *** elastic strain energy
      sedEl = ZERO
      DO i = 1 , ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i)-epsPl(i))
      END DO
      sedEl    = sedEl * HALF
      ustatev(nStatev) = sigy
      do i=1,ncomp
         ustatev(i+1) = epsPl(i)
      end do

 990  CONTINUE
c
      return
      end

*deck,usermat3d  parallel           user                             jmgerken
      subroutine usermat3d(matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, F0, F1,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
!----------------------------------------------------------------------
!-  user material constitutive routine
!-    - 3D Neo-Hookean hyperelastic material
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,PARAMETER :: RK=KIND(1.D0) ! real kind
! external procedures
      EXTERNAL polarRU ! polar decomposition subroutine
      EXTERNAL trans_matrx_6 ! 6D voigt rotation matrix subroutine
	  DOUBLE PRECISION DoubleDot
      EXTERNAL transform_cy_cart_vector
! arguments
      INTEGER,INTENT(IN) :: matID      ! material ID number
      INTEGER,INTENT(IN) :: elemID      ! element ID number
      INTEGER,INTENT(IN) :: kDomIntPt    ! integration point
      INTEGER,INTENT(IN) :: kLayer      ! layer number
      INTEGER,INTENT(IN) :: kSectPt      ! section number
      INTEGER,INTENT(IN) :: ldstep     ! load step number
      INTEGER,INTENT(IN) :: isubst     ! substep number
      INTEGER,INTENT(OUT) :: keycut    ! time step reduction
      INTEGER,INTENT(IN) :: ncomp    ! no. stress components
      INTEGER,INTENT(IN) :: nDirect    ! no. direct stress components
      INTEGER,INTENT(IN) :: nShear    ! no. shear stress components
      INTEGER,INTENT(IN) :: nStatev    ! size of state variable array
      INTEGER,INTENT(IN) :: nProp    ! size of property array
      REAL(RK),INTENT(IN) :: Time ! time at begining of increment
      REAL(RK),INTENT(IN) :: dTime ! time increment
      REAL(RK),INTENT(IN) :: Temp ! temperature at begining of increment
      REAL(RK),INTENT(IN) :: dTemp ! temperature increment
      REAL(RK),INTENT(INOUT) :: sedEl ! stored (elastic) energy
      REAL(RK),INTENT(INOUT) :: sedPl ! dissipated (plastic) energy
      REAL(RK),INTENT(INOUT) :: epseq ! equivalent plastic strain
      REAL(RK),INTENT(INOUT) :: epsZZ ! plane stress thickness strain
      REAL(RK),INTENT(INOUT) :: stress(ncomp) !stress
      REAL(RK),INTENT(INOUT) :: ustatev(nStatev) !state variables
      REAL(RK),INTENT(OUT) :: dsdePl(ncomp,ncomp) !tangent stiffness
      REAL(RK),INTENT(IN) :: Strain(ncomp) !strain at begining of inc
      REAL(RK),INTENT(IN) :: dStrain(ncomp) !strain increment
      REAL(RK),INTENT(IN) :: epsPl(ncomp) !work space
      REAL(RK),INTENT(IN) :: prop(nProp) ! mat props, user input
      REAL(RK),INTENT(IN) :: coords(3) ! coordinates of material point
      REAL(RK),INTENT(IN) :: F1(3,3) ! deformation gradient at t+dt
      REAL(RK),INTENT(IN) :: F0(3,3) ! deformation gradient at t
      REAL(RK),INTENT(OUT) :: tsstif(2) ! transverse shear stiffness
      REAL(RK),INTENT(IN) :: var0(*),var1,var2,var3,var4 !place holder
      REAL(RK),INTENT(IN) :: var5,var6,var7,var8 !place holder
! local
      REAL(RK) :: R(3,3), Q(6,6) !rotations
      REAL(RK) :: U(3,3) !stretch
      REAL(RK) :: Bbar(6) ! vector form isochoric left cauchy green
      REAL(RK) :: Gd(3),outerGd(6),Fg(6),Fg33(3,3),Gd_cy(3)! growth tensor
      REAL(RK) :: s(6),s_e(6),cauchystress_fibre(6)
     &	  ,cauchystress_fibre33(3,3),cauchystress_fibre_rotate(3,3) ! Cauchy stress
      REAL(RK) :: M(6,6) ! material jacobian
      REAL(RK) :: G,K! shear & bulk modulus
      REAL(RK) :: I1bar3,detF,detF23 ! deformation invariants
      REAL(RK) :: C1,C2,C3 ! coefficient to tangent stiffness
	 
      REAL(RK) :: Fr_e(6),Fgr_e(3,3),Fgr_e_inv(3,3),G_e_33(3,3)
     &	  ,Ftot_e(3,3),Fel_e(3,3),detFel_e,Fel_iso(3,3)!elastic tensor
	 
      REAL(RK) :: Ftot_c(3,3),Fgr_c(3,3),G_c(6),G_c_33(3,3),Fr_c(3,3)
     &	  ,Fgr_c_inv(3,3),Cel_c(6),Fel_c(3,3),I4,Ftot_c_direction(3),
     &	  outer_direction_Ftot_c(6),
     &    Fel_c_inv(3,3),
     &	  Cel_c_tot(6),Fg_c(6)!collagene tensor
	 
      REAL(RK) :: a_0_cart(3),norme_a_0_cart,outer_a_0(6)
     &	  ,outer_a_0_33(3,3),a_0_cy(3)!fibre direction a0
      REAL(RK) :: Fgr_direction(3),norme_Fgr_direction,a_0(3) !Fgr*a0
      REAL(RK) :: outer_direction_Fgr(6),outer_direction_Fgr_33(3,3) ! Fgr*a0 outer Fgr*a0
	 
      REAL(RK) :: C_direction(3,3,3,3),tangent_stiffness(3,3,3,3)
     &	  ,CC_p(3,3,3,3),cc66(6,6),cc66_totale(6,6),cc66_jau(6,6) !piola transform, voight
	 
      REAL(RK) :: rho_elastine_initiale,parametre_rho
     &	  ,rho_elastine,fibre_stress,fibre_homeo_stress
     &	 ,rate_rho_collagene,remodeling_turnover_term !update mass density term
	 
      REAL(RK) :: velocity_stress(3,3),Cauchygreen_fois_remodeling(3,3),
     &   remodeling_velocity(3,3)
	  
	  
      REAL(RK),PARAMETER :: I(6) = (/1._RK,1._RK,1._RK,0._RK, 
     &      	                        0._RK,0._RK/) ! matrice d'identité
      REAL(RK), PARAMETER :: pi   =  3.14159265358979324_RK  
	  
      REAL(RK),PARAMETER ::Identity(3,3)=reshape((/1._RK,0._RK,0._RK,
     &                                      0._RK,1._RK,0._RK,
     &                                      0._RK,0._RK,1._RK/),(/3,3/))
	 
      REAL(RK) :: phi_e,phi_c,rho_totale_initiale,Dmax 
     &	  ,time_damage ,turnover_time,gain_parameter_collagene
     &	 ,ANG,k1c,k2c,dep_stretch_c,stretch_remodeling
     &	 ,turnover_time_collagene,rho_collagene,fibre_stress_sans_rot
     &   ,ang_axial,ang_cir,ang_diag,lamda_e_rr,lamda_e_theta,lamda_e_zz
	 
      REAL(RK) :: innerF1,innerF_el_c
     &	  ,lamda_elastique,terme_stiffness1,terme_stiffness2,terme_exp_1
     &	 ,terme_exp_2,stiffness_c, softness,rate_remodeling
     &	 ,I4_tot,theta,G_e(6),RT_cauchystress(3,3)
	 
      INTEGER :: z,n
	  
      REAL(RK) :: G_e_33_cy(3,3),dw_dlamda,dw2_dlamda2,norme_Gd,
     &     norme_a_0,lamda_m,norme_Ftot_direction,detFgr_e
     &	 ,detFtot_e,detFtot_c

      REAL(RK) :: Gd_cir(3),Gd_axial(3),Gd_c(3),Gd_a(3),norme_Gdc
     &	  ,norme_Gda
	 
      REAL(RK) :: fr_x,fr_y,fr_z,fc_x,fc_y,fc_z,fa_x,fa_y,fa_z
     &	  ,fdiagonal_plus_x,fdiagonal_plus_y,fdiagonal_plus_z
     &	 ,fdiagonal_minus_x,fdiagonal_minus_y, fdiagonal_minus_z
     &	 ,fiber_matrix(5,3),rotation_matrix(3,3)
      REAL(RK) :: Ge11,Ge12,Ge13,Ge21,Ge22,Ge23,Ge31,Ge32,Ge33	 
      REAL(RK), PARAMETER :: sigma_actmax = 0.054_RK
	  REAL(RK), PARAMETER :: lamda_max = 1.1_RK
      REAL(RK), PARAMETER :: lamda_zero_m = 1.1_RK
      REAL(RK), PARAMETER :: lamda_zero = 0.8_RK
      REAL(RK), PARAMETER :: rho_initiale = 1050._RK	
C       integer, dimension(30) :: myArray = [3640,3643,3646,3649,3652,3655
C      & ,3667,3670,3673,3676,3679,3682,3694,3697,3700,3703
C      & ,3706,3709,3721,3724,3727,3730,3733,3736
C      & ,3748,3751,3754,3757,3760,3763]	
      integer,PARAMETER :: myArray(30) =(/3642,3645,3648,3651,3654,3657
     & ,3669,3672,3675,3678,3681,3684,3696,3699,3702,3705
     & ,3708,3711,3723,3726,3729,3732,3735,3738
     & ,3750,3753,3756,3759,3762,3765/)
!----------------------------------------------------------------------
! deformation gradient determinant
      detF = F1(1,1) * (F1(2,2)*F1(3,3) - F1(2,3)*F1(3,2))
     &     + F1(1,2) * (F1(2,3)*F1(3,1) - F1(2,1)*F1(3,3))
     &     + F1(1,3) * (F1(2,1)*F1(3,2) - F1(2,2)*F1(3,1))
  
!define property parameter
      G = prop(1)
	  K = prop(2)
      rho_elastine_initiale   = prop(3)
      phi_c                   = prop(4)
      rho_totale_initiale     = prop(5)
      Dmax                    = prop(6)
      time_damage             = prop(7)
      turnover_time           = prop(8)
      gain_parameter_collagene= prop(9)
      ANG                     = prop(10)*pi/180._RK
      k1c                     = prop(11)
      k2c                     = prop(12)
      dep_stretch_c           = prop(13)
      turnover_time_collagene           = prop(14)
C 	  Ldamage                 == prop(15)
C 	  k1m                == prop(16)   
C 	  k2m                == prop(17)           
      ang_axial               = prop(18)*pi/180._RK
      ang_cir                 = prop(19)*pi/180._RK
      ang_diag                = prop(20)*pi/180._RK
	  lamda_e_rr              = prop(21)
	  lamda_e_theta           = prop(22)
	  lamda_e_zz              = prop(23)
	  fr_x                    = prop(24)
	  fr_y                    = prop(25)
	  fr_z                    = prop(26)
	  fc_x                    = prop(27)
	  fc_y                    = prop(28)
	  fc_z                    = prop(29)
	  fa_x                    = prop(30)
	  fa_y                    = prop(31)
	  fa_z                    = prop(32)
	  Ge11                    = prop(33)
	  Ge12                    = prop(34)
	  Ge13                    = prop(35)
	  Ge21                    = prop(36)
	  Ge22                    = prop(37)
	  Ge23                    = prop(38)
	  Ge31                    = prop(39)
	  Ge32                    = prop(40)
	  Ge33                    = prop(41)

C 	  G = prop(1) * ustatev(6)/rho_elastine_initiale	  
C       write(*,*),G, time, ustatev(6)


c*** define direction and update direction

	  Gd_c(1) = fc_x
      Gd_c(2) = fc_y
      Gd_c(3) = fc_z

	  Gd_a(1) = fa_x
      Gd_a(2) = fa_y
      Gd_a(3) = fa_z
	  
	  Gd(1) = fr_x
      Gd(2) = fr_y
      Gd(3) = fr_z
	  
C       if(isubst/=1)then
C 	  Gd_c = matmul(F1,Gd_c)
	  norme_Gdc = sqrt(Gd_c(1)**2+Gd_c(2)**2+Gd_c(3)**2)
	  Gd_c(:) = Gd_c(:)/norme_Gdc
	  
C       Gd_a = matmul(F1,Gd_a)
	  norme_Gda = sqrt(Gd_a(1)**2+Gd_a(2)**2+Gd_a(3)**2)
	  Gd_a(:) = Gd_a(:)/norme_Gda

C 	  Gd = matmul(F1,Gd)
	  norme_Gd = sqrt(Gd(1)**2+Gd(2)**2+Gd(3)**2)
	  Gd(:) = Gd(:)/norme_Gd
C    	  end if
	  
	  ! define five fibre angle(axial, circonférentiel,diagonal,-diagonal,circonférentiel(SMCS))
      fdiagonal_plus_x =  Gd_a(1) * cos(ang_diag) + Gd_c(1) * 
     &	  sin(ang_diag)
	  fdiagonal_plus_y =  Gd_a(2) * cos(ang_diag) + Gd_c(2) * 
     &	  sin(ang_diag)
	  fdiagonal_plus_z =  Gd_a(3) * cos(ang_diag) + Gd_c(3) * 
     &	  sin(ang_diag)
	  
      fdiagonal_minus_x =  Gd_a(1) * cos(ang_diag) - Gd_c(1) * 
     &	  sin(ang_diag)
	  fdiagonal_minus_y =  Gd_a(2) * cos(ang_diag) - Gd_c(2) * 
     &	  sin(ang_diag)
	  fdiagonal_minus_z =  Gd_a(3) * cos(ang_diag) - Gd_c(3) * 
     &	  sin(ang_diag)
	  
! fiber information
       fiber_matrix(1,1) =	 Gd_a(1)
	   fiber_matrix(1,2) =	 Gd_a(2)
	   fiber_matrix(1,3) =	 Gd_a(3)
       fiber_matrix(2,1) =	 Gd_c(1)
	   fiber_matrix(2,2) =	 Gd_c(2)
	   fiber_matrix(2,3) =	 Gd_c(3)
       fiber_matrix(3,1) =	 fdiagonal_plus_x
	   fiber_matrix(3,2) =	 fdiagonal_plus_y
	   fiber_matrix(3,3) =	 fdiagonal_plus_z	   
       fiber_matrix(4,1) =	 fdiagonal_minus_x
	   fiber_matrix(4,2) =	 fdiagonal_minus_y
	   fiber_matrix(4,3) =	 fdiagonal_minus_z		   
       fiber_matrix(5,1) =	 Gd_c(1)
	   fiber_matrix(5,2) =	 Gd_c(2)
	   fiber_matrix(5,3) =	 Gd_c(3)
	   
! rotation matrix
       rotation_matrix(1,1) =	 Gd(1)
	   rotation_matrix(1,2) =	 Gd(2)
	   rotation_matrix(1,3) =	 Gd(3)
       rotation_matrix(2,1) =	 Gd_c(1)
	   rotation_matrix(2,2) =	 Gd_c(2)
	   rotation_matrix(2,3) =	 Gd_c(3)
       rotation_matrix(3,1) =	 Gd_a(1)
	   rotation_matrix(3,2) =	 Gd_a(2)
	   rotation_matrix(3,3) =	 Gd_a(3)
	   
	   G_e_33(1,1) = Ge11
	   G_e_33(1,2) = Ge12
	   G_e_33(1,3) = Ge13
	   G_e_33(2,1) = Ge21
	   G_e_33(2,2) = Ge22
	   G_e_33(2,3) = Ge23
	   G_e_33(3,1) = Ge31
	   G_e_33(3,2) = Ge32
	   G_e_33(3,3) = Ge33
      
C ! deposition stretch elastin in local coordinate	  
C       G_e_33_cy(1,1)=lamda_e_rr
C 	  G_e_33_cy(1,2)=0._RK
C 	  G_e_33_cy(1,3)=0._RK
C 	  G_e_33_cy(2,1)=0._RK
C 	  G_e_33_cy(2,2)=lamda_e_theta
C 	  G_e_33_cy(2,3)=0._RK
C 	  G_e_33_cy(3,1)=0._RK
C 	  G_e_33_cy(3,2)=0._RK
C 	  G_e_33_cy(3,3)=lamda_e_zz

C ! deposition stretch elastin in global coordinate	  
C       G_e_33 = matmul(transpose(rotation_matrix),
C      &	  matmul(G_e_33_cy,rotation_matrix))
	  


c*** outer product growth direction
      call outer_product(outerGd,Gd,Gd)

! update rho_totale(t)/rho_initiale(t0)(exclude SMCS)
      if(isubst==1)then
	  ustatev(12) = 1._RK
      else
      ustatev(12) = (ustatev(6)+ustatev(7)+ustatev(8)
     &	  +ustatev(9)+ustatev(10))/rho_totale_initiale 	
      end if

!calculate normalized collagen mass density
      ustatev(18) = (ustatev(7)+ustatev(8)
     &	  +ustatev(9)+ustatev(10))/651._RK
c*** deformation matrice growth ,identique for all the constituent

      Fg(:) = I(:) + ustatev(12)*outerGd(:) - outerGd(:)
C       Fg(:) = I(:)
      ustatev(19) = sqrt(coords(1)**2._RK + coords(2)**2._RK)	  
c*******************************************************************************************     
c***************part elastine

c*** we suppose for the elastine, Fr_e=I
      Fr_e(:) = I(:)
c*** deformation Fgr_elastine
 
      call Deformation_growth_remodeling(Fgr_e,Fr_e,Fg)
      call inverse(Fgr_e_inv(:,:),Fgr_e(:,:),1._RK)
	  
c*** deformation Ftot_e
C       if(isubst==1)then
C 	  Ftot_e(:,:) = matmul(Identity,G_e_33_cy)  
C 	  else
C 	  Ftot_e(:,:) = matmul(F1,G_e_33_cy) 
C 	  end if
	  Ftot_e(:,:) = matmul(F1,G_e_33)  
C       Ftot_e(:,:) = F1(:,:)
		!save Ftot_e in state variable
	  ustatev(21) = Ftot_e(1,1)
	  ustatev(22) = Ftot_e(1,2)
	  ustatev(23) = Ftot_e(1,3)
	  ustatev(24) = Ftot_e(2,1)
	  ustatev(25) = Ftot_e(2,2)
	  ustatev(26) = Ftot_e(2,3)
	  ustatev(27) = Ftot_e(3,1)
	  ustatev(28) = Ftot_e(3,2)
	  ustatev(29) = Ftot_e(3,3)					
c*** deformation Fel_elastine
      Fel_e(:,:) = matmul(Ftot_e,Fgr_e_inv)  


c*** deformation gradient élastique élastine determinant
      
      call determinant(Fel_e,detFel_e)
	  
c***isochoric Fel
      Fel_iso(:,:) = detFel_e**(-1._RK/3._RK)*Fel_e(:,:)
	  
! isochoric left Cauchy Green

      call Left_cauchygreen(Bbar,Fel_iso)
      
      I1bar3 = (Bbar(1) + Bbar(2) + Bbar(3))/3._RK ! 1st invariant / 3
 
!Cauchy stress
      s_e(1:3) = (G*(Bbar(1:3)-I1bar3)/ detFel_e + K*(detFel_e-1._RK))
     &	  * ustatev(6)
      s_e(4:6) = (G*(Bbar(4:6)) / detFel_e)
     &	  * ustatev(6)

!-------------------------------------------------------------------------------	  
	  !material Jacobian, spatial frame
      C1 = 2._RK/3._RK*G/detFel_e*ustatev(6)
      C2 = 3._RK/4._RK*C1	  
	  M(:,:) = 0._RK    
      ! shear modulus terms
      M(1,1) = C1 * (Bbar(1) + I1bar3)
      M(2,1) =-C1 * (Bbar(1) + Bbar(2) - I1bar3)
      M(3,1) =-C1 * (Bbar(1) + Bbar(3) - I1bar3)
      M(4,1) = C1 * Bbar(4) / 2._RK
      M(5,1) =-C1 * Bbar(5)
      M(6,1) = C1 * Bbar(6) / 2._RK

      M(1,2) = M(2,1)
      M(2,2) = C1 * (Bbar(2) + I1bar3)
      M(3,2) =-C1 * (Bbar(2) + Bbar(3) - I1bar3)
      M(4,2) = C1 * Bbar(4) / 2._RK
      M(5,2) = C1 * Bbar(5) / 2._RK
      M(6,2) =-C1 * Bbar(6)

      M(1,3) = M(3,1)
      M(2,3) = M(3,2)
      M(3,3) = C1 * (Bbar(3) + I1bar3)
      M(4,3) =-C1 * Bbar(4)
      M(5,3) = C1 * Bbar(5) / 2._RK
      M(6,3) = C1 * Bbar(6) / 2._RK

      M(1,4) = M(4,1)
      M(2,4) = M(4,2)
      M(3,4) = M(4,3)
      M(4,4) = C2 * (Bbar(1) + Bbar(2))
      M(5,4) = C2 * Bbar(6)
      M(6,4) = C2 * Bbar(5)

      M(1,5) = M(5,1)
      M(2,5) = M(5,2)
      M(3,5) = M(5,3)
      M(4,5) = M(5,4)
      M(5,5) = C2 * (Bbar(2) + Bbar(3))
      M(6,5) = C2 * Bbar(4)

      M(1,6) = M(6,1)
      M(2,6) = M(6,2)
      M(3,6) = M(6,3)
      M(4,6) = M(6,4)
      M(5,6) = M(6,5)
      M(6,6) = C2 * (Bbar(1) + Bbar(3))
      
      ! bulk modulus terms
      M(1:3,1:3) = M(1:3,1:3) + K*(2._RK*detFel_e - 1._RK)* ustatev(6)

c*******************************************************************************************	  
c************************part collagène

!initialize the term for the tangent stiffness
      cc66_totale(:,:) = 0._RK
	  s(:) = s_e(:)
! calculate cauchy stress for five fibers
      DO z = 1,5

          if(z == 5)then
		     k1c = prop(16)
             k2c = prop(17)
          end if
		  
c*** fibre direction		  
         !change the angle of the collagene and the SMC

          a_0(1)=fiber_matrix(z,1)
          a_0(2)=fiber_matrix(z,2)
          a_0(3)=fiber_matrix(z,3)
		  
C           if(isubst/=1)then
C 	      a_0 = matmul(F1,a_0)
C 	      norme_a_0 = sqrt(a_0(1)**2+a_0(2)**2+a_0(3)**2)
C 	      a_0(:) = a_0(:)/norme_a_0
C           end if
	  
	  
C 
     	  

          call outer_product(outer_a_0,a_0,a_0)
          call change_shape_six_to_three(outer_a_0_33,outer_a_0)	  					

!!!!!!! calculate homeostatic stress
          if(isubst ==1)then
          G_c(:) = dep_stretch_c*outer_a_0(:)+(1._RK/sqrt(dep_stretch_c)
     &	       )*(I(:)-outer_a_0(:))

          call change_shape_six_to_three(G_c_33,G_c)
          Ftot_c(:,:) = matmul(Identity,G_c_33)

          Fr_c(:,:) = Identity(:,:)

          

C 	c*** deformation growth and remodeling Fgr_c

          call change_shape_six_to_three(Fg33,Fg)	
          Fgr_c(:,:)=matmul(Fr_c,Fg33)																
          call Inverse(Fgr_c_inv(:,:),Fgr_c(:,:),1._RK)


C 	c*** deformation élastique collagène
          Fel_c = matmul(Ftot_c,Fgr_c_inv)

          call Right_cauchygreen(Cel_c,Fel_c)


C 	c*** growth direction fibre

          Fgr_direction(:) = matmul(Fgr_c,a_0)

          norme_Fgr_direction = Fgr_direction(1)**2._RK
     &      +Fgr_direction(2)**2._RK
     &      +Fgr_direction(3)**2._RK
		 
          call outer_product(outer_direction_Fgr,
     &      Fgr_direction,Fgr_direction)
	      outer_direction_Fgr(:) = outer_direction_Fgr(:)
     &		  /norme_Fgr_direction

		  
          I4 = Cel_c(1)*outer_direction_Fgr(1) 
     &   +Cel_c(2)*outer_direction_Fgr(2)
     &   +Cel_c(3)*outer_direction_Fgr(3)
     &   +2._RK*Cel_c(4)*outer_direction_Fgr(4)
     &   +2._RK*Cel_c(5)*outer_direction_Fgr(5)
     &   +2._RK*Cel_c(6)*outer_direction_Fgr(6)  
	      
	 
          call Right_cauchygreen(Cel_c_tot,Identity)
		  
          I4_tot = Cel_c_tot(1)*outer_a_0(1) 
     &   + Cel_c_tot(2)*outer_a_0(2)
     &   + Cel_c_tot(3)*outer_a_0(3)
     &   + 2._RK*Cel_c_tot(4)*outer_a_0(4)
     &   + 2._RK*Cel_c_tot(5)*outer_a_0(5)
     &   + 2._RK*Cel_c_tot(6)*outer_a_0(6) 
          if (I4-1._RK<0._RK .AND. z/=5) then
		  k1c = k1c*0.075_RK
          endif
    
          Ftot_c_direction(:) = matmul(Ftot_c,a_0)
		  
		  norme_Ftot_direction = Ftot_c_direction(1)**2._RK
     &      +Ftot_c_direction(2)**2._RK
     &      +Ftot_c_direction(3)**2._RK
	 
	      lamda_m =sqrt(norme_Ftot_direction)
          call outer_product(outer_direction_Ftot_c,
     &      Ftot_c_direction,Ftot_c_direction)
	 
	      if (z/=5) then
          cauchystress_fibre(:) = (2._RK *k1c*exp(k2c*(I4-1._RK)**2._RK)
     &	  *(I4-1._RK)*outer_direction_Ftot_c(:)/norme_Fgr_direction)
     &	 *ustatev(6+z)/detF
	      else
          cauchystress_fibre(:) =((2._RK *k1c*exp(k2c*(I4-1._RK)**2._RK)
     &	  *(I4-1._RK)*outer_direction_Ftot_c(:)/norme_Fgr_direction)
     &	+sigma_actmax/(rho_initiale*I4_tot)*(1-(lamda_max-lamda_m
     &   *lamda_zero_m)**2/(lamda_max-lamda_zero)**2)
     &   *outer_direction_Ftot_c(:)/norme_Ftot_direction)
     &	 *ustatev(6+z)/detF
          endif   
	 
          call change_shape_six_to_three(cauchystress_fibre33
     &	  ,cauchystress_fibre)	  
          call polarRU(Identity,R,U) ! rotation from polar decomposition
		 
          RT_cauchystress(:,:) = matmul(transpose(R)
     &		  ,cauchystress_fibre33)
          cauchystress_fibre_rotate(:,:) = matmul(RT_cauchystress,R)
          
          fibre_stress = DoubleDot(outer_a_0_33
     &		  ,cauchystress_fibre_rotate)
        ! homeostatic stress = cauchy stress when time=0

          ustatev(12+z) = fibre_stress
C 		  if(elemID==1)then
C 		  write(*,*)'homeostress',elemId,  ustatev(13)
C           endif
          endif
    
	 
	 

C 	c*** deposition stretch collagène
	      G_c(:) = dep_stretch_c*outer_a_0(:)+(1._RK/sqrt(dep_stretch_c))
     &	       *(I(:)-outer_a_0(:))

          call change_shape_six_to_three(G_c_33,G_c)
          Ftot_c(:,:) = matmul(F1,G_c_33)
          if (isubst == 1) then
		  Fr_c(:,:) = Identity(:,:)
		  else
          Fr_c(:,:) = ustatev(z) * outer_a_0_33(:,:) 
     &    + 1._RK/SQRT(ustatev(z))*(Identity(:,:)-outer_a_0_33(:,:))
	      end if

          

C 	c*** deformation growth and remodeling Fgr_c
C           Fg(:) = I(:)
          call change_shape_six_to_three(Fg33,Fg)	
          Fgr_c(:,:)=matmul(Fr_c,Fg33)																
          call Inverse(Fgr_c_inv(:,:),Fgr_c(:,:),1._RK)


C 	c*** deformation élastique collagène
          Fel_c = matmul(Ftot_c,Fgr_c_inv)

          call Right_cauchygreen(Cel_c,Fel_c)


C 	c*** growth direction fibre

          Fgr_direction(:) = matmul(Fgr_c,a_0)

          norme_Fgr_direction = Fgr_direction(1)**2._RK
     &      +Fgr_direction(2)**2._RK
     &      +Fgr_direction(3)**2._RK
		 
          call outer_product(outer_direction_Fgr,
     &      Fgr_direction,Fgr_direction)
	      outer_direction_Fgr(:) = outer_direction_Fgr(:)
     &		  /norme_Fgr_direction

		  
           I4 = Cel_c(1)*outer_direction_Fgr(1) 
     &   +Cel_c(2)*outer_direction_Fgr(2)
     &   +Cel_c(3)*outer_direction_Fgr(3)
     &   +2._RK*Cel_c(4)*outer_direction_Fgr(4)
     &   +2._RK*Cel_c(5)*outer_direction_Fgr(5)
     &   +2._RK*Cel_c(6)*outer_direction_Fgr(6) 

          call Right_cauchygreen(Cel_c_tot,Ftot_c)
		  
          I4_tot = Cel_c_tot(1)*outer_a_0(1) 
     &   + Cel_c_tot(2)*outer_a_0(2)
     &   + Cel_c_tot(3)*outer_a_0(3)
     &   + 2._RK*Cel_c_tot(4)*outer_a_0(4)
     &   + 2._RK*Cel_c_tot(5)*outer_a_0(5)
     &   + 2._RK*Cel_c_tot(6)*outer_a_0(6)  	 
	 

         Ftot_c_direction(:) = matmul(Ftot_c,a_0)
		  
		  norme_Ftot_direction = Ftot_c_direction(1)**2._RK
     &      +Ftot_c_direction(2)**2._RK
     &      +Ftot_c_direction(3)**2._RK
	 
	      lamda_m =sqrt(norme_Ftot_direction)
		 
          call outer_product(outer_direction_Ftot_c,
     &      Ftot_c_direction,Ftot_c_direction)
	 
															
							  
	 
	      if (z/=5) then
          cauchystress_fibre(:) = (2._RK *k1c*exp(k2c*(I4-1._RK)**2._RK)
     &	  *(I4-1._RK)*outer_direction_Ftot_c(:)/norme_Fgr_direction)
     &	 *ustatev(6+z)/detF
	      else
          cauchystress_fibre(:) =((2._RK *k1c*exp(k2c*(I4-1._RK)**2._RK)
     &	  *(I4-1._RK)*outer_direction_Ftot_c(:)/norme_Fgr_direction)
     &	+sigma_actmax/(rho_initiale*I4_tot)*(1-(lamda_max-lamda_m
     &   *lamda_zero_m)**2/(lamda_max-lamda_zero)**2)
     &   *outer_direction_Ftot_c(:)/norme_Ftot_direction)
     &	 *ustatev(6+z)/detF
          endif

          !update cauchy stress
          
          s(:) = s(:) + cauchystress_fibre(:)
		  

 !-------------------------------------------------------------------------------  
! calculate tangent stiffness for fiber 

	!material Jacobian, spatial frame
	      if (z/=5) then
          C3 = (8._RK*k1c*k2c*(I4-1._RK)**2+4._RK*k1c)
     &	  *exp(k2c*(I4-1._RK)**2._RK)*
     &       (ustatev(6+z)/detF)/(norme_Fgr_direction**2)
	      else
          C3 = ((8._RK*k1c*k2c*(I4-1._RK)**2+4._RK*k1c)
     &	  *exp(k2c*(I4-1._RK)**2._RK)/(norme_Fgr_direction**2)-
     &	 2._RK*(sigma_actmax/(rho_initiale*I4_tot**2))
     &	 *(1-(lamda_max-lamda_m*lamda_zero_m)**2
     &	 /(lamda_max-lamda_zero)**2))
     &       *(ustatev(6+z))/detF
          endif				 
	!stiffness_direction
	      call TensorProduct(C_direction,outer_a_0_33
     &	  ,outer_a_0_33,1._RK)  
          tangent_stiffness(:,:,:,:)=C3*C_direction(:,:,:,:)
  
	! piola transform
          call piola_transorm(CC_p,F1,tangent_stiffness,detF) ! revoir use the F1 or Fel_c
          call Voight(cc66,CC_p) 
 	
	!jausmann tangant stiffness
          cc66(1,1)= cc66(1,1) +2._RK*cauchystress_fibre(1)
          cc66(1,4)= cc66(1,4) +cauchystress_fibre(4)
          cc66(1,6)= cc66(1,6) +cauchystress_fibre(6)
          
          cc66(2,2)= cc66(2,2) +2._RK*cauchystress_fibre(2)
          cc66(2,4)= cc66(2,4) +cauchystress_fibre(4)
          cc66(2,5)= cc66(2,5) +cauchystress_fibre(5)
          
          cc66(3,3)= cc66(3,3) +2._RK*cauchystress_fibre(3)
          cc66(3,5)= cc66(3,5) +cauchystress_fibre(5)
          cc66(3,6)= cc66(3,6) +cauchystress_fibre(6) 
          
          cc66(4,1)= cc66(4,1) +cauchystress_fibre(4)
          cc66(4,2)= cc66(4,2) +cauchystress_fibre(4)
          cc66(4,4)= cc66(4,4) +(cauchystress_fibre(2)
     &	  +cauchystress_fibre(1))/2._RK
          cc66(4,5)= cc66(4,5) +cauchystress_fibre(6)/2._RK
          cc66(4,6)= cc66(4,6) +cauchystress_fibre(5)/2._RK
          
          cc66(5,2)= cc66(5,2) +cauchystress_fibre(5)
          cc66(5,3)= cc66(5,3) +cauchystress_fibre(5)
          cc66(5,4)= cc66(5,4) +cauchystress_fibre(6)/2._RK
          cc66(5,5)= cc66(5,5) +(cauchystress_fibre(2)
     &	  +cauchystress_fibre(3))/2._RK
          cc66(5,6)= cc66(5,6) +cauchystress_fibre(4)/2._RK

          cc66(6,1)= cc66(6,1) +cauchystress_fibre(6)
          cc66(6,3)= cc66(6,3) +cauchystress_fibre(6)
          cc66(6,4)= cc66(6,4) +cauchystress_fibre(5)/2._RK
          cc66(6,5)= cc66(6,5) +cauchystress_fibre(4)/2._RK
          cc66(6,6)= cc66(6,6) +(cauchystress_fibre(1)
     &	  +cauchystress_fibre(3))/2._RK	  
  
          cc66_totale(:,:)= cc66_totale(:,:) + cc66(:,:)
!-------------------------------------------------------------------------------
	!calculate rate mass density collagene
          call change_shape_six_to_three(cauchystress_fibre33
     &	  ,cauchystress_fibre)	  
          call polarRU(F1,R,U) ! rotation from polar decomposition
		 
          RT_cauchystress(:,:) = matmul(transpose(R)
     &		  ,cauchystress_fibre33)
          cauchystress_fibre_rotate(:,:) = matmul(RT_cauchystress,R)
          
          fibre_stress = DoubleDot(outer_a_0_33
     &		  ,cauchystress_fibre_rotate)

			 
          fibre_homeo_stress = ustatev(12+z)

		 
C 		   !rate rho_collagene
          if (z/=5) then
          rate_rho_collagene = gain_parameter_collagene * ustatev(6+z)
     &    * (fibre_stress-fibre_homeo_stress) / fibre_homeo_stress
	      else if(z==5) then
		  rate_rho_collagene = 0._RK
          end if
		  
c*******************************************************************************************	  
c************************part update 
	  
C 	!calculate rate remodeling of collagene	  

          remodeling_turnover_term =  rate_rho_collagene/ustatev(6+z)
     &	  + 1._RK/turnover_time_collagene


          innerF1= I4_tot

          innerF_el_c = Cel_c(1)*outer_a_0(1) 
     &   + Cel_c(2)*outer_a_0(2)
     &   + Cel_c(3)*outer_a_0(3)
     &   + 2._RK*Cel_c(4)*outer_a_0(4)
     &   + 2._RK*Cel_c(5)*outer_a_0(5)
     &   + 2._RK*Cel_c(6)*outer_a_0(6)  

          lamda_elastique = sqrt(innerF_el_c)

! dérivée de strain energie dw/dlamda +lamda*dw2/dlamda2

C           ! dw/dlamda
C 		  dw_dlamda = 2._RK*k1c*(lamda_elastique**3-lamda_elastique)
C      &		  *exp(k2c*(lamda_elastique**2._RK-1._RK)**2._RK)
	 
C 	      ! dw2/dlamda2
C 		  dw2_dlamda2 = (2._RK*k1c*(3._RK*lamda_elastique**2-1._RK) 
C      &	        + 8._RK*k1c*k2c*lamda_elastique*(lamda_elastique**3
C      &	        -lamda_elastique)*(innerF_el_c-1._RK))
C      &	        *exp(k2c*(lamda_elastique**2._RK-1._RK)**2._RK)
	           ! dw/dlamda
		  dw_dlamda = 2._RK*k1c*(lamda_elastique**3-lamda_elastique)
     &		  *exp(k2c*(lamda_elastique**2._RK-1._RK)**2._RK)
	 
	      ! dw2/dlamda2
		  dw2_dlamda2 = (2._RK*k1c*(3._RK*lamda_elastique**2-1._RK) 
     &	        + 8._RK*k1c*k2c*(lamda_elastique**3
     &	        -lamda_elastique)**2)
     &	        *exp(k2c*(lamda_elastique**2._RK-1._RK)**2._RK)
	 
	      stiffness_c = (dw_dlamda + lamda_elastique*dw2_dlamda2)
     &		  *ustatev(6+z)/detF
C 	      stiffness_c = (4._RK*k2c*innerF_el_c*(innerF_el_c-1._RK)**2._RK
C      &		  +4._RK*innerF_el_c-2._RK)*k1c*lamda_elastique
C      &        *EXP(k2c*(innerF_el_c-1._RK)**2._RK)*ustatev(6+z)
		  
          softness = sqrt(innerF1)/(innerF_el_c*stiffness_c)


          rate_remodeling = remodeling_turnover_term
     &      *softness*(fibre_stress-fibre_homeo_stress)
          
C           write(*,*),z,cc66(2,2),cc66(3,3),cc66(5,5)
!update remodeling and mass density
C           if (isubst/=1) then	
          stretch_remodeling = ustatev(z) + dTime * rate_remodeling
          ustatev(z) = stretch_remodeling
		  
C           end if
		  
C           if ( isubst/=1) then				  
          rho_collagene = ustatev(6+z) + dTime * rate_rho_collagene 
          ustatev(6+z) = rho_collagene

	 
C           end if


      end do	  
	    
	  
!-------------------------------------------------------------------------------	  
	  	  !update stiffness
		  
      M(:,:) = M(:,:)+cc66_totale(:,:)

!-------------------------------------------------------------------------------
C !update mass density elastin with locale in intervalle 0-height/2

C       call damage_locale_elastin(parametre_rho,time
C      &	  ,turnover_time,Dmax,time_damage,coords(3),prop(15))
C       rho_elastine = rho_elastine_initiale * parametre_rho
C 	  ustatev(6) = rho_elastine
	  
      if (isubst==1) then	
	  parametre_rho = 0
      else
C       parametre_rho = -ustatev(6)/turnover_time-Dmax/time_damage
C      &	  *exp(-time/time_damage)
C      &	 *rho_elastine_initiale
	 
      parametre_rho = -Dmax/time_damage
     &	  *exp(-time/time_damage-0.5*(coords(3)/prop(15))**2)
     &	 *rho_elastine_initiale

      end if
      rho_elastine = ustatev(6) +dTime * parametre_rho
	  ustatev(6) = rho_elastine
	  
 
C       do n=1,30
C 	     if(elemID == myArray(n)) then 
C 		    parametre_rho = -ustatev(6)/turnover_time-Dmax/time_damage
C      &	  *exp(-time/time_damage)
C      &	 *rho_elastine_initiale
C 	        rho_elastine = ustatev(6) +dTime * parametre_rho
C 	        ustatev(6) = rho_elastine	  
C C 	     write(*,*),elemID
C 	     else
C 		    parametre_rho = 0
C 			rho_elastine = ustatev(6) +dTime * parametre_rho
C 	        ustatev(6) = rho_elastine	  
C          endif
C       end do


	  

c*******************************************************************************************	  
c************************part output
!Change from spatial to corotated frame
      call polarRU(F1,R,U) ! rotation from polar decomposition
      call trans_matrx_6(Q,transpose(R)) ! dim 6 rotation matrix
      stress(:) = matmul(Q,s)
      dsdePl(:,:) = matmul(Q,matmul(M,transpose(Q)))
      
C !stored and dissipated energy
C       sedEl = G/2._RK*(I1bar3*3._RK-3._RK) + K/2._RK*(detF-1._RK)**2
C       sedEl = sedEl/detF
C       sedPl = 0._RK
      
! transverse shear stiffness
      tsstif(:) = dsdePl(6,6)
      
      end subroutine usermat3d
!-------------------------------------------------------------------------------
      subroutine polarRU(F,R,U)
!----------------------------------------------------------------------
!-  Polar decomposition, F=RU
!    Reference:
!    JC Simo & TJR Hughes, Computational Inelasticity, 1998, Springer
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,PARAMETER :: RK=KIND(1.D0) ! real kind
! arguments
      REAL(RK),INTENT(IN) :: F(3,3) ! matrix to decompose
      REAL(RK),INTENT(OUT) :: R(3,3) !rotation
      REAL(RK),INTENT(OUT) :: U(3,3) !stretch
! local
      REAL(RK), PARAMETER :: pi   =  3.14159265358979324_RK
      REAL(RK),PARAMETER ::I(3,3)=reshape((/1._RK,0._RK,0._RK,
     &                                      0._RK,1._RK,0._RK,
     &                                      0._RK,0._RK,1._RK/),(/3,3/))
      REAL(RK) :: C(3,3),CC(3,3) ! right Cauchy-Green, squared
      REAL(RK) :: Ui(3,3) ! U inverse
      REAL(RK) :: IC, IIC, IIIC ! invariants of C
      REAL(RK) :: IU, IIU, IIIU ! invariants of U
      REAL(RK) :: l1,l2,l3 ! eigenvalues 
      REAL(RK) :: p,q,m,n,t,D ! constants
!----------------------------------------------------------------------
!right Cauchy-Green
      C(:,:) = matmul(transpose(F),F)
      CC(:,:) = matmul(C,C)
      
      !invariants
      IC = C(1,1) + C(2,2) + C(3,3)
      IIC = 0.5_RK*(IC**2  - (CC(1,1) + CC(2,2) + CC(3,3)))
      IIIC = C(1,1) * (C(2,2)*C(3,3) - C(2,3)*C(3,2))
     &     + C(1,2) * (C(2,3)*C(3,1) - C(2,1)*C(3,3))
     &     + C(1,3) * (C(2,1)*C(3,2) - C(2,2)*C(3,1))
      
! eigenvalues of sqrt(C)
      p = IIC - (IC**2)/3._RK
      q = -(2._RK/27._RK)*IC**3+IC*IIC/3._RK-IIIC
      if(abs(p)<epsilon(1._RK))then
        l1 = sqrt( abs(-abs(q)**(1._RK/3._RK) + IC/3._RK) )
        l2 = l1
        l3 = l2
      else
        m = 2._RK*sqrt(abs(p)/3._RK)
        n = 3._RK*q/(m*p)
        if(abs(n)>1._RK) n = sign(1._RK, n)
        t = atan2(sqrt(1._RK-n**2),n)/3._RK
        l1 = sqrt( abs(m*cos(t) + IC/3._RK) )
        l2 = sqrt( abs(m*cos(t+2._RK/3._RK*pi) + IC/3._RK) )
        l3 = sqrt( abs(m*cos(t+4._RK/3._RK*pi) + IC/3._RK) )
      endif
      
!stretch and inverse
       !invariants
       IU = l1 + l2 + l3
       IIU = l1*l2 + l1*l3 + l2*l3
       IIIU = l1*l2*l3
       D = IU*IIU-IIIU
       
       U(:,:)  = (-CC(:,:) + (IU**2-IIU)*C(:,:) + IU*IIIU*I(:,:))/D
       Ui(:,:) = (C(:,:) - IU*U(:,:) + IIU*I(:,:))/IIIU

! Rotation
      R(:,:) = matmul(F,Ui)
        
      end subroutine polarRU
!----------------------------------------------------------------------
      SUBROUTINE trans_matrx_6(Q,R)
!----------------------------------------------------------------------
! Construct the 6x6 Voigt rotation matrix
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,PARAMETER :: RK=KIND(1.D0) ! real kind
! arguments
      REAL(RK), INTENT(OUT) :: Q(6,6)
      REAL(RK), INTENT(IN) :: R(3,3)
!----------------------------------------------------------------------
      Q(1,1) = R(1,1)**2
      Q(2,1) = R(2,1)**2
      Q(3,1) = R(3,1)**2
      Q(4,1) = R(1,1)*R(2,1)
      Q(5,1) = R(2,1)*R(3,1)
      Q(6,1) = R(3,1)*R(1,1)

      Q(1,2) = R(1,2)**2
      Q(2,2) = R(2,2)**2
      Q(3,2) = R(3,2)**2
      Q(4,2) = R(1,2)*R(2,2)
      Q(5,2) = R(2,2)*R(3,2)
      Q(6,2) = R(3,2)*R(1,2)

      Q(1,3) = R(1,3)**2
      Q(2,3) = R(2,3)**2
      Q(3,3) = R(3,3)**2
      Q(4,3) = R(1,3)*R(2,3)
      Q(5,3) = R(2,3)*R(3,3)
      Q(6,3) = R(3,3)*R(1,3)

      Q(1,4) = R(1,1)*R(1,2)*2._RK
      Q(2,4) = R(2,1)*R(2,2)*2._RK
      Q(3,4) = R(3,1)*R(3,2)*2._RK
      Q(4,4) = R(1,1)*R(2,2)+R(1,2)*R(2,1)
      Q(5,4) = R(2,1)*R(3,2)+R(2,2)*R(3,1)
      Q(6,4) = R(3,1)*R(1,2)+R(3,2)*R(1,1)

      Q(1,5) = R(1,2)*R(1,3)*2._RK
      Q(2,5) = R(2,2)*R(2,3)*2._RK
      Q(3,5) = R(3,2)*R(3,3)*2._RK
      Q(4,5) = R(1,2)*R(2,3)+R(1,3)*R(2,2)
      Q(5,5) = R(2,2)*R(3,3)+R(2,3)*R(3,2)
      Q(6,5) = R(3,2)*R(1,3)+R(3,3)*R(1,2)

      Q(1,6) = R(1,1)*R(1,3)*2._RK
      Q(2,6) = R(2,1)*R(2,3)*2._RK
      Q(3,6) = R(3,1)*R(3,3)*2._RK
      Q(4,6) = R(1,1)*R(2,3)+R(1,3)*R(2,1)
      Q(5,6) = R(2,1)*R(3,3)+R(2,3)*R(3,1)
      Q(6,6) = R(3,1)*R(1,3)+R(3,3)*R(1,1)
      END SUBROUTINE trans_matrx_6
	  
!-------------------------------------------------------------------------------	  
	  ! Returns Inverse of the matrix
      SUBROUTINE Inverse(TensRes,Tensr,fact)
      DOUBLE PRECISION TensRes(3,3),Tensr(3,3),fact,d
      INTEGER i,j
      
      TensRes(1,1)=Tensr(2,2)*Tensr(3,3)-Tensr(2,3)*Tensr(3,2)
      TensRes(1,2)=Tensr(1,3)*Tensr(3,2)-Tensr(3,3)*Tensr(1,2)
      TensRes(1,3)=Tensr(1,2)*Tensr(2,3)-Tensr(2,2)*Tensr(1,3)
      TensRes(2,1)=Tensr(2,3)*Tensr(3,1)-Tensr(3,3)*Tensr(2,1)
      TensRes(2,2)=Tensr(1,1)*Tensr(3,3)-Tensr(1,3)*Tensr(3,1)
      TensRes(2,3)=Tensr(1,3)*Tensr(2,1)-Tensr(2,3)*Tensr(1,1)
      TensRes(3,1)=Tensr(2,1)*Tensr(3,2)-Tensr(2,2)*Tensr(3,1)
      TensRes(3,2)=Tensr(1,2)*Tensr(3,1)-Tensr(1,1)*Tensr(3,2)
      TensRes(3,3)=Tensr(2,2)*Tensr(1,1)-Tensr(2,1)*Tensr(1,2)
      
c***      d=determinant(Tensr)
      call determinant(Tensr,d)
      DO j=1,3
		DO i=1,3
			TensRes(i,j)=TensRes(i,j)*fact/d
		END DO
      END DO
      
      RETURN
      END SUBROUTINE Inverse
      
!-------------------------------------------------------------------------------
      subroutine determinant(T,detT)
      DOUBLE PRECISION detT,T(3,3)
      detT = T(1,1) * (T(2,2)*T(3,3) - T(2,3)*T(3,2))
     &     + T(1,2) * (T(2,3)*T(3,1) - T(2,1)*T(3,3))
     &     + T(1,3) * (T(2,1)*T(3,2) - T(2,2)*T(3,1))
      RETURN
      END 
!-------------------------------------------------------------------------------
      SUBROUTINE Right_cauchygreen(Cbar,F)
      DOUBLE PRECISION Cbar(6),F(3,3)
      integer i
      
      Cbar(1) = F(1,1)*F(1,1)+F(2,1)*F(2,1)
     &      +F(3,1)*F(3,1)
      Cbar(2) = F(1,2)*F(1,2)+F(2,2)*F(2,2)
     &      +F(3,2)*F(3,2)
      Cbar(3) = F(1,3)*F(1,3)+F(2,3)*F(2,3)
     &      +F(3,3)*F(3,3)
      Cbar(4) = F(1,2)*F(1,1)+F(2,2)*F(2,1)
     &      +F(3,2)*F(3,1)
      Cbar(5) = F(1,3)*F(1,2)+F(2,3)*F(2,2)
     &      +F(3,3)*F(3,2)
      Cbar(6) = F(1,3)*F(1,1)+F(2,3)*F(2,1)
     &      +F(3,3)*F(3,1)
      
      RETURN
      END SUBROUTINE Right_cauchygreen
	  
      SUBROUTINE outer_product(outer,vectora,vectorb)

      DOUBLE PRECISION outer(6),vectora(3),vectorb(3)
      
      outer(1)= vectora(1)*vectorb(1)
      outer(2)= vectora(2)*vectorb(2)
      outer(3)= vectora(3)*vectorb(3)
      outer(4)= vectora(1)*vectorb(2)
      outer(5)= vectora(2)*vectorb(3)
      outer(6)= vectora(3)*vectorb(1)     

      RETURN 
      END SUBROUTINE outer_product
      SUBROUTINE Left_cauchygreen(Bbar,F)
      DOUBLE PRECISION Bbar(6),F(3,3)
      
      Bbar(1) = F(1,1)*F(1,1)+F(1,2)*F(1,2)
     &      +F(1,3)*F(1,3)
      Bbar(2) = F(2,1)*F(2,1)+F(2,2)*F(2,2)
     &      +F(2,3)*F(2,3)
      Bbar(3) = F(3,1)*F(3,1)+F(3,2)*F(3,2)
     &      +F(3,3)*F(3,3)
      Bbar(4) = F(2,1)*F(1,1)+F(2,2)*F(1,2)
     &      +F(2,3)*F(1,3)
      Bbar(5) = F(3,1)*F(2,1)+F(3,2)*F(2,2)
     &      +F(3,3)*F(2,3)
      Bbar(6) = F(3,1)*F(1,1)+F(3,2)*F(1,2)
     &      +F(3,3)*F(1,3)           
      RETURN
      END SUBROUTINE Left_cauchygreen
	  
      SUBROUTINE change_shape_six_to_three(vector33,vector6)

      DOUBLE PRECISION vector33(3,3),vector6(6)
      
      vector33(1,1) = vector6(1)
      vector33(2,2) = vector6(2)
      vector33(3,3) = vector6(3)
      vector33(1,2) = vector6(4)
      vector33(2,3) = vector6(5)
      vector33(1,3) = vector6(6)
      vector33(2,1) = vector6(4)
      vector33(3,2) = vector6(5)
      vector33(3,1) = vector6(6)     
	  
      RETURN
      END SUBROUTINE change_shape_six_to_three
	  
      SUBROUTINE DoubleDot4(C,A,B,fact)
      DOUBLE PRECISION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3) 
      DOUBLE PRECISION fact 
      INTEGER i,j,k,l,m,n 
      do i=1,3 
        do j=1,3 
          do k=1,3 
            do l=1,3 
               C(i,j,k,l)=0.0d0 
               do m=1,3 
                 do n=1,3 
                    C(i,j,k,l)=C(i,j,k,l)+fact*A(i,j,m,n)*B(m,n,k,l) 
                 end do 
               end do 
            end do 
          end do 
        end do 
      end do 
      END
	
      SUBROUTINE TensorProduct(TensRes,Tens1,Tens2,fact) 
      DOUBLE PRECISION TensRes(3,3,3,3),Tens1(3,3),Tens2(3,3),fact 
      INTEGER i,j,k,l 
      do i=1,3 
         do j=1,3 
            do k=1,3 
               do l=1,3 
                  TensRes(i,j,k,l)=fact*Tens1(i,j)*Tens2(k,l) 
               end do 
            end do 
         end do 
      end do 
      END
	  
      SUBROUTINE Voight(Tss,Ttttt) 
      Double Precision Tss(6,6),Ttttt(3,3,3,3) 
      Tss(1,1)=Ttttt(1,1,1,1) 
      Tss(1,2)=Ttttt(1,1,2,2) 
      Tss(1,3)=Ttttt(1,1,3,3) 
      Tss(1,4)=Ttttt(1,1,1,2) 
      Tss(1,5)=Ttttt(1,1,2,3) 
      Tss(1,6)=Ttttt(1,1,1,3) 
      Tss(2,1)=Ttttt(2,2,1,1) 
      Tss(2,2)=Ttttt(2,2,2,2) 
      Tss(2,3)=Ttttt(2,2,3,3) 
      Tss(2,4)=Ttttt(2,2,1,2) 
      Tss(2,5)=Ttttt(2,2,2,3) 
      Tss(2,6)=Ttttt(2,2,1,3) 
      Tss(3,1)=Ttttt(3,3,1,1) 
      Tss(3,2)=Ttttt(3,3,2,2) 
      Tss(3,3)=Ttttt(3,3,3,3) 
      Tss(3,4)=Ttttt(3,3,1,2) 
      Tss(3,5)=Ttttt(3,3,2,3) 
      Tss(3,6)=Ttttt(3,3,1,3) 
      Tss(4,1)=Ttttt(1,2,1,1) 
      Tss(4,2)=Ttttt(1,2,2,2) 
      Tss(4,3)=Ttttt(1,2,3,3) 
      Tss(4,4)=Ttttt(1,2,1,2) 
      Tss(4,5)=Ttttt(1,2,2,3) 
      Tss(4,6)=Ttttt(1,2,1,3) 
      Tss(5,1)=Ttttt(2,3,1,1) 
      Tss(5,2)=Ttttt(2,3,2,2) 
      Tss(5,3)=Ttttt(2,3,3,3) 
      Tss(5,4)=Ttttt(2,3,1,2) 
      Tss(5,5)=Ttttt(2,3,2,3) 
      Tss(5,6)=Ttttt(2,3,1,3) 
      Tss(6,1)=Ttttt(1,3,1,1) 
      Tss(6,2)=Ttttt(1,3,2,2) 
      Tss(6,3)=Ttttt(1,3,3,3) 
      Tss(6,4)=Ttttt(1,3,1,2) 
      Tss(6,5)=Ttttt(1,3,2,3) 
      Tss(6,6)=Ttttt(1,3,1,3) 
      end 
	  
      subroutine piola_transorm(CC,F,stiffness,fact)
      Double Precision CC(3,3,3,3),F(3,3),stiffness(3,3,3,3),fact
      INTEGER i,j,k,l,ic,jc,kc,lc 
      do i=1,3 
         do j=1,3 
            do k=1,3 
               do l=1,3 
                  CC(i,j,k,l)=0.0d0
                  do ic=1,3 
                     do jc=1,3 
                        do kc=1,3 
                           do lc=1,3 
                              CC(i,j,k,l)=CC(i,j,k,l)+ 
     & F(i,ic)*F(j,jc)*F(k,kc)*F(l,lc)*stiffness(ic,jc,kc,lc)/fact

                           end do 
                        end do 
                     end do 
                  end do 
               end do 
            end do 
         end do 
      end do 
      end	  
	  
      subroutine gauss(x,b,a,n)
      Double Precision b(9),a(9,9),n,x(9),cte
      INTEGER i,j,k
      x(:) = 0
      cte = 0
      DO i=1,n-1

         DO j =i+2,n
            cte = -a(j,i)/a(i,i)
			write(*,*),'cte',cte
            DO k = 1,n
                a(j,k) = a(j,k) + a(i,k)*cte
            end do
            b(j) = b(j) + b(i)*cte
         end do
      end do 
      x(n-1)=b(n-1)/a(n-1,n-1)
	  write(*,*),'x(n-1)',x(n-1)
      DO i =n-1,1,-1
         DO j =i+1,n
            b(i) = b(i) - a(i,j)*x(j)
         end do
         x(i) = b(i)/a(i,i)
         
      end do 
      end
	  
      subroutine remodeling_term_coeff(x,c,Aa)
      Double Precision a(9,9),b(9),c(3,3),Aa(3,3,3,3),a1(9),x(3,3)
      INTEGER i,j,k,l
	  b = (/c(1,1),c(1,2),c(1,3),c(2,1),c(2,2),c(2,3),c(3,1)
     &	  ,c(3,2),c(3,3) /)
	  write(*,*),'b',b
      a(:,:)=0
      DO i =1,3
           DO j  =1,3
              a1(:)=0
              DO k =1,3
                  DO l =1,3
                     a1((k-1)*3+l) = Aa(i,j,k,l)
				  end do
              end do
              a((i-1)*3+j,:) = a1(:)
           end do
      end do
      write(*,*),'a',a
	  
      call gauss(x,b,a,9)
      x = reshape( x, (/ 3, 3 /), order = (/ 2, 1 /))
	  write(*,*),'x',x
      end
	  
      Function DoubleDot(Tens1,Tens2) 
      DOUBLE PRECISION Tens1(3,3),Tens2(3,3),DoubleDot 
      INTEGER i,j 
      DoubleDot=0.0d0 
      do i=1,3 
         do j=1,3 
            DoubleDot=DoubleDot+Tens1(i,j)*Tens2(i,j) 
         end do 
      end do 
      return 
      END 	

      SUBROUTINE Deformation_growth_remodeling(Fgr,Fg,Fr)
      DOUBLE PRECISION Fgr(3,3),Fg(6),Fr(6)
      
      Fgr(1,1) = Fg(1) * Fr(1) + Fg(4) * Fr(4) + Fg(6) * Fr(6)
      Fgr(2,2) = Fg(4) * Fr(4) + Fg(2) * Fr(2) + Fg(5) * Fr(5)
      Fgr(3,3) = Fg(6) * Fr(6) + Fg(5) * Fr(5) + Fg(3) * Fr(3)
      Fgr(1,2) = Fg(1) * Fr(4) + Fg(4) * Fr(2) + Fg(6) * Fr(5)
      Fgr(2,3) = Fg(4) * Fr(6) + Fg(2) * Fr(5) + Fg(5) * Fr(3)
      Fgr(1,3) = Fg(1) * Fr(6) + Fg(4) * Fr(5) + Fg(6) * Fr(3)
      Fgr(2,1) = Fg(4) * Fr(1) + Fg(2) * Fr(4) + Fg(5) * Fr(6)
      Fgr(3,2) = Fg(6) * Fr(4) + Fg(5) * Fr(2) + Fg(3) * Fr(5)
      Fgr(3,1) = Fg(6) * Fr(1) + Fg(5) * Fr(4) + Fg(3) * Fr(6)           
      RETURN
      END SUBROUTINE Deformation_growth_remodeling	  
	  
      SUBROUTINE transform_cy_cart_vector(Gd,Gd_cy,theta)
      DOUBLE PRECISION Gd(3),Gd_cy(3),theta
      
	  Gd(1) = COS(theta)*Gd_cy(1)-SIN(theta)*Gd_cy(2)
	  Gd(2) = SIN(theta)*Gd_cy(1)+COS(theta)*Gd_cy(2)
	  Gd(3) = Gd_cy(3)     
      RETURN
      END SUBROUTINE transform_cy_cart_vector	
	  
      SUBROUTINE transform_cy_cart_matrix(G_e_33,G_e_33_cy,theta)
      DOUBLE PRECISION G_e_33(3,3),G_e_33_cy(3,3),theta,transform_a(3,3)
     &	  ,transform_b(3,3),G_e_33_transform(3,3)
	 
	  transform_a(1,1)=COS(theta)
	  transform_a(1,2)=-SIN(theta)
	  transform_a(1,3)=0
	  transform_a(2,1)=SIN(theta)
	  transform_a(2,2)=COS(theta)
	  transform_a(2,3)=0
	  transform_a(3,1)=0
	  transform_a(3,2)=0
	  transform_a(3,3)=1	  
	  
      transform_b(1,1)=COS(theta)
	  transform_b(1,2)=SIN(theta)
	  transform_b(1,3)=0
	  transform_b(2,1)=-SIN(theta)
	  transform_b(2,2)=COS(theta)
	  transform_b(2,3)=0
	  transform_b(3,1)=0
	  transform_b(3,2)=0
	  transform_b(3,3)=1	  
      
      G_e_33_transform = matmul(transform_a,G_e_33_cy)
	  G_e_33 = matmul(G_e_33_transform,transform_b)
        
      RETURN
      END SUBROUTINE transform_cy_cart_matrix	
	  
      SUBROUTINE damage_locale_elastin(parametre_rho,time
     &	  ,turnover_time,Dmax,time_damage,z,L_damage)
      DOUBLE PRECISION parametre_rho,time,turnover_time
     &	  ,Dmax,time_damage,z,L_damage
      
C       parametre_rho = exp(-time/turnover_time) + Dmax 
C      &	  * (turnover_time/(time_damage-turnover_time) ) 
C      &	 * (exp(-time/turnover_time)-exp(-time/time_damage))*
C      &	 exp(-0.5*(z/L_damage)**2)    
	 
      parametre_rho = exp(-time/turnover_time)*(1-Dmax 
     &	  * (turnover_time/(turnover_time-time_damage))
     &	 *exp(-0.5*(z/L_damage)**2)) + Dmax 
     &	  * (turnover_time/(turnover_time-time_damage) ) 
     &	 * exp(-time/time_damage-0.5*(z/L_damage)**2)
      RETURN
      END SUBROUTINE damage_locale_elastin	
	  
	  