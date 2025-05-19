!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.50-g
!
! SUBROUTINE
!   AllocateMem
!
! FILE
!   AllocateMem
!
! DESCRIPTION
!   Allocation of memory for Global State variables and other Variables
!  
!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! AUTHORS
!   mfstep/ERSEM team
!
! CHANGE_LOG
!   ---
!
! COPYING
!   
!   Copyright (C) 2004 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!
  subroutine AllocateMem
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  use (import) other modules
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem
  use mem
  use mem_Param
  use bio_bfm, only: calc_sigma_depth
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  integer:: status
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Start the allocation of pelagic state global
  ! matrix and pointers
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-"
  write(LOGUNIT,*) "# Allocating State Variables and Rates array ..."

#ifndef NOT_STANDALONE
 
     allocate(D3STATE(1:NO_BOXES,1:NO_D3_BOX_STATES),stat=status)
     if (status /=0)call error_msg_prn(ALLOC,"AllocateMem","D3STATE" )
     D3STATE= ZERO
     allocate(D3SOURCE(1:NO_BOXES,1:NO_D3_BOX_STATES,1:NO_D3_BOX_STATES),stat=status)
     if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D3SOURCE")
     D3SOURCE = ZERO
     allocate(D3SINK(1:NO_BOXES,1:NO_D3_BOX_STATES,1:NO_D3_BOX_STATES)&
      & ,stat=status)
     if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D3SINK")
     D3SINK = ZERO
     allocate(D3STATETYPE(1:NO_D3_BOX_STATES ),stat=status)
     if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem","D3STATETYPE")
     D3STATETYPE = ZERO

 
     allocate(D2STATE(1:NO_BOXES_XY,1:NO_D2_BOX_STATES),stat=status)
     if (status /=0)call error_msg_prn(ALLOC,"AllocateMem","D2STATE" )
     D2STATE= ZERO
     allocate(D2SOURCE(1:NO_BOXES_XY,1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES),stat=status)
     if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D2SOURCE")
     D2SOURCE = ZERO
     allocate(D2SINK(1:NO_BOXES_XY,1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES)&
      & ,stat=status)
     if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D2SINK")
     D2SINK = ZERO
     allocate(D2STATETYPE(1:NO_D2_BOX_STATES ),stat=status)
     if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem","D2STATETYPE")
     D2STATETYPE = ZERO

 
     allocate(D3DIAGNOS(1:NO_BOXES,1:NO_D3_BOX_DIAGNOSS),stat=status)
     if (status /=0)call error_msg_prn(ALLOC,"AllocateMem","D3DIAGNOS" )
     D3DIAGNOS= ZERO

 
     allocate(D2DIAGNOS(1:NO_BOXES_XY,1:NO_D2_BOX_DIAGNOSS),stat=status)
     if (status /=0)call error_msg_prn(ALLOC,"AllocateMem","D2DIAGNOS" )
     D2DIAGNOS= ZERO

 
     allocate(D3DIAGNOS_PRF(1:40,1:NO_D3_BOX_DIAGNOSS_PRF),stat=status)
     if (status /=0)call error_msg_prn(ALLOC,"AllocateMem","D3DIAGNOS_PRF"&
      & )
     D3DIAGNOS_PRF= ZERO

#endif


    allocate(CoupledtoBDc(1:iiPhytoPlankton ),stat=status); CoupledtoBDc = 0


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocation of Pelagic variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    R9x => D3STATE(:,ppR9x); R9x=ZERO
    O2o => D3STATE(:,ppO2o); O2o=ZERO
    N1p => D3STATE(:,ppN1p); N1p=ZERO
    N3n => D3STATE(:,ppN3n); N3n=ZERO
    N4n => D3STATE(:,ppN4n); N4n=ZERO
    N5s => D3STATE(:,ppN5s); N5s=ZERO
    N6r => D3STATE(:,ppN6r); N6r=ZERO
    B1c => D3STATE(:,ppB1c); B1c=ZERO
    B1n => D3STATE(:,ppB1n); B1n=ZERO
    B1p => D3STATE(:,ppB1p); B1p=ZERO
    Bac => D3STATE(:,ppBac); Bac=ZERO
    P1c => D3STATE(:,ppP1c); P1c=ZERO
    P1n => D3STATE(:,ppP1n); P1n=ZERO
    P1p => D3STATE(:,ppP1p); P1p=ZERO
    P1l => D3STATE(:,ppP1l); P1l=ZERO
    P1s => D3STATE(:,ppP1s); P1s=ZERO
    P2c => D3STATE(:,ppP2c); P2c=ZERO
    P2n => D3STATE(:,ppP2n); P2n=ZERO
    P2p => D3STATE(:,ppP2p); P2p=ZERO
    P2l => D3STATE(:,ppP2l); P2l=ZERO
    P3c => D3STATE(:,ppP3c); P3c=ZERO
    P3n => D3STATE(:,ppP3n); P3n=ZERO
    P3p => D3STATE(:,ppP3p); P3p=ZERO
    P3l => D3STATE(:,ppP3l); P3l=ZERO
    P4c => D3STATE(:,ppP4c); P4c=ZERO
    P4n => D3STATE(:,ppP4n); P4n=ZERO
    P4p => D3STATE(:,ppP4p); P4p=ZERO
    P4l => D3STATE(:,ppP4l); P4l=ZERO
    P5c => D3STATE(:,ppP5c); P5c=ZERO
    P5n => D3STATE(:,ppP5n); P5n=ZERO
    P5p => D3STATE(:,ppP5p); P5p=ZERO
    P5l => D3STATE(:,ppP5l); P5l=ZERO
    P5s => D3STATE(:,ppP5s); P5s=ZERO
    P6c => D3STATE(:,ppP6c); P6c=ZERO
    P6n => D3STATE(:,ppP6n); P6n=ZERO
    P6p => D3STATE(:,ppP6p); P6p=ZERO
    P6l => D3STATE(:,ppP6l); P6l=ZERO
    Pcc => D3STATE(:,ppPcc); Pcc=ZERO
    R1c => D3STATE(:,ppR1c); R1c=ZERO
    R1n => D3STATE(:,ppR1n); R1n=ZERO
    R1p => D3STATE(:,ppR1p); R1p=ZERO
    R2c => D3STATE(:,ppR2c); R2c=ZERO
    R2n => D3STATE(:,ppR2n); R2n=ZERO
    R3c => D3STATE(:,ppR3c); R3c=ZERO
    R6c => D3STATE(:,ppR6c); R6c=ZERO
    R6n => D3STATE(:,ppR6n); R6n=ZERO
    R6p => D3STATE(:,ppR6p); R6p=ZERO
    R6s => D3STATE(:,ppR6s); R6s=ZERO
    RZc => D3STATE(:,ppRZc); RZc=ZERO
    O3c => D3STATE(:,ppO3c); O3c=ZERO
    O3h => D3STATE(:,ppO3h); O3h=ZERO
    Z3c => D3STATE(:,ppZ3c); Z3c=ZERO
    Z4c => D3STATE(:,ppZ4c); Z4c=ZERO
    Z2c => D3STATE(:,ppZ2c); Z2c=ZERO
    Z5c => D3STATE(:,ppZ5c); Z5c=ZERO
    Z6c => D3STATE(:,ppZ6c); Z6c=ZERO


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocation of Benthic variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    BP1c => D2STATE(:,ppBP1c); BP1c=ZERO
    BP1n => D2STATE(:,ppBP1n); BP1n=ZERO
    BP1p => D2STATE(:,ppBP1p); BP1p=ZERO
    BP1l => D2STATE(:,ppBP1l); BP1l=ZERO
    BP1s => D2STATE(:,ppBP1s); BP1s=ZERO
    Y1c => D2STATE(:,ppY1c); Y1c=ZERO
    Y1n => D2STATE(:,ppY1n); Y1n=ZERO
    Y1p => D2STATE(:,ppY1p); Y1p=ZERO
    Y2c => D2STATE(:,ppY2c); Y2c=ZERO
    Y2n => D2STATE(:,ppY2n); Y2n=ZERO
    Y2p => D2STATE(:,ppY2p); Y2p=ZERO
    Y4c => D2STATE(:,ppY4c); Y4c=ZERO
    Y4n => D2STATE(:,ppY4n); Y4n=ZERO
    Y4p => D2STATE(:,ppY4p); Y4p=ZERO
    Y5c => D2STATE(:,ppY5c); Y5c=ZERO
    Y5n => D2STATE(:,ppY5n); Y5n=ZERO
    Y5p => D2STATE(:,ppY5p); Y5p=ZERO
    Yy3c => D2STATE(:,ppYy3c); Yy3c=ZERO
    Yy3n => D2STATE(:,ppYy3n); Yy3n=ZERO
    Yy3p => D2STATE(:,ppYy3p); Yy3p=ZERO
    Y3c => D2STATE(:,ppY3c); Y3c=ZERO
    Y3n => D2STATE(:,ppY3n); Y3n=ZERO
    Y3p => D2STATE(:,ppY3p); Y3p=ZERO
    Ys3c => D2STATE(:,ppYs3c); Ys3c=ZERO
    Q6c => D2STATE(:,ppQ6c); Q6c=ZERO
    Q6n => D2STATE(:,ppQ6n); Q6n=ZERO
    Q6p => D2STATE(:,ppQ6p); Q6p=ZERO
    Q6s => D2STATE(:,ppQ6s); Q6s=ZERO
    Q9x => D2STATE(:,ppQ9x); Q9x=ZERO
    Qp9x => D2STATE(:,ppQp9x); Qp9x=ZERO
    QSx => D2STATE(:,ppQSx); QSx=ZERO
    Qun => D2STATE(:,ppQun); Qun=ZERO
    Q1un => D2STATE(:,ppQ1un); Q1un=ZERO
    Q2un => D2STATE(:,ppQ2un); Q2un=ZERO
    Q2c => D2STATE(:,ppQ2c); Q2c=ZERO
    Q2n => D2STATE(:,ppQ2n); Q2n=ZERO
    Q12c => D2STATE(:,ppQ12c); Q12c=ZERO
    Q12n => D2STATE(:,ppQ12n); Q12n=ZERO
    Q1c => D2STATE(:,ppQ1c); Q1c=ZERO
    Q1n => D2STATE(:,ppQ1n); Q1n=ZERO
    Q1p => D2STATE(:,ppQ1p); Q1p=ZERO
    Q11c => D2STATE(:,ppQ11c); Q11c=ZERO
    Q11n => D2STATE(:,ppQ11n); Q11n=ZERO
    Q11p => D2STATE(:,ppQ11p); Q11p=ZERO
    Q21c => D2STATE(:,ppQ21c); Q21c=ZERO
    Q21n => D2STATE(:,ppQ21n); Q21n=ZERO
    Q21p => D2STATE(:,ppQ21p); Q21p=ZERO
    H1c => D2STATE(:,ppH1c); H1c=ZERO
    H1n => D2STATE(:,ppH1n); H1n=ZERO
    H1p => D2STATE(:,ppH1p); H1p=ZERO
    H2c => D2STATE(:,ppH2c); H2c=ZERO
    H2n => D2STATE(:,ppH2n); H2n=ZERO
    H2p => D2STATE(:,ppH2p); H2p=ZERO
    H3c => D2STATE(:,ppH3c); H3c=ZERO
    H3n => D2STATE(:,ppH3n); H3n=ZERO
    H3p => D2STATE(:,ppH3p); H3p=ZERO
    HNc => D2STATE(:,ppHNc); HNc=ZERO
    HNn => D2STATE(:,ppHNn); HNn=ZERO
    HNp => D2STATE(:,ppHNp); HNp=ZERO
    Hac => D2STATE(:,ppHac); Hac=ZERO
    Kp1p => D2STATE(:,ppKp1p); Kp1p=ZERO
    Kp3n => D2STATE(:,ppKp3n); Kp3n=ZERO
    Kp4n => D2STATE(:,ppKp4n); Kp4n=ZERO
    Kn4n => D2STATE(:,ppKn4n); Kn4n=ZERO
    Kp5s => D2STATE(:,ppKp5s); Kp5s=ZERO
    Qpun => D2STATE(:,ppQpun); Qpun=ZERO
    K1p => D2STATE(:,ppK1p); K1p=ZERO
    K11p => D2STATE(:,ppK11p); K11p=ZERO
    K21p => D2STATE(:,ppK21p); K21p=ZERO
    K4n => D2STATE(:,ppK4n); K4n=ZERO
    K14n => D2STATE(:,ppK14n); K14n=ZERO
    K24n => D2STATE(:,ppK24n); K24n=ZERO
    K5s => D2STATE(:,ppK5s); K5s=ZERO
    K15s => D2STATE(:,ppK15s); K15s=ZERO
    K6r => D2STATE(:,ppK6r); K6r=ZERO
    K16r => D2STATE(:,ppK16r); K16r=ZERO
    K26r => D2STATE(:,ppK26r); K26r=ZERO
    K3n => D2STATE(:,ppK3n); K3n=ZERO
    K13n => D2STATE(:,ppK13n); K13n=ZERO
    K23n => D2STATE(:,ppK23n); K23n=ZERO
    G2o => D2STATE(:,ppG2o); G2o=ZERO
    G4n => D2STATE(:,ppG4n); G4n=ZERO
    Dfm => D2STATE(:,ppDfm); Dfm=ZERO
    Dcm => D2STATE(:,ppDcm); Dcm=ZERO
    Dlm => D2STATE(:,ppDlm); Dlm=ZERO
    D1m => D2STATE(:,ppD1m); D1m=ZERO
    D2m => D2STATE(:,ppD2m); D2m=ZERO
    D6m => D2STATE(:,ppD6m); D6m=ZERO
    D7m => D2STATE(:,ppD7m); D7m=ZERO
    D8m => D2STATE(:,ppD8m); D8m=ZERO
    D9m => D2STATE(:,ppD9m); D9m=ZERO
    DSm => D2STATE(:,ppDSm); DSm=ZERO
    DH2m => D2STATE(:,ppDH2m); DH2m=ZERO
    DH3m => D2STATE(:,ppDH3m); DH3m=ZERO
    irri_bio => D2STATE(:,ppirri_bio); irri_bio=ZERO
    G3c => D2STATE(:,ppG3c); G3c=ZERO
    G3h => D2STATE(:,ppG3h); G3h=ZERO
    G13c => D2STATE(:,ppG13c); G13c=ZERO
    G13h => D2STATE(:,ppG13h); G13h=ZERO
    G23c => D2STATE(:,ppG23c); G23c=ZERO
    G23h => D2STATE(:,ppG23h); G23h=ZERO



    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !Allocations(s) of and assigning values to alternative Z-axis

    allocate(seddepth(1:40),stat=status)
    if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "seddepth")

    call calc_sigma_depth(40,2.0D+00,0.30D+00,seddepth)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    allocate(prmdepth(1:40),stat=status)
    if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "prmdepth")

    call calc_sigma_depth(40,2.0D+00,0.05D+00,prmdepth)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of other pelagic variables which can be outputted
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-"
  write(LOGUNIT,*) "# Allocating Other Global Variables .."

    ppiNPI(iiP1)=72
    ppiNPI(iiP2)=73
    ppiNPI(iiP3)=74
    ppiNPI(iiP4)=75
    ppiNPI(iiP5)=76
    ppiNPI(iiP6)=77
    ppsugPI(iiP1)=78
    ppsugPI(iiP2)=79
    ppsugPI(iiP3)=80
    ppsugPI(iiP4)=81
    ppsugPI(iiP5)=82
    ppsugPI(iiP6)=83
    ppsunPI(iiP1)=84
    ppsunPI(iiP2)=85
    ppsunPI(iiP3)=86
    ppsunPI(iiP4)=87
    ppsunPI(iiP5)=88
    ppsunPI(iiP6)=89
    ppsdoPI(iiP1)=90
    ppsdoPI(iiP2)=91
    ppsdoPI(iiP3)=92
    ppsdoPI(iiP4)=93
    ppsdoPI(iiP5)=94
    ppsdoPI(iiP6)=95
    ppqpPc(iiP1)=96
    ppqpPc(iiP2)=97
    ppqpPc(iiP3)=98
    ppqpPc(iiP4)=99
    ppqpPc(iiP5)=100
    ppqpPc(iiP6)=101
    ppqnPc(iiP1)=102
    ppqnPc(iiP2)=103
    ppqnPc(iiP3)=104
    ppqnPc(iiP4)=105
    ppqnPc(iiP5)=106
    ppqnPc(iiP6)=107
    ppqsPc(iiP1)=108
    ppqsPc(iiP2)=109
    ppqsPc(iiP3)=110
    ppqsPc(iiP4)=111
    ppqsPc(iiP5)=112
    ppqsPc(iiP6)=113
    ppqlPc(iiP1)=114
    ppqlPc(iiP2)=115
    ppqlPc(iiP3)=116
    ppqlPc(iiP4)=117
    ppqlPc(iiP5)=118
    ppqlPc(iiP6)=119
    ppqpZc(iiZ3)=120
    ppqpZc(iiZ4)=121
    ppqpZc(iiZ2)=122
    ppqnZc(iiZ3)=123
    ppqnZc(iiZ4)=124
    ppqnZc(iiZ2)=125
    ppqp_mz(iiZ5)=126
    ppqp_mz(iiZ6)=127
    ppqn_mz(iiZ5)=128
    ppqn_mz(iiZ6)=129
    ppflPIR1n(iiP1)=130
    ppflPIR1n(iiP2)=131
    ppflPIR1n(iiP3)=132
    ppflPIR1n(iiP4)=133
    ppflPIR1n(iiP5)=134
    ppflPIR1n(iiP6)=135
    ppflPIR1p(iiP1)=136
    ppflPIR1p(iiP2)=137
    ppflPIR1p(iiP3)=138
    ppflPIR1p(iiP4)=139
    ppflPIR1p(iiP5)=140
    ppflPIR1p(iiP6)=141
    ppflPIR6n(iiP1)=142
    ppflPIR6n(iiP2)=143
    ppflPIR6n(iiP3)=144
    ppflPIR6n(iiP4)=145
    ppflPIR6n(iiP5)=146
    ppflPIR6n(iiP6)=147
    ppflPIR6p(iiP1)=148
    ppflPIR6p(iiP2)=149
    ppflPIR6p(iiP3)=150
    ppflPIR6p(iiP4)=151
    ppflPIR6p(iiP5)=152
    ppflPIR6p(iiP6)=153
    ppflPIR6s(iiP1)=154
    ppflPIR6s(iiP2)=155
    ppflPIR6s(iiP3)=156
    ppflPIR6s(iiP4)=157
    ppflPIR6s(iiP5)=158
    ppflPIR6s(iiP6)=159
    ppfr_lim_PI_n(iiP1)=160
    ppfr_lim_PI_n(iiP2)=161
    ppfr_lim_PI_n(iiP3)=162
    ppfr_lim_PI_n(iiP4)=163
    ppfr_lim_PI_n(iiP5)=164
    ppfr_lim_PI_n(iiP6)=165
    ppfr_lim_PI_p(iiP1)=166
    ppfr_lim_PI_p(iiP2)=167
    ppfr_lim_PI_p(iiP3)=168
    ppfr_lim_PI_p(iiP4)=169
    ppfr_lim_PI_p(iiP5)=170
    ppfr_lim_PI_p(iiP6)=171
    ppfl_xgrazing_PIc(iiP1)=172
    ppfl_xgrazing_PIc(iiP2)=173
    ppfl_xgrazing_PIc(iiP3)=174
    ppfl_xgrazing_PIc(iiP4)=175
    ppfl_xgrazing_PIc(iiP5)=176
    ppfl_xgrazing_PIc(iiP6)=177
    ppsediPI(iiP1)=178
    ppsediPI(iiP2)=179
    ppsediPI(iiP3)=180
    ppsediPI(iiP4)=181
    ppsediPI(iiP5)=182
    ppsediPI(iiP6)=183
    ppsediMiZ(iiZ5)=184
    ppsediMiZ(iiZ6)=185
    ppsediMeZ(iiZ3)=186
    ppsediMeZ(iiZ4)=187
    ppsediMeZ(iiZ2)=188
    ppPI_dw(iiP1)=189
    ppPI_dw(iiP2)=190
    ppPI_dw(iiP3)=191
    ppPI_dw(iiP4)=192
    ppPI_dw(iiP5)=193
    ppPI_dw(iiP6)=194
    ppeiPI(iiP1)=195
    ppeiPI(iiP2)=196
    ppeiPI(iiP3)=197
    ppeiPI(iiP4)=198
    ppeiPI(iiP5)=199
    ppeiPI(iiP6)=200
    ppEPLi(iiP1)=201
    ppEPLi(iiP2)=202
    ppEPLi(iiP3)=203
    ppEPLi(iiP4)=204
    ppEPLi(iiP5)=205
    ppEPLi(iiP6)=206

    ETW => D3DIAGNOS(:,ppETW); ETW=ZERO
    ESW => D3DIAGNOS(:,ppESW); ESW=ZERO
    EIR => D3DIAGNOS(:,ppEIR); EIR=ZERO
    ERHO => D3DIAGNOS(:,ppERHO); ERHO=ZERO
    ESS => D3DIAGNOS(:,ppESS); ESS=ZERO
    cmO2o => D3DIAGNOS(:,ppcmO2o); cmO2o=ZERO
    XO2o => D3DIAGNOS(:,ppXO2o); XO2o=ZERO
    eO2mO2 => D3DIAGNOS(:,ppeO2mO2); eO2mO2=ZERO
    Chla => D3DIAGNOS(:,ppChla); Chla=ZERO
    ws_out => D3DIAGNOS(:,ppws_out); ws_out=ZERO
    flPTN6r => D3DIAGNOS(:,ppflPTN6r); flPTN6r=ZERO
    sN4N3n => D3DIAGNOS(:,ppsN4N3n); sN4N3n=ZERO
    flN4N3n => D3DIAGNOS(:,ppflN4N3n); flN4N3n=ZERO
    flN3N4n => D3DIAGNOS(:,ppflN3N4n); flN3N4n=ZERO
    flN3O4n => D3DIAGNOS(:,ppflN3O4n); flN3O4n=ZERO
    flBaZTc => D3DIAGNOS(:,ppflBaZTc); flBaZTc=ZERO
    pNaNt => D3DIAGNOS(:,pppNaNt); pNaNt=ZERO
    flR1O3c => D3DIAGNOS(:,ppflR1O3c); flR1O3c=ZERO
    flR1N4n => D3DIAGNOS(:,ppflR1N4n); flR1N4n=ZERO
    flB1N4n => D3DIAGNOS(:,ppflB1N4n); flB1N4n=ZERO
    flPIN4n => D3DIAGNOS(:,ppflPIN4n); flPIN4n=ZERO
    rmB1c => D3DIAGNOS(:,pprmB1c); rmB1c=ZERO
    flR3R2c => D3DIAGNOS(:,ppflR3R2c); flR3R2c=ZERO
    p_xfree_R2 => D3DIAGNOS(:,ppp_xfree_R2); p_xfree_R2=ZERO
    fr_lim_BN_n => D3DIAGNOS(:,ppfr_lim_BN_n); fr_lim_BN_n=ZERO
    fr_lim_B1_n => D3DIAGNOS(:,ppfr_lim_B1_n); fr_lim_B1_n=ZERO
    fr_lim_B1_p => D3DIAGNOS(:,ppfr_lim_B1_p); fr_lim_B1_p=ZERO
    fl_xgrazing_B1c => D3DIAGNOS(:,ppfl_xgrazing_B1c); fl_xgrazing_B1c=ZERO
    rml => D3DIAGNOS(:,pprml); rml=ZERO
    qpR6c => D3DIAGNOS(:,ppqpR6c); qpR6c=ZERO
    qnR6c => D3DIAGNOS(:,ppqnR6c); qnR6c=ZERO
    qsR6c => D3DIAGNOS(:,ppqsR6c); qsR6c=ZERO
    qpB1c => D3DIAGNOS(:,ppqpB1c); qpB1c=ZERO
    qnB1c => D3DIAGNOS(:,ppqnB1c); qnB1c=ZERO
    pMIupZ4 => D3DIAGNOS(:,pppMIupZ4); pMIupZ4=ZERO
    rnetPTc => D3DIAGNOS(:,pprnetPTc); rnetPTc=ZERO
    flnDIn => D3DIAGNOS(:,ppflnDIn); flnDIn=ZERO
    flnDIp => D3DIAGNOS(:,ppflnDIp); flnDIp=ZERO
    Nun => D3DIAGNOS(:,ppNun); Nun=ZERO
    Nup => D3DIAGNOS(:,ppNup); Nup=ZERO
    sediR2 => D3DIAGNOS(:,ppsediR2); sediR2=ZERO
    sediR6 => D3DIAGNOS(:,ppsediR6); sediR6=ZERO
    sediR6s => D3DIAGNOS(:,ppsediR6s); sediR6s=ZERO
    sediRZ => D3DIAGNOS(:,ppsediRZ); sediRZ=ZERO
    sediR9 => D3DIAGNOS(:,ppsediR9); sediR9=ZERO
    sediB1 => D3DIAGNOS(:,ppsediB1); sediB1=ZERO
    POC => D3DIAGNOS(:,ppPOC); POC=ZERO
    PON => D3DIAGNOS(:,ppPON); PON=ZERO
    POP => D3DIAGNOS(:,ppPOP); POP=ZERO
    PTi => D3DIAGNOS(:,ppPTi); PTi=ZERO
    limnuti => D3DIAGNOS(:,pplimnuti); limnuti=ZERO
    O2o_vr => D3DIAGNOS(:,ppO2o_vr); O2o_vr=ZERO
    N1p_vr => D3DIAGNOS(:,ppN1p_vr); N1p_vr=ZERO
    N3n_vr => D3DIAGNOS(:,ppN3n_vr); N3n_vr=ZERO
    Chla_vr => D3DIAGNOS(:,ppChla_vr); Chla_vr=ZERO
    R6c_vr => D3DIAGNOS(:,ppR6c_vr); R6c_vr=ZERO
    xSizeMA_m => D3DIAGNOS(:,ppxSizeMA_m); xSizeMA_m=ZERO
    xSizeMA_d => D3DIAGNOS(:,ppxSizeMA_d); xSizeMA_d=ZERO
    qR2P1 => D3DIAGNOS(:,ppqR2P1); qR2P1=ZERO
    xEPS => D3DIAGNOS(:,ppxEPS); xEPS=ZERO
    xEPS_0 => D3DIAGNOS(:,ppxEPS_0); xEPS_0=ZERO
    xEPS_ESS => D3DIAGNOS(:,ppxEPS_ESS); xEPS_ESS=ZERO
    xEPS_Chl => D3DIAGNOS(:,ppxEPS_Chl); xEPS_Chl=ZERO
    pCO2 => D3DIAGNOS(:,pppCO2); pCO2=ZERO
    CO2 => D3DIAGNOS(:,ppCO2); CO2=ZERO
    HCO3 => D3DIAGNOS(:,ppHCO3); HCO3=ZERO
    CO3 => D3DIAGNOS(:,ppCO3); CO3=ZERO
    pH => D3DIAGNOS(:,pppH); pH=ZERO
    Ac => D3DIAGNOS(:,ppAc); Ac=ZERO
    CAc => D3DIAGNOS(:,ppCAc); CAc=ZERO
    DIC => D3DIAGNOS(:,ppDIC); DIC=ZERO

    iNPI => D3DIAGNOS(:,ppiNPI(iiP1): ppiNPI(iiP6))
    iNPI=ZERO
    sugPI => D3DIAGNOS(:,ppsugPI(iiP1): ppsugPI(iiP6))
    sugPI=ZERO
    sunPI => D3DIAGNOS(:,ppsunPI(iiP1): ppsunPI(iiP6))
    sunPI=ZERO
    sdoPI => D3DIAGNOS(:,ppsdoPI(iiP1): ppsdoPI(iiP6))
    sdoPI=ZERO
    qpPc => D3DIAGNOS(:,ppqpPc(iiP1): ppqpPc(iiP6))
    qpPc=ZERO
    qnPc => D3DIAGNOS(:,ppqnPc(iiP1): ppqnPc(iiP6))
    qnPc=ZERO
    qsPc => D3DIAGNOS(:,ppqsPc(iiP1): ppqsPc(iiP6))
    qsPc=ZERO
    qlPc => D3DIAGNOS(:,ppqlPc(iiP1): ppqlPc(iiP6))
    qlPc=ZERO
    qpZc => D3DIAGNOS(:,ppqpZc(iiZ3): ppqpZc(iiZ2))
    qpZc=ZERO
    qnZc => D3DIAGNOS(:,ppqnZc(iiZ3): ppqnZc(iiZ2))
    qnZc=ZERO
    qp_mz => D3DIAGNOS(:,ppqp_mz(iiZ5): ppqp_mz(iiZ6))
    qp_mz=ZERO
    qn_mz => D3DIAGNOS(:,ppqn_mz(iiZ5): ppqn_mz(iiZ6))
    qn_mz=ZERO
    flPIR1n => D3DIAGNOS(:,ppflPIR1n(iiP1): ppflPIR1n(iiP6))
    flPIR1n=ZERO
    flPIR1p => D3DIAGNOS(:,ppflPIR1p(iiP1): ppflPIR1p(iiP6))
    flPIR1p=ZERO
    flPIR6n => D3DIAGNOS(:,ppflPIR6n(iiP1): ppflPIR6n(iiP6))
    flPIR6n=ZERO
    flPIR6p => D3DIAGNOS(:,ppflPIR6p(iiP1): ppflPIR6p(iiP6))
    flPIR6p=ZERO
    flPIR6s => D3DIAGNOS(:,ppflPIR6s(iiP1): ppflPIR6s(iiP6))
    flPIR6s=ZERO
    fr_lim_PI_n => D3DIAGNOS(:,ppfr_lim_PI_n(iiP1): ppfr_lim_PI_n(iiP6))
    fr_lim_PI_n=ZERO
    fr_lim_PI_p => D3DIAGNOS(:,ppfr_lim_PI_p(iiP1): ppfr_lim_PI_p(iiP6))
    fr_lim_PI_p=ZERO
    fl_xgrazing_PIc => D3DIAGNOS(:,ppfl_xgrazing_PIc(iiP1):&
     & ppfl_xgrazing_PIc(iiP6))
    fl_xgrazing_PIc=ZERO
    sediPI => D3DIAGNOS(:,ppsediPI(iiP1): ppsediPI(iiP6))
    sediPI=ZERO
    sediMiZ => D3DIAGNOS(:,ppsediMiZ(iiZ5): ppsediMiZ(iiZ6))
    sediMiZ=ZERO
    sediMeZ => D3DIAGNOS(:,ppsediMeZ(iiZ3): ppsediMeZ(iiZ2))
    sediMeZ=ZERO
    PI_dw => D3DIAGNOS(:,ppPI_dw(iiP1): ppPI_dw(iiP6))
    PI_dw=ZERO
    eiPI => D3DIAGNOS(:,ppeiPI(iiP1): ppeiPI(iiP6))
    eiPI=ZERO
    EPLi => D3DIAGNOS(:,ppEPLi(iiP1): ppEPLi(iiP6))
    EPLi=ZERO

    PrsM1p => D3DIAGNOS_PRF(:,ppPrsM1p); PrsM1p=ZERO
    PrM1p => D3DIAGNOS_PRF(:,ppPrM1p); PrM1p=ZERO
    PrM3n => D3DIAGNOS_PRF(:,ppPrM3n); PrM3n=ZERO
    PrM4n => D3DIAGNOS_PRF(:,ppPrM4n); PrM4n=ZERO
    PreM5s => D3DIAGNOS_PRF(:,ppPreM5s); PreM5s=ZERO
    PrM5s => D3DIAGNOS_PRF(:,ppPrM5s); PrM5s=ZERO
    PrM6r => D3DIAGNOS_PRF(:,ppPrM6r); PrM6r=ZERO
    PrQun => D3DIAGNOS_PRF(:,ppPrQun); PrQun=ZERO
    PrDIC => D3DIAGNOS_PRF(:,ppPrDIC); PrDIC=ZERO
    PrAc => D3DIAGNOS_PRF(:,ppPrAc); PrAc=ZERO
    PrpH => D3DIAGNOS_PRF(:,ppPrpH); PrpH=ZERO
    PrBDc => D3DIAGNOS_PRF(:,ppPrBDc); PrBDc=ZERO
    PrBChlC => D3DIAGNOS_PRF(:,ppPrBChlC); PrBChlC=ZERO


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of other benthic variables which can be outputted 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    ppjbotBPc(iiBP1)=212
    ppsunBI(iiBP1)=213
    ppsugBI(iiBP1)=214
    ppjBTQIc(iiQ1)=215
    ppjBTQIc(iiQ11)=216
    ppjBTQIc(iiQ21)=217
    ppjBTQIn(iiQ1)=218
    ppjBTQIn(iiQ11)=219
    ppjBTQIn(iiQ21)=220
    ppjBTQIp(iiQ1)=221
    ppjBTQIp(iiQ11)=222
    ppjBTQIp(iiQ21)=223
    ppjQIBTc(iiQ1)=224
    ppjQIBTc(iiQ11)=225
    ppjQIBTc(iiQ21)=226
    ppjQIBTn(iiQ1)=227
    ppjQIBTn(iiQ11)=228
    ppjQIBTn(iiQ21)=229
    ppjQIBTp(iiQ1)=230
    ppjQIBTp(iiQ11)=231
    ppjQIBTp(iiQ21)=232
    ppjnetHIc(iiH1)=233
    ppjnetHIc(iiH2)=234
    ppjnetHIc(iiH3)=235
    ppjnetHIc(iiHN)=236
    ppjugYIc(iiY1)=237
    ppjugYIc(iiY2)=238
    ppjugYIc(iiY4)=239
    ppjugYIc(iiY5)=240
    ppjnetYIc(iiY1)=241
    ppjnetYIc(iiY2)=242
    ppjnetYIc(iiY4)=243
    ppjnetYIc(iiY5)=244
    ppjmYIc(iiY1)=245
    ppjmYIc(iiY2)=246
    ppjmYIc(iiY4)=247
    ppjmYIc(iiY5)=248
    ppjmY3c(iiYy3)=249
    ppjmY3c(iiY3)=250
    ppjrrY3c(iiYy3)=251
    ppjrrY3c(iiY3)=252
    ppjrrYIc(iiY1)=253
    ppjrrYIc(iiY2)=254
    ppjrrYIc(iiY4)=255
    ppjrrYIc(iiY5)=256
    pprugY3c(iiYy3)=257
    pprugY3c(iiY3)=258
    ppefsatY3(iiYy3)=259
    ppefsatY3(iiY3)=260
    ppfr_lim_HI_n(iiH1)=261
    ppfr_lim_HI_n(iiH2)=262
    ppfr_lim_HI_n(iiH3)=263
    ppfr_lim_HI_n(iiHN)=264
    ppfr_lim_HI_p(iiH1)=265
    ppfr_lim_HI_p(iiH2)=266
    ppfr_lim_HI_p(iiH3)=267
    ppfr_lim_HI_p(iiHN)=268
    ppfr_lim_HI_o(iiH1)=269
    ppfr_lim_HI_o(iiH2)=270
    ppfr_lim_HI_o(iiH3)=271
    ppfr_lim_HI_o(iiHN)=272
    ppfr_lim_BPI_n(iiBP1)=273
    ppfr_lim_BPI_p(iiBP1)=274
    ppZI_Fc(iiZ5)=275
    ppZI_Fc(iiZ6)=276
    ppZI_Fn(iiZ5)=277
    ppZI_Fn(iiZ6)=278
    ppZI_Fp(iiZ5)=279
    ppZI_Fp(iiZ6)=280
    ppjPIY3c(iiP1)=281
    ppjPIY3c(iiP2)=282
    ppjPIY3c(iiP3)=283
    ppjPIY3c(iiP4)=284
    ppjPIY3c(iiP5)=285
    ppjPIY3c(iiP6)=286
    ppjP6Y3c(iiYy3)=287
    ppjP6Y3c(iiY3)=288
    ppjZEY3c(iiZ3)=289
    ppjZEY3c(iiZ4)=290
    ppjZEY3c(iiZ2)=291
    ppjPIQ6s(iiBP1)=292
    ppPI_Benc(iiP1)=293
    ppPI_Benc(iiP2)=294
    ppPI_Benc(iiP3)=295
    ppPI_Benc(iiP4)=296
    ppPI_Benc(iiP5)=297
    ppPI_Benc(iiP6)=298
    ppPI_Benn(iiP1)=299
    ppPI_Benn(iiP2)=300
    ppPI_Benn(iiP3)=301
    ppPI_Benn(iiP4)=302
    ppPI_Benn(iiP5)=303
    ppPI_Benn(iiP6)=304
    ppPI_Benp(iiP1)=305
    ppPI_Benp(iiP2)=306
    ppPI_Benp(iiP3)=307
    ppPI_Benp(iiP4)=308
    ppPI_Benp(iiP5)=309
    ppPI_Benp(iiP6)=310
    ppPI_Benl(iiP1)=311
    ppPI_Benl(iiP2)=312
    ppPI_Benl(iiP3)=313
    ppPI_Benl(iiP4)=314
    ppPI_Benl(iiP5)=315
    ppPI_Benl(iiP6)=316
    ppPI_Bens(iiP1)=317
    ppPI_Bens(iiP2)=318
    ppPI_Bens(iiP3)=319
    ppPI_Bens(iiP4)=320
    ppPI_Bens(iiP5)=321
    ppPI_Bens(iiP6)=322
    ppZE_Benc(iiZ3)=323
    ppZE_Benc(iiZ4)=324
    ppZE_Benc(iiZ2)=325
    ppZE_Benn(iiZ3)=326
    ppZE_Benn(iiZ4)=327
    ppZE_Benn(iiZ2)=328
    ppZE_Benp(iiZ3)=329
    ppZE_Benp(iiZ4)=330
    ppZE_Benp(iiZ2)=331
    pppuPIY3(iiP1)=332
    pppuPIY3(iiP2)=333
    pppuPIY3(iiP3)=334
    pppuPIY3(iiP4)=335
    pppuPIY3(iiP5)=336
    pppuPIY3(iiP6)=337
    pppuZEY3(iiZ3)=338
    pppuZEY3(iiZ4)=339
    pppuZEY3(iiZ2)=340
    pppuP6Y3(iiYy3)=341
    pppuP6Y3(iiY3)=342
    ppsediPI_Ben(iiP1)=343
    ppsediPI_Ben(iiP2)=344
    ppsediPI_Ben(iiP3)=345
    ppsediPI_Ben(iiP4)=346
    ppsediPI_Ben(iiP5)=347
    ppsediPI_Ben(iiP6)=348
    ppsediZE_Ben(iiZ3)=349
    ppsediZE_Ben(iiZ4)=350
    ppsediZE_Ben(iiZ2)=351
    ppjnetPIc(iiP1)=352
    ppjnetPIc(iiP2)=353
    ppjnetPIc(iiP3)=354
    ppjnetPIc(iiP4)=355
    ppjnetPIc(iiP5)=356
    ppjnetPIc(iiP6)=357

    EIRr => D2DIAGNOS(:,ppEIRr); EIRr=ZERO
    EUWIND => D2DIAGNOS(:,ppEUWIND); EUWIND=ZERO
    EVWIND => D2DIAGNOS(:,ppEVWIND); EVWIND=ZERO
    dry_z => D2DIAGNOS(:,ppdry_z); dry_z=ZERO
    ETAUB => D2DIAGNOS(:,ppETAUB); ETAUB=ZERO
    EUCURR_LEVEL => D2DIAGNOS(:,ppEUCURR_LEVEL); EUCURR_LEVEL=ZERO
    EVCURR_LEVEL => D2DIAGNOS(:,ppEVCURR_LEVEL); EVCURR_LEVEL=ZERO
    jrESS => D2DIAGNOS(:,ppjrESS); jrESS=ZERO
    shiftDlm => D2DIAGNOS(:,ppshiftDlm); shiftDlm=ZERO
    shiftD1m => D2DIAGNOS(:,ppshiftD1m); shiftD1m=ZERO
    shiftD2m => D2DIAGNOS(:,ppshiftD2m); shiftD2m=ZERO
    G2_xavail_o => D2DIAGNOS(:,ppG2_xavail_o); G2_xavail_o=ZERO
    jO2Y2o => D2DIAGNOS(:,ppjO2Y2o); jO2Y2o=ZERO
    jcrrBTo => D2DIAGNOS(:,ppjcrrBTo); jcrrBTo=ZERO
    rrBTo => D2DIAGNOS(:,pprrBTo); rrBTo=ZERO
    reK6o => D2DIAGNOS(:,ppreK6o); reK6o=ZERO
    rrDTo => D2DIAGNOS(:,pprrDTo); rrDTo=ZERO
    rrATo => D2DIAGNOS(:,pprrATo); rrATo=ZERO
    reBTo => D2DIAGNOS(:,ppreBTo); reBTo=ZERO
    ruBPc => D2DIAGNOS(:,ppruBPc); ruBPc=ZERO
    ruBPn => D2DIAGNOS(:,ppruBPn); ruBPn=ZERO
    ruBPp => D2DIAGNOS(:,ppruBPp); ruBPp=ZERO
    ruBPs => D2DIAGNOS(:,ppruBPs); ruBPs=ZERO
    ruBPn3 => D2DIAGNOS(:,ppruBPn3); ruBPn3=ZERO
    ruBTc => D2DIAGNOS(:,ppruBTc); ruBTc=ZERO
    ruBTn => D2DIAGNOS(:,ppruBTn); ruBTn=ZERO
    ruBTp => D2DIAGNOS(:,ppruBTp); ruBTp=ZERO
    ruBTs => D2DIAGNOS(:,ppruBTs); ruBTs=ZERO
    reBTc => D2DIAGNOS(:,ppreBTc); reBTc=ZERO
    reBTn => D2DIAGNOS(:,ppreBTn); reBTn=ZERO
    reBTp => D2DIAGNOS(:,ppreBTp); reBTp=ZERO
    reBTs => D2DIAGNOS(:,ppreBTs); reBTs=ZERO
    reDTn => D2DIAGNOS(:,ppreDTn); reDTn=ZERO
    reDTp => D2DIAGNOS(:,ppreDTp); reDTp=ZERO
    reATc => D2DIAGNOS(:,ppreATc); reATc=ZERO
    reATn => D2DIAGNOS(:,ppreATn); reATn=ZERO
    reATp => D2DIAGNOS(:,ppreATp); reATp=ZERO
    reATs => D2DIAGNOS(:,ppreATs); reATs=ZERO
    xEPS_Sedi => D2DIAGNOS(:,ppxEPS_Sedi); xEPS_Sedi=ZERO
    jQuBTn => D2DIAGNOS(:,ppjQuBTn); jQuBTn=ZERO
    jBTQun => D2DIAGNOS(:,ppjBTQun); jBTQun=ZERO
    jrrPTc => D2DIAGNOS(:,ppjrrPTc); jrrPTc=ZERO
    jrrMec => D2DIAGNOS(:,ppjrrMec); jrrMec=ZERO
    jrrMic => D2DIAGNOS(:,ppjrrMic); jrrMic=ZERO
    irrenh => D2DIAGNOS(:,ppirrenh); irrenh=ZERO
    turenh => D2DIAGNOS(:,ppturenh); turenh=ZERO
    pxturinD1 => D2DIAGNOS(:,pppxturinD1); pxturinD1=ZERO
    jG2K3o => D2DIAGNOS(:,ppjG2K3o); jG2K3o=ZERO
    jG2K7o => D2DIAGNOS(:,ppjG2K7o); jG2K7o=ZERO
    M1p => D2DIAGNOS(:,ppM1p); M1p=ZERO
    M11p => D2DIAGNOS(:,ppM11p); M11p=ZERO
    M21p => D2DIAGNOS(:,ppM21p); M21p=ZERO
    M4n => D2DIAGNOS(:,ppM4n); M4n=ZERO
    M14n => D2DIAGNOS(:,ppM14n); M14n=ZERO
    M24n => D2DIAGNOS(:,ppM24n); M24n=ZERO
    M3n => D2DIAGNOS(:,ppM3n); M3n=ZERO
    M5s => D2DIAGNOS(:,ppM5s); M5s=ZERO
    M15s => D2DIAGNOS(:,ppM15s); M15s=ZERO
    M6r => D2DIAGNOS(:,ppM6r); M6r=ZERO
    Mp1p => D2DIAGNOS(:,ppMp1p); Mp1p=ZERO
    Mp3n => D2DIAGNOS(:,ppMp3n); Mp3n=ZERO
    Mp4n => D2DIAGNOS(:,ppMp4n); Mp4n=ZERO
    Mp5s => D2DIAGNOS(:,ppMp5s); Mp5s=ZERO
    fr_lim_Ha_n => D2DIAGNOS(:,ppfr_lim_Ha_n); fr_lim_Ha_n=ZERO
    fr_lim_Ha_o => D2DIAGNOS(:,ppfr_lim_Ha_o); fr_lim_Ha_o=ZERO
    cNIBTc => D2DIAGNOS(:,ppcNIBTc); cNIBTc=ZERO
    RI_Fc => D2DIAGNOS(:,ppRI_Fc); RI_Fc=ZERO
    RI_Fn => D2DIAGNOS(:,ppRI_Fn); RI_Fn=ZERO
    RI_Fp => D2DIAGNOS(:,ppRI_Fp); RI_Fp=ZERO
    RI_Fs => D2DIAGNOS(:,ppRI_Fs); RI_Fs=ZERO
    jPTY3c => D2DIAGNOS(:,ppjPTY3c); jPTY3c=ZERO
    jPTY3n => D2DIAGNOS(:,ppjPTY3n); jPTY3n=ZERO
    jPTY3p => D2DIAGNOS(:,ppjPTY3p); jPTY3p=ZERO
    jPTZTn => D2DIAGNOS(:,ppjPTZTn); jPTZTn=ZERO
    jPTZTp => D2DIAGNOS(:,ppjPTZTp); jPTZTp=ZERO
    jZIY3c => D2DIAGNOS(:,ppjZIY3c); jZIY3c=ZERO
    jRIQIc => D2DIAGNOS(:,ppjRIQIc); jRIQIc=ZERO
    jRIQIs => D2DIAGNOS(:,ppjRIQIs); jRIQIs=ZERO
    jY3RIc => D2DIAGNOS(:,ppjY3RIc); jY3RIc=ZERO
    jY3RIs => D2DIAGNOS(:,ppjY3RIs); jY3RIs=ZERO
    jY3O3c => D2DIAGNOS(:,ppjY3O3c); jY3O3c=ZERO
    jY3N4n => D2DIAGNOS(:,ppjY3N4n); jY3N4n=ZERO
    jY3N1p => D2DIAGNOS(:,ppjY3N1p); jY3N1p=ZERO
    jbotQ6c => D2DIAGNOS(:,ppjbotQ6c); jbotQ6c=ZERO
    jbotQ6n => D2DIAGNOS(:,ppjbotQ6n); jbotQ6n=ZERO
    jbotQ6p => D2DIAGNOS(:,ppjbotQ6p); jbotQ6p=ZERO
    jbotQ6s => D2DIAGNOS(:,ppjbotQ6s); jbotQ6s=ZERO
    jbotQ2c => D2DIAGNOS(:,ppjbotQ2c); jbotQ2c=ZERO
    jbotQ2n => D2DIAGNOS(:,ppjbotQ2n); jbotQ2n=ZERO
    Depth_Ben => D2DIAGNOS(:,ppDepth_Ben); Depth_Ben=ZERO
    EIR_Ben => D2DIAGNOS(:,ppEIR_Ben); EIR_Ben=ZERO
    ETW_Ben => D2DIAGNOS(:,ppETW_Ben); ETW_Ben=ZERO
    ERHO_Ben => D2DIAGNOS(:,ppERHO_Ben); ERHO_Ben=ZERO
    ESW_Ben => D2DIAGNOS(:,ppESW_Ben); ESW_Ben=ZERO
    Ru_Benn => D2DIAGNOS(:,ppRu_Benn); Ru_Benn=ZERO
    R3_Benc => D2DIAGNOS(:,ppR3_Benc); R3_Benc=ZERO
    R3_Benn => D2DIAGNOS(:,ppR3_Benn); R3_Benn=ZERO
    R3_Benp => D2DIAGNOS(:,ppR3_Benp); R3_Benp=ZERO
    R2_Benc => D2DIAGNOS(:,ppR2_Benc); R2_Benc=ZERO
    R2_Benn => D2DIAGNOS(:,ppR2_Benn); R2_Benn=ZERO
    O2o_Ben => D2DIAGNOS(:,ppO2o_Ben); O2o_Ben=ZERO
    cmO2o_Ben => D2DIAGNOS(:,ppcmO2o_Ben); cmO2o_Ben=ZERO
    N1p_Ben => D2DIAGNOS(:,ppN1p_Ben); N1p_Ben=ZERO
    N3n_Ben => D2DIAGNOS(:,ppN3n_Ben); N3n_Ben=ZERO
    N4n_Ben => D2DIAGNOS(:,ppN4n_Ben); N4n_Ben=ZERO
    N5s_Ben => D2DIAGNOS(:,ppN5s_Ben); N5s_Ben=ZERO
    N6r_Ben => D2DIAGNOS(:,ppN6r_Ben); N6r_Ben=ZERO
    sediR6_Ben => D2DIAGNOS(:,ppsediR6_Ben); sediR6_Ben=ZERO
    sediR2_Ben => D2DIAGNOS(:,ppsediR2_Ben); sediR2_Ben=ZERO
    efilP6Y3 => D2DIAGNOS(:,ppefilP6Y3); efilP6Y3=ZERO
    efilPART => D2DIAGNOS(:,ppefilPART); efilPART=ZERO
    ctfPm2c => D2DIAGNOS(:,ppctfPm2c); ctfPm2c=ZERO
    ctfZem2c => D2DIAGNOS(:,ppctfZem2c); ctfZem2c=ZERO
    ctfZim2c => D2DIAGNOS(:,ppctfZim2c); ctfZim2c=ZERO
    sK4K3 => D2DIAGNOS(:,ppsK4K3); sK4K3=ZERO
    jK4K3n => D2DIAGNOS(:,ppjK4K3n); jK4K3n=ZERO
    jKuK4n => D2DIAGNOS(:,ppjKuK4n); jKuK4n=ZERO
    jK3G4n => D2DIAGNOS(:,ppjK3G4n); jK3G4n=ZERO
    jK23G4n => D2DIAGNOS(:,ppjK23G4n); jK23G4n=ZERO
    jK31K21p => D2DIAGNOS(:,ppjK31K21p); jK31K21p=ZERO
    jK34K24n => D2DIAGNOS(:,ppjK34K24n); jK34K24n=ZERO
    jK13K3n => D2DIAGNOS(:,ppjK13K3n); jK13K3n=ZERO
    jK23K13n => D2DIAGNOS(:,ppjK23K13n); jK23K13n=ZERO
    jK25K15s => D2DIAGNOS(:,ppjK25K15s); jK25K15s=ZERO
    jK36K26r => D2DIAGNOS(:,ppjK36K26r); jK36K26r=ZERO
    totPELc => D2DIAGNOS(:,pptotPELc); totPELc=ZERO
    totPELn => D2DIAGNOS(:,pptotPELn); totPELn=ZERO
    totPELp => D2DIAGNOS(:,pptotPELp); totPELp=ZERO
    totPELs => D2DIAGNOS(:,pptotPELs); totPELs=ZERO
    totBENc => D2DIAGNOS(:,pptotBENc); totBENc=ZERO
    totBENn => D2DIAGNOS(:,pptotBENn); totBENn=ZERO
    totBENp => D2DIAGNOS(:,pptotBENp); totBENp=ZERO
    totBENs => D2DIAGNOS(:,pptotBENs); totBENs=ZERO
    totSYSc => D2DIAGNOS(:,pptotSYSc); totSYSc=ZERO
    totSYSn => D2DIAGNOS(:,pptotSYSn); totSYSn=ZERO
    totSYSp => D2DIAGNOS(:,pptotSYSp); totSYSp=ZERO
    totSYSs => D2DIAGNOS(:,pptotSYSs); totSYSs=ZERO
    jtPelc => D2DIAGNOS(:,ppjtPelc); jtPelc=ZERO
    jtPeln => D2DIAGNOS(:,ppjtPeln); jtPeln=ZERO
    jtPelp => D2DIAGNOS(:,ppjtPelp); jtPelp=ZERO
    jtPels => D2DIAGNOS(:,ppjtPels); jtPels=ZERO
    jtBenc => D2DIAGNOS(:,ppjtBenc); jtBenc=ZERO
    jtBenn => D2DIAGNOS(:,ppjtBenn); jtBenn=ZERO
    jtBenp => D2DIAGNOS(:,ppjtBenp); jtBenp=ZERO
    jtBens => D2DIAGNOS(:,ppjtBens); jtBens=ZERO
    jtotbenpelc => D2DIAGNOS(:,ppjtotbenpelc); jtotbenpelc=ZERO
    jtotbenpeln => D2DIAGNOS(:,ppjtotbenpeln); jtotbenpeln=ZERO
    jtotbenpelp => D2DIAGNOS(:,ppjtotbenpelp); jtotbenpelp=ZERO
    jtotbenpels => D2DIAGNOS(:,ppjtotbenpels); jtotbenpels=ZERO
    jupPELc => D2DIAGNOS(:,ppjupPELc); jupPELc=ZERO
    jminPELc => D2DIAGNOS(:,ppjminPELc); jminPELc=ZERO
    jupBENc => D2DIAGNOS(:,ppjupBENc); jupBENc=ZERO
    jminBENc => D2DIAGNOS(:,ppjminBENc); jminBENc=ZERO
    totPELInc => D2DIAGNOS(:,pptotPELInc); totPELInc=ZERO
    totBENInc => D2DIAGNOS(:,pptotBENInc); totBENInc=ZERO
    totPELh => D2DIAGNOS(:,pptotPELh); totPELh=ZERO
    totBENh => D2DIAGNOS(:,pptotBENh); totBENh=ZERO
    jsdoMesoc => D2DIAGNOS(:,ppjsdoMesoc); jsdoMesoc=ZERO
    jrsMicroc => D2DIAGNOS(:,ppjrsMicroc); jrsMicroc=ZERO
    jnetMeZc => D2DIAGNOS(:,ppjnetMeZc); jnetMeZc=ZERO
    jnetMiZc => D2DIAGNOS(:,ppjnetMiZc); jnetMiZc=ZERO
    jnetB1c => D2DIAGNOS(:,ppjnetB1c); jnetB1c=ZERO
    jnetY3c => D2DIAGNOS(:,ppjnetY3c); jnetY3c=ZERO
    jspaY3c => D2DIAGNOS(:,ppjspaY3c); jspaY3c=ZERO
    jnetYy3c => D2DIAGNOS(:,ppjnetYy3c); jnetYy3c=ZERO
    jnetBPTc => D2DIAGNOS(:,ppjnetBPTc); jnetBPTc=ZERO
    jPLO3c => D2DIAGNOS(:,ppjPLO3c); jPLO3c=ZERO
    jZIR6n => D2DIAGNOS(:,ppjZIR6n); jZIR6n=ZERO
    jZIR6p => D2DIAGNOS(:,ppjZIR6p); jZIR6p=ZERO
    jZIDIn => D2DIAGNOS(:,ppjZIDIn); jZIDIn=ZERO
    jZIDIp => D2DIAGNOS(:,ppjZIDIp); jZIDIp=ZERO
    jCaCO3Y3c => D2DIAGNOS(:,ppjCaCO3Y3c); jCaCO3Y3c=ZERO
    Output2d_1 => D2DIAGNOS(:,ppOutput2d_1); Output2d_1=ZERO
    Output2d_2 => D2DIAGNOS(:,ppOutput2d_2); Output2d_2=ZERO
    Output2d_3 => D2DIAGNOS(:,ppOutput2d_3); Output2d_3=ZERO
    Output2d_4 => D2DIAGNOS(:,ppOutput2d_4); Output2d_4=ZERO
    jPelFishInput => D2DIAGNOS(:,ppjPelFishInput); jPelFishInput=ZERO
    jBenFishInput => D2DIAGNOS(:,ppjBenFishInput); jBenFishInput=ZERO
    SdTdzTh => D2DIAGNOS(:,ppSdTdzTh); SdTdzTh=ZERO
    TMLd => D2DIAGNOS(:,ppTMLd); TMLd=ZERO
    BMLd => D2DIAGNOS(:,ppBMLd); BMLd=ZERO
    Hs_out => D2DIAGNOS(:,ppHs_out); Hs_out=ZERO
    Tz_out => D2DIAGNOS(:,ppTz_out); Tz_out=ZERO
    u_orb_out => D2DIAGNOS(:,ppu_orb_out); u_orb_out=ZERO
    TauW => D2DIAGNOS(:,ppTauW); TauW=ZERO
    TauC => D2DIAGNOS(:,ppTauC); TauC=ZERO
    TauBed => D2DIAGNOS(:,ppTauBed); TauBed=ZERO
    eta_out => D2DIAGNOS(:,ppeta_out); eta_out=ZERO
    DICae => D2DIAGNOS(:,ppDICae); DICae=ZERO
    DICan => D2DIAGNOS(:,ppDICan); DICan=ZERO
    O3c_Ben => D2DIAGNOS(:,ppO3c_Ben); O3c_Ben=ZERO
    O3h_Ben => D2DIAGNOS(:,ppO3h_Ben); O3h_Ben=ZERO
    ctO3m2h => D2DIAGNOS(:,ppctO3m2h); ctO3m2h=ZERO
    CAcae => D2DIAGNOS(:,ppCAcae); CAcae=ZERO
    CAcan => D2DIAGNOS(:,ppCAcan); CAcan=ZERO
    Acae => D2DIAGNOS(:,ppAcae); Acae=ZERO
    Acan => D2DIAGNOS(:,ppAcan); Acan=ZERO
    pHae => D2DIAGNOS(:,pppHae); pHae=ZERO
    pHdn => D2DIAGNOS(:,pppHdn); pHdn=ZERO
    pHan => D2DIAGNOS(:,pppHan); pHan=ZERO
    pCO2ae => D2DIAGNOS(:,pppCO2ae); pCO2ae=ZERO
    pCO2an => D2DIAGNOS(:,pppCO2an); pCO2an=ZERO
    CO2ae => D2DIAGNOS(:,ppCO2ae); CO2ae=ZERO
    HCO3ae => D2DIAGNOS(:,ppHCO3ae); HCO3ae=ZERO
    CO3ae => D2DIAGNOS(:,ppCO3ae); CO3ae=ZERO
    DIC_Ben => D2DIAGNOS(:,ppDIC_Ben); DIC_Ben=ZERO
    CO2an => D2DIAGNOS(:,ppCO2an); CO2an=ZERO
    HCO3an => D2DIAGNOS(:,ppHCO3an); HCO3an=ZERO
    CO3an => D2DIAGNOS(:,ppCO3an); CO3an=ZERO
    jG33G23h => D2DIAGNOS(:,ppjG33G23h); jG33G23h=ZERO
    jG33G23c => D2DIAGNOS(:,ppjG33G23c); jG33G23c=ZERO

    jbotBPc => D2DIAGNOS(:,ppjbotBPc(iiBP1): ppjbotBPc(iiBP1))
    jbotBPc=ZERO
    sunBI => D2DIAGNOS(:,ppsunBI(iiBP1): ppsunBI(iiBP1))
    sunBI=ZERO
    sugBI => D2DIAGNOS(:,ppsugBI(iiBP1): ppsugBI(iiBP1))
    sugBI=ZERO
    jBTQIc => D2DIAGNOS(:,ppjBTQIc(iiQ1): ppjBTQIc(iiQ21))
    jBTQIc=ZERO
    jBTQIn => D2DIAGNOS(:,ppjBTQIn(iiQ1): ppjBTQIn(iiQ21))
    jBTQIn=ZERO
    jBTQIp => D2DIAGNOS(:,ppjBTQIp(iiQ1): ppjBTQIp(iiQ21))
    jBTQIp=ZERO
    jQIBTc => D2DIAGNOS(:,ppjQIBTc(iiQ1): ppjQIBTc(iiQ21))
    jQIBTc=ZERO
    jQIBTn => D2DIAGNOS(:,ppjQIBTn(iiQ1): ppjQIBTn(iiQ21))
    jQIBTn=ZERO
    jQIBTp => D2DIAGNOS(:,ppjQIBTp(iiQ1): ppjQIBTp(iiQ21))
    jQIBTp=ZERO
    jnetHIc => D2DIAGNOS(:,ppjnetHIc(iiH1): ppjnetHIc(iiHN))
    jnetHIc=ZERO
    jugYIc => D2DIAGNOS(:,ppjugYIc(iiY1): ppjugYIc(iiY5))
    jugYIc=ZERO
    jnetYIc => D2DIAGNOS(:,ppjnetYIc(iiY1): ppjnetYIc(iiY5))
    jnetYIc=ZERO
    jmYIc => D2DIAGNOS(:,ppjmYIc(iiY1): ppjmYIc(iiY5))
    jmYIc=ZERO
    jmY3c => D2DIAGNOS(:,ppjmY3c(iiYy3): ppjmY3c(iiY3))
    jmY3c=ZERO
    jrrY3c => D2DIAGNOS(:,ppjrrY3c(iiYy3): ppjrrY3c(iiY3))
    jrrY3c=ZERO
    jrrYIc => D2DIAGNOS(:,ppjrrYIc(iiY1): ppjrrYIc(iiY5))
    jrrYIc=ZERO
    rugY3c => D2DIAGNOS(:,pprugY3c(iiYy3): pprugY3c(iiY3))
    rugY3c=ZERO
    efsatY3 => D2DIAGNOS(:,ppefsatY3(iiYy3): ppefsatY3(iiY3))
    efsatY3=ZERO
    fr_lim_HI_n => D2DIAGNOS(:,ppfr_lim_HI_n(iiH1): ppfr_lim_HI_n(iiHN))
    fr_lim_HI_n=ZERO
    fr_lim_HI_p => D2DIAGNOS(:,ppfr_lim_HI_p(iiH1): ppfr_lim_HI_p(iiHN))
    fr_lim_HI_p=ZERO
    fr_lim_HI_o => D2DIAGNOS(:,ppfr_lim_HI_o(iiH1): ppfr_lim_HI_o(iiHN))
    fr_lim_HI_o=ZERO
    fr_lim_BPI_n => D2DIAGNOS(:,ppfr_lim_BPI_n(iiBP1):&
     & ppfr_lim_BPI_n(iiBP1))
    fr_lim_BPI_n=ZERO
    fr_lim_BPI_p => D2DIAGNOS(:,ppfr_lim_BPI_p(iiBP1):&
     & ppfr_lim_BPI_p(iiBP1))
    fr_lim_BPI_p=ZERO
    ZI_Fc => D2DIAGNOS(:,ppZI_Fc(iiZ5): ppZI_Fc(iiZ6))
    ZI_Fc=ZERO
    ZI_Fn => D2DIAGNOS(:,ppZI_Fn(iiZ5): ppZI_Fn(iiZ6))
    ZI_Fn=ZERO
    ZI_Fp => D2DIAGNOS(:,ppZI_Fp(iiZ5): ppZI_Fp(iiZ6))
    ZI_Fp=ZERO
    jPIY3c => D2DIAGNOS(:,ppjPIY3c(iiP1): ppjPIY3c(iiP6))
    jPIY3c=ZERO
    jP6Y3c => D2DIAGNOS(:,ppjP6Y3c(iiYy3): ppjP6Y3c(iiY3))
    jP6Y3c=ZERO
    jZEY3c => D2DIAGNOS(:,ppjZEY3c(iiZ3): ppjZEY3c(iiZ2))
    jZEY3c=ZERO
    jPIQ6s => D2DIAGNOS(:,ppjPIQ6s(iiBP1): ppjPIQ6s(iiBP1))
    jPIQ6s=ZERO
    PI_Benc => D2DIAGNOS(:,ppPI_Benc(iiP1): ppPI_Benc(iiP6))
    PI_Benc=ZERO
    PI_Benn => D2DIAGNOS(:,ppPI_Benn(iiP1): ppPI_Benn(iiP6))
    PI_Benn=ZERO
    PI_Benp => D2DIAGNOS(:,ppPI_Benp(iiP1): ppPI_Benp(iiP6))
    PI_Benp=ZERO
    PI_Benl => D2DIAGNOS(:,ppPI_Benl(iiP1): ppPI_Benl(iiP6))
    PI_Benl=ZERO
    PI_Bens => D2DIAGNOS(:,ppPI_Bens(iiP1): ppPI_Bens(iiP6))
    PI_Bens=ZERO
    ZE_Benc => D2DIAGNOS(:,ppZE_Benc(iiZ3): ppZE_Benc(iiZ2))
    ZE_Benc=ZERO
    ZE_Benn => D2DIAGNOS(:,ppZE_Benn(iiZ3): ppZE_Benn(iiZ2))
    ZE_Benn=ZERO
    ZE_Benp => D2DIAGNOS(:,ppZE_Benp(iiZ3): ppZE_Benp(iiZ2))
    ZE_Benp=ZERO
    puPIY3 => D2DIAGNOS(:,pppuPIY3(iiP1): pppuPIY3(iiP6))
    puPIY3=ZERO
    puZEY3 => D2DIAGNOS(:,pppuZEY3(iiZ3): pppuZEY3(iiZ2))
    puZEY3=ZERO
    puP6Y3 => D2DIAGNOS(:,pppuP6Y3(iiYy3): pppuP6Y3(iiY3))
    puP6Y3=ZERO
    sediPI_Ben => D2DIAGNOS(:,ppsediPI_Ben(iiP1): ppsediPI_Ben(iiP6))
    sediPI_Ben=ZERO
    sediZE_Ben => D2DIAGNOS(:,ppsediZE_Ben(iiZ3): ppsediZE_Ben(iiZ2))
    sediZE_Ben=ZERO
    jnetPIc => D2DIAGNOS(:,ppjnetPIc(iiP1): ppjnetPIc(iiP6))
    jnetPIc=ZERO


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    PELSURFACE => D2DIAGNOS(:,357+1:357+56); PELSURFACE=ZERO
    jsurR9x => D2DIAGNOS(:,357+ppR9x); jsurR9x=ZERO
    jsurO2o => D2DIAGNOS(:,357+ppO2o); jsurO2o=ZERO
    jsurN1p => D2DIAGNOS(:,357+ppN1p); jsurN1p=ZERO
    jsurN3n => D2DIAGNOS(:,357+ppN3n); jsurN3n=ZERO
    jsurN4n => D2DIAGNOS(:,357+ppN4n); jsurN4n=ZERO
    jsurN5s => D2DIAGNOS(:,357+ppN5s); jsurN5s=ZERO
    jsurN6r => D2DIAGNOS(:,357+ppN6r); jsurN6r=ZERO
    jsurB1c => D2DIAGNOS(:,357+ppB1c); jsurB1c=ZERO
    jsurB1n => D2DIAGNOS(:,357+ppB1n); jsurB1n=ZERO
    jsurB1p => D2DIAGNOS(:,357+ppB1p); jsurB1p=ZERO
    jsurBac => D2DIAGNOS(:,357+ppBac); jsurBac=ZERO
    jsurP1c => D2DIAGNOS(:,357+ppP1c); jsurP1c=ZERO
    jsurP1n => D2DIAGNOS(:,357+ppP1n); jsurP1n=ZERO
    jsurP1p => D2DIAGNOS(:,357+ppP1p); jsurP1p=ZERO
    jsurP1l => D2DIAGNOS(:,357+ppP1l); jsurP1l=ZERO
    jsurP1s => D2DIAGNOS(:,357+ppP1s); jsurP1s=ZERO
    jsurP2c => D2DIAGNOS(:,357+ppP2c); jsurP2c=ZERO
    jsurP2n => D2DIAGNOS(:,357+ppP2n); jsurP2n=ZERO
    jsurP2p => D2DIAGNOS(:,357+ppP2p); jsurP2p=ZERO
    jsurP2l => D2DIAGNOS(:,357+ppP2l); jsurP2l=ZERO
    jsurP3c => D2DIAGNOS(:,357+ppP3c); jsurP3c=ZERO
    jsurP3n => D2DIAGNOS(:,357+ppP3n); jsurP3n=ZERO
    jsurP3p => D2DIAGNOS(:,357+ppP3p); jsurP3p=ZERO
    jsurP3l => D2DIAGNOS(:,357+ppP3l); jsurP3l=ZERO
    jsurP4c => D2DIAGNOS(:,357+ppP4c); jsurP4c=ZERO
    jsurP4n => D2DIAGNOS(:,357+ppP4n); jsurP4n=ZERO
    jsurP4p => D2DIAGNOS(:,357+ppP4p); jsurP4p=ZERO
    jsurP4l => D2DIAGNOS(:,357+ppP4l); jsurP4l=ZERO
    jsurP5c => D2DIAGNOS(:,357+ppP5c); jsurP5c=ZERO
    jsurP5n => D2DIAGNOS(:,357+ppP5n); jsurP5n=ZERO
    jsurP5p => D2DIAGNOS(:,357+ppP5p); jsurP5p=ZERO
    jsurP5l => D2DIAGNOS(:,357+ppP5l); jsurP5l=ZERO
    jsurP5s => D2DIAGNOS(:,357+ppP5s); jsurP5s=ZERO
    jsurP6c => D2DIAGNOS(:,357+ppP6c); jsurP6c=ZERO
    jsurP6n => D2DIAGNOS(:,357+ppP6n); jsurP6n=ZERO
    jsurP6p => D2DIAGNOS(:,357+ppP6p); jsurP6p=ZERO
    jsurP6l => D2DIAGNOS(:,357+ppP6l); jsurP6l=ZERO
    jsurPcc => D2DIAGNOS(:,357+ppPcc); jsurPcc=ZERO
    jsurR1c => D2DIAGNOS(:,357+ppR1c); jsurR1c=ZERO
    jsurR1n => D2DIAGNOS(:,357+ppR1n); jsurR1n=ZERO
    jsurR1p => D2DIAGNOS(:,357+ppR1p); jsurR1p=ZERO
    jsurR2c => D2DIAGNOS(:,357+ppR2c); jsurR2c=ZERO
    jsurR2n => D2DIAGNOS(:,357+ppR2n); jsurR2n=ZERO
    jsurR3c => D2DIAGNOS(:,357+ppR3c); jsurR3c=ZERO
    jsurR6c => D2DIAGNOS(:,357+ppR6c); jsurR6c=ZERO
    jsurR6n => D2DIAGNOS(:,357+ppR6n); jsurR6n=ZERO
    jsurR6p => D2DIAGNOS(:,357+ppR6p); jsurR6p=ZERO
    jsurR6s => D2DIAGNOS(:,357+ppR6s); jsurR6s=ZERO
    jsurRZc => D2DIAGNOS(:,357+ppRZc); jsurRZc=ZERO
    jsurO3c => D2DIAGNOS(:,357+ppO3c); jsurO3c=ZERO
    jsurO3h => D2DIAGNOS(:,357+ppO3h); jsurO3h=ZERO
    jsurZ3c => D2DIAGNOS(:,357+ppZ3c); jsurZ3c=ZERO
    jsurZ4c => D2DIAGNOS(:,357+ppZ4c); jsurZ4c=ZERO
    jsurZ2c => D2DIAGNOS(:,357+ppZ2c); jsurZ2c=ZERO
    jsurZ5c => D2DIAGNOS(:,357+ppZ5c); jsurZ5c=ZERO
    jsurZ6c => D2DIAGNOS(:,357+ppZ6c); jsurZ6c=ZERO

    PELBOTTOM => D2DIAGNOS(:,413+1:413+56); PELBOTTOM=ZERO
    jbotR9x => D2DIAGNOS(:,413+ppR9x); jbotR9x=ZERO
    jbotO2o => D2DIAGNOS(:,413+ppO2o); jbotO2o=ZERO
    jbotN1p => D2DIAGNOS(:,413+ppN1p); jbotN1p=ZERO
    jbotN3n => D2DIAGNOS(:,413+ppN3n); jbotN3n=ZERO
    jbotN4n => D2DIAGNOS(:,413+ppN4n); jbotN4n=ZERO
    jbotN5s => D2DIAGNOS(:,413+ppN5s); jbotN5s=ZERO
    jbotN6r => D2DIAGNOS(:,413+ppN6r); jbotN6r=ZERO
    jbotB1c => D2DIAGNOS(:,413+ppB1c); jbotB1c=ZERO
    jbotB1n => D2DIAGNOS(:,413+ppB1n); jbotB1n=ZERO
    jbotB1p => D2DIAGNOS(:,413+ppB1p); jbotB1p=ZERO
    jbotBac => D2DIAGNOS(:,413+ppBac); jbotBac=ZERO
    jbotP1c => D2DIAGNOS(:,413+ppP1c); jbotP1c=ZERO
    jbotP1n => D2DIAGNOS(:,413+ppP1n); jbotP1n=ZERO
    jbotP1p => D2DIAGNOS(:,413+ppP1p); jbotP1p=ZERO
    jbotP1l => D2DIAGNOS(:,413+ppP1l); jbotP1l=ZERO
    jbotP1s => D2DIAGNOS(:,413+ppP1s); jbotP1s=ZERO
    jbotP2c => D2DIAGNOS(:,413+ppP2c); jbotP2c=ZERO
    jbotP2n => D2DIAGNOS(:,413+ppP2n); jbotP2n=ZERO
    jbotP2p => D2DIAGNOS(:,413+ppP2p); jbotP2p=ZERO
    jbotP2l => D2DIAGNOS(:,413+ppP2l); jbotP2l=ZERO
    jbotP3c => D2DIAGNOS(:,413+ppP3c); jbotP3c=ZERO
    jbotP3n => D2DIAGNOS(:,413+ppP3n); jbotP3n=ZERO
    jbotP3p => D2DIAGNOS(:,413+ppP3p); jbotP3p=ZERO
    jbotP3l => D2DIAGNOS(:,413+ppP3l); jbotP3l=ZERO
    jbotP4c => D2DIAGNOS(:,413+ppP4c); jbotP4c=ZERO
    jbotP4n => D2DIAGNOS(:,413+ppP4n); jbotP4n=ZERO
    jbotP4p => D2DIAGNOS(:,413+ppP4p); jbotP4p=ZERO
    jbotP4l => D2DIAGNOS(:,413+ppP4l); jbotP4l=ZERO
    jbotP5c => D2DIAGNOS(:,413+ppP5c); jbotP5c=ZERO
    jbotP5n => D2DIAGNOS(:,413+ppP5n); jbotP5n=ZERO
    jbotP5p => D2DIAGNOS(:,413+ppP5p); jbotP5p=ZERO
    jbotP5l => D2DIAGNOS(:,413+ppP5l); jbotP5l=ZERO
    jbotP5s => D2DIAGNOS(:,413+ppP5s); jbotP5s=ZERO
    jbotP6c => D2DIAGNOS(:,413+ppP6c); jbotP6c=ZERO
    jbotP6n => D2DIAGNOS(:,413+ppP6n); jbotP6n=ZERO
    jbotP6p => D2DIAGNOS(:,413+ppP6p); jbotP6p=ZERO
    jbotP6l => D2DIAGNOS(:,413+ppP6l); jbotP6l=ZERO
    jbotPcc => D2DIAGNOS(:,413+ppPcc); jbotPcc=ZERO
    jbotR1c => D2DIAGNOS(:,413+ppR1c); jbotR1c=ZERO
    jbotR1n => D2DIAGNOS(:,413+ppR1n); jbotR1n=ZERO
    jbotR1p => D2DIAGNOS(:,413+ppR1p); jbotR1p=ZERO
    jbotR2c => D2DIAGNOS(:,413+ppR2c); jbotR2c=ZERO
    jbotR2n => D2DIAGNOS(:,413+ppR2n); jbotR2n=ZERO
    jbotR3c => D2DIAGNOS(:,413+ppR3c); jbotR3c=ZERO
    jbotR6c => D2DIAGNOS(:,413+ppR6c); jbotR6c=ZERO
    jbotR6n => D2DIAGNOS(:,413+ppR6n); jbotR6n=ZERO
    jbotR6p => D2DIAGNOS(:,413+ppR6p); jbotR6p=ZERO
    jbotR6s => D2DIAGNOS(:,413+ppR6s); jbotR6s=ZERO
    jbotRZc => D2DIAGNOS(:,413+ppRZc); jbotRZc=ZERO
    jbotO3c => D2DIAGNOS(:,413+ppO3c); jbotO3c=ZERO
    jbotO3h => D2DIAGNOS(:,413+ppO3h); jbotO3h=ZERO
    jbotZ3c => D2DIAGNOS(:,413+ppZ3c); jbotZ3c=ZERO
    jbotZ4c => D2DIAGNOS(:,413+ppZ4c); jbotZ4c=ZERO
    jbotZ2c => D2DIAGNOS(:,413+ppZ2c); jbotZ2c=ZERO
    jbotZ5c => D2DIAGNOS(:,413+ppZ5c); jbotZ5c=ZERO
    jbotZ6c => D2DIAGNOS(:,413+ppZ6c); jbotZ6c=ZERO

!    3d-state-field-alloc-pointer river
    allocate(iiPELSINKREF(1:NO_D3_BOX_STATES ),stat=status)
    if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem","iiPELSINKREF")
    iiPELSINKREF = 0





    allocate(OCDepth(1:NO_BOXES), stat=status); OCDepth = ZERO
    allocate(Depth(1:NO_BOXES), stat=status); Depth = ZERO
    allocate(ABIO_eps(1:NO_BOXES), stat=status); ABIO_eps = ZERO


    allocate(CO2_Ben(1:NO_BOXES_XY), stat=status); CO2_Ben = ZERO
    allocate(HCO3_Ben(1:NO_BOXES_XY), stat=status); HCO3_Ben = ZERO
    allocate(CO3_Ben(1:NO_BOXES_XY), stat=status); CO3_Ben = ZERO

    allocate(KPO4(1:NO_BOXES_XY), stat=status); KPO4 = 0
    allocate(KPO4sh(1:NO_BOXES_XY), stat=status); KPO4sh = 0
    allocate(KNH4p(1:NO_BOXES_XY), stat=status); KNH4p = 0
    allocate(KNH4(1:NO_BOXES_XY), stat=status); KNH4 = 0
    allocate(KNO3(1:NO_BOXES_XY), stat=status); KNO3 = 0
    allocate(KNO3E(1:NO_BOXES_XY), stat=status); KNO3E = 0
    allocate(KRED(1:NO_BOXES_XY), stat=status); KRED = 0
    allocate(KSIO3eq(1:NO_BOXES_XY), stat=status); KSIO3eq = 0
    allocate(KSIO3(1:NO_BOXES_XY), stat=status); KSIO3 = 0
    allocate(KQ1c(1:NO_BOXES_XY), stat=status); KQ1c = 0
    allocate(KQun(1:NO_BOXES_XY), stat=status); KQun = 0
    allocate(sw_CalcPhyto( 1:NO_BOXES_XY,1:iiPhytoPlankton),stat=status)
    sw_CalcPhyto = 0
    allocate(PelBoxAbove(1:NO_BOXES_XY), stat=status); PelBoxAbove = 0
    allocate(KCO2(1:NO_BOXES_XY), stat=status); KCO2 = 0
    allocate(KHPLUS(1:NO_BOXES_XY), stat=status); KHPLUS = 0

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of vars for calculation of combined fluxes for output
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    allocate(flx_calc_nr(0:58),stat=status)
    allocate(flx_CalcIn(1:58),stat=status)
    allocate(flx_option(1:58),stat=status)
    allocate(flx_t(1:440),stat=status)
    allocate(flx_SS(1:440),stat=status)
    allocate(flx_states(1:440),stat=status)
    allocate(flx_ostates(1:440),stat=status)
    flx_calc_nr(0)=0
    flx_cal_ben_start=22


    ! ruPTc=P.c <- O3c        (flux):
    flx_calc_nr(1)= 7; flx_CalcIn(1)=iiPel; flx_option(1)=0
    flx_t(1)=+1.00;flx_SS(1)=0; flx_states(1)=ppP1c;flx_ostates(1)=ppO3c
    flx_t(2)=+1.00;flx_SS(2)=0; flx_states(2)=ppP2c;flx_ostates(2)=ppO3c
    flx_t(3)=+1.00;flx_SS(3)=0; flx_states(3)=ppP3c;flx_ostates(3)=ppO3c
    flx_t(4)=+1.00;flx_SS(4)=0; flx_states(4)=ppP4c;flx_ostates(4)=ppO3c
    flx_t(5)=+1.00;flx_SS(5)=0; flx_states(5)=ppP5c;flx_ostates(5)=ppO3c
    flx_t(6)=+1.00;flx_SS(6)=0; flx_states(6)=ppP6c;flx_ostates(6)=ppO3c
    flx_t(7)=+1.00;flx_SS(7)=0; flx_states(7)=ppPcc;flx_ostates(7)=ppO3c

    ! ruPTn=P.n <- N3n+N4n        (flux):
    flx_calc_nr(2)= 19; flx_CalcIn(2)=iiPel; flx_option(2)=0
    flx_t(8)=+1.00;flx_SS(8)=0; flx_states(8)=ppP1n;flx_ostates(8)=ppN3n
    flx_t(9)=+1.00;flx_SS(9)=0; flx_states(9)=ppP1n;flx_ostates(9)=ppN4n
    flx_t(10)=+1.00;flx_SS(10)=0; flx_states(10)=ppP2n;flx_ostates(10)=ppN3n
    flx_t(11)=+1.00;flx_SS(11)=0; flx_states(11)=ppP2n;flx_ostates(11)=ppN4n
    flx_t(12)=+1.00;flx_SS(12)=0; flx_states(12)=ppP3n;flx_ostates(12)=ppN3n
    flx_t(13)=+1.00;flx_SS(13)=0; flx_states(13)=ppP3n;flx_ostates(13)=ppN4n
    flx_t(14)=+1.00;flx_SS(14)=0; flx_states(14)=ppP4n;flx_ostates(14)=ppN3n
    flx_t(15)=+1.00;flx_SS(15)=0; flx_states(15)=ppP4n;flx_ostates(15)=ppN4n
    flx_t(16)=+1.00;flx_SS(16)=0; flx_states(16)=ppP5n;flx_ostates(16)=ppN3n
    flx_t(17)=+1.00;flx_SS(17)=0; flx_states(17)=ppP5n;flx_ostates(17)=ppN4n
    flx_t(18)=+1.00;flx_SS(18)=0; flx_states(18)=ppP6n;flx_ostates(18)=ppN3n
    flx_t(19)=+1.00;flx_SS(19)=0; flx_states(19)=ppP6n;flx_ostates(19)=ppN4n

    ! ruPTP=P.p <- N1p        (flux):
    flx_calc_nr(3)= 25; flx_CalcIn(3)=iiPel; flx_option(3)=0
    flx_t(20)=+1.00;flx_SS(20)=0; flx_states(20)=ppP1p;flx_ostates(20)=ppN1p
    flx_t(21)=+1.00;flx_SS(21)=0; flx_states(21)=ppP2p;flx_ostates(21)=ppN1p
    flx_t(22)=+1.00;flx_SS(22)=0; flx_states(22)=ppP3p;flx_ostates(22)=ppN1p
    flx_t(23)=+1.00;flx_SS(23)=0; flx_states(23)=ppP4p;flx_ostates(23)=ppN1p
    flx_t(24)=+1.00;flx_SS(24)=0; flx_states(24)=ppP5p;flx_ostates(24)=ppN1p
    flx_t(25)=+1.00;flx_SS(25)=0; flx_states(25)=ppP6p;flx_ostates(25)=ppN1p

    ! exPP=(P.c->R1c+R2c+R6c)        (flux):
    flx_calc_nr(4)= 46; flx_CalcIn(4)=iiPel; flx_option(4)=0
    flx_t(26)=+1.00;flx_SS(26)=1; flx_states(26)=ppP1c;flx_ostates(26)=ppR1c
    flx_t(27)=+1.00;flx_SS(27)=1; flx_states(27)=ppP1c;flx_ostates(27)=ppR2c
    flx_t(28)=+1.00;flx_SS(28)=1; flx_states(28)=ppP1c;flx_ostates(28)=ppR6c
    flx_t(29)=+1.00;flx_SS(29)=1; flx_states(29)=ppP2c;flx_ostates(29)=ppR1c
    flx_t(30)=+1.00;flx_SS(30)=1; flx_states(30)=ppP2c;flx_ostates(30)=ppR2c
    flx_t(31)=+1.00;flx_SS(31)=1; flx_states(31)=ppP2c;flx_ostates(31)=ppR6c
    flx_t(32)=+1.00;flx_SS(32)=1; flx_states(32)=ppP3c;flx_ostates(32)=ppR1c
    flx_t(33)=+1.00;flx_SS(33)=1; flx_states(33)=ppP3c;flx_ostates(33)=ppR2c
    flx_t(34)=+1.00;flx_SS(34)=1; flx_states(34)=ppP3c;flx_ostates(34)=ppR6c
    flx_t(35)=+1.00;flx_SS(35)=1; flx_states(35)=ppP4c;flx_ostates(35)=ppR1c
    flx_t(36)=+1.00;flx_SS(36)=1; flx_states(36)=ppP4c;flx_ostates(36)=ppR2c
    flx_t(37)=+1.00;flx_SS(37)=1; flx_states(37)=ppP4c;flx_ostates(37)=ppR6c
    flx_t(38)=+1.00;flx_SS(38)=1; flx_states(38)=ppP5c;flx_ostates(38)=ppR1c
    flx_t(39)=+1.00;flx_SS(39)=1; flx_states(39)=ppP5c;flx_ostates(39)=ppR2c
    flx_t(40)=+1.00;flx_SS(40)=1; flx_states(40)=ppP5c;flx_ostates(40)=ppR6c
    flx_t(41)=+1.00;flx_SS(41)=1; flx_states(41)=ppP6c;flx_ostates(41)=ppR1c
    flx_t(42)=+1.00;flx_SS(42)=1; flx_states(42)=ppP6c;flx_ostates(42)=ppR2c
    flx_t(43)=+1.00;flx_SS(43)=1; flx_states(43)=ppP6c;flx_ostates(43)=ppR6c
    flx_t(44)=+1.00;flx_SS(44)=1; flx_states(44)=ppPcc;flx_ostates(44)=ppR1c
    flx_t(45)=+1.00;flx_SS(45)=1; flx_states(45)=ppPcc;flx_ostates(45)=ppR2c
    flx_t(46)=+1.00;flx_SS(46)=1; flx_states(46)=ppPcc;flx_ostates(46)=ppR6c

    ! resPP=(P.c->O3c)        (flux):
    flx_calc_nr(5)= 53; flx_CalcIn(5)=iiPel; flx_option(5)=0
    flx_t(47)=+1.00;flx_SS(47)=1; flx_states(47)=ppP1c;flx_ostates(47)=ppO3c
    flx_t(48)=+1.00;flx_SS(48)=1; flx_states(48)=ppP2c;flx_ostates(48)=ppO3c
    flx_t(49)=+1.00;flx_SS(49)=1; flx_states(49)=ppP3c;flx_ostates(49)=ppO3c
    flx_t(50)=+1.00;flx_SS(50)=1; flx_states(50)=ppP4c;flx_ostates(50)=ppO3c
    flx_t(51)=+1.00;flx_SS(51)=1; flx_states(51)=ppP5c;flx_ostates(51)=ppO3c
    flx_t(52)=+1.00;flx_SS(52)=1; flx_states(52)=ppP6c;flx_ostates(52)=ppO3c
    flx_t(53)=+1.00;flx_SS(53)=1; flx_states(53)=ppPcc;flx_ostates(53)=ppO3c

    ! ruZTc=(Z.c<-P.c+B1c+Z.c)-(Z.c->R1c+R6c)        (flux):
    flx_calc_nr(6)= 123; flx_CalcIn(6)=iiPel; flx_option(6)=0
    flx_t(54)=-1.00;flx_SS(54)=1; flx_states(54)=ppZ3c;flx_ostates(54)=ppR1c
    flx_t(55)=+1.00;flx_SS(55)=0; flx_states(55)=ppZ3c;flx_ostates(55)=ppB1c
    flx_t(56)=+1.00;flx_SS(56)=0; flx_states(56)=ppZ3c;flx_ostates(56)=ppP1c
    flx_t(57)=+1.00;flx_SS(57)=0; flx_states(57)=ppZ3c;flx_ostates(57)=ppP2c
    flx_t(58)=+1.00;flx_SS(58)=0; flx_states(58)=ppZ3c;flx_ostates(58)=ppP3c
    flx_t(59)=+1.00;flx_SS(59)=0; flx_states(59)=ppZ3c;flx_ostates(59)=ppP4c
    flx_t(60)=+1.00;flx_SS(60)=0; flx_states(60)=ppZ3c;flx_ostates(60)=ppP5c
    flx_t(61)=+1.00;flx_SS(61)=0; flx_states(61)=ppZ3c;flx_ostates(61)=ppP6c
    flx_t(62)=+1.00;flx_SS(62)=0; flx_states(62)=ppZ3c;flx_ostates(62)=ppPcc
    flx_t(63)=+1.00;flx_SS(63)=0; flx_states(63)=ppZ3c;flx_ostates(63)=ppZ3c
    flx_t(64)=+1.00;flx_SS(64)=0; flx_states(64)=ppZ3c;flx_ostates(64)=ppZ4c
    flx_t(65)=+1.00;flx_SS(65)=0; flx_states(65)=ppZ3c;flx_ostates(65)=ppZ2c
    flx_t(66)=+1.00;flx_SS(66)=0; flx_states(66)=ppZ3c;flx_ostates(66)=ppZ5c
    flx_t(67)=+1.00;flx_SS(67)=0; flx_states(67)=ppZ3c;flx_ostates(67)=ppZ6c
    flx_t(68)=-1.00;flx_SS(68)=1; flx_states(68)=ppZ4c;flx_ostates(68)=ppR1c
    flx_t(69)=+1.00;flx_SS(69)=0; flx_states(69)=ppZ4c;flx_ostates(69)=ppB1c
    flx_t(70)=+1.00;flx_SS(70)=0; flx_states(70)=ppZ4c;flx_ostates(70)=ppP1c
    flx_t(71)=+1.00;flx_SS(71)=0; flx_states(71)=ppZ4c;flx_ostates(71)=ppP2c
    flx_t(72)=+1.00;flx_SS(72)=0; flx_states(72)=ppZ4c;flx_ostates(72)=ppP3c
    flx_t(73)=+1.00;flx_SS(73)=0; flx_states(73)=ppZ4c;flx_ostates(73)=ppP4c
    flx_t(74)=+1.00;flx_SS(74)=0; flx_states(74)=ppZ4c;flx_ostates(74)=ppP5c
    flx_t(75)=+1.00;flx_SS(75)=0; flx_states(75)=ppZ4c;flx_ostates(75)=ppP6c
    flx_t(76)=+1.00;flx_SS(76)=0; flx_states(76)=ppZ4c;flx_ostates(76)=ppPcc
    flx_t(77)=+1.00;flx_SS(77)=0; flx_states(77)=ppZ4c;flx_ostates(77)=ppZ3c
    flx_t(78)=+1.00;flx_SS(78)=0; flx_states(78)=ppZ4c;flx_ostates(78)=ppZ4c
    flx_t(79)=+1.00;flx_SS(79)=0; flx_states(79)=ppZ4c;flx_ostates(79)=ppZ2c
    flx_t(80)=+1.00;flx_SS(80)=0; flx_states(80)=ppZ4c;flx_ostates(80)=ppZ5c
    flx_t(81)=+1.00;flx_SS(81)=0; flx_states(81)=ppZ4c;flx_ostates(81)=ppZ6c
    flx_t(82)=-1.00;flx_SS(82)=1; flx_states(82)=ppZ2c;flx_ostates(82)=ppR1c
    flx_t(83)=+1.00;flx_SS(83)=0; flx_states(83)=ppZ2c;flx_ostates(83)=ppB1c
    flx_t(84)=+1.00;flx_SS(84)=0; flx_states(84)=ppZ2c;flx_ostates(84)=ppP1c
    flx_t(85)=+1.00;flx_SS(85)=0; flx_states(85)=ppZ2c;flx_ostates(85)=ppP2c
    flx_t(86)=+1.00;flx_SS(86)=0; flx_states(86)=ppZ2c;flx_ostates(86)=ppP3c
    flx_t(87)=+1.00;flx_SS(87)=0; flx_states(87)=ppZ2c;flx_ostates(87)=ppP4c
    flx_t(88)=+1.00;flx_SS(88)=0; flx_states(88)=ppZ2c;flx_ostates(88)=ppP5c
    flx_t(89)=+1.00;flx_SS(89)=0; flx_states(89)=ppZ2c;flx_ostates(89)=ppP6c
    flx_t(90)=+1.00;flx_SS(90)=0; flx_states(90)=ppZ2c;flx_ostates(90)=ppPcc
    flx_t(91)=+1.00;flx_SS(91)=0; flx_states(91)=ppZ2c;flx_ostates(91)=ppZ3c
    flx_t(92)=+1.00;flx_SS(92)=0; flx_states(92)=ppZ2c;flx_ostates(92)=ppZ4c
    flx_t(93)=+1.00;flx_SS(93)=0; flx_states(93)=ppZ2c;flx_ostates(93)=ppZ2c
    flx_t(94)=+1.00;flx_SS(94)=0; flx_states(94)=ppZ2c;flx_ostates(94)=ppZ5c
    flx_t(95)=+1.00;flx_SS(95)=0; flx_states(95)=ppZ2c;flx_ostates(95)=ppZ6c
    flx_t(96)=-1.00;flx_SS(96)=1; flx_states(96)=ppZ5c;flx_ostates(96)=ppR1c
    flx_t(97)=+1.00;flx_SS(97)=0; flx_states(97)=ppZ5c;flx_ostates(97)=ppB1c
    flx_t(98)=+1.00;flx_SS(98)=0; flx_states(98)=ppZ5c;flx_ostates(98)=ppP1c
    flx_t(99)=+1.00;flx_SS(99)=0; flx_states(99)=ppZ5c;flx_ostates(99)=ppP2c
    flx_t(100)=+1.00;flx_SS(100)=0; flx_states(100)=ppZ5c
    flx_ostates(100)=ppP3c
    flx_t(101)=+1.00;flx_SS(101)=0; flx_states(101)=ppZ5c
    flx_ostates(101)=ppP4c
    flx_t(102)=+1.00;flx_SS(102)=0; flx_states(102)=ppZ5c
    flx_ostates(102)=ppP5c
    flx_t(103)=+1.00;flx_SS(103)=0; flx_states(103)=ppZ5c
    flx_ostates(103)=ppP6c
    flx_t(104)=+1.00;flx_SS(104)=0; flx_states(104)=ppZ5c
    flx_ostates(104)=ppPcc
    flx_t(105)=+1.00;flx_SS(105)=0; flx_states(105)=ppZ5c
    flx_ostates(105)=ppZ3c
    flx_t(106)=+1.00;flx_SS(106)=0; flx_states(106)=ppZ5c
    flx_ostates(106)=ppZ4c
    flx_t(107)=+1.00;flx_SS(107)=0; flx_states(107)=ppZ5c
    flx_ostates(107)=ppZ2c
    flx_t(108)=+1.00;flx_SS(108)=0; flx_states(108)=ppZ5c
    flx_ostates(108)=ppZ5c
    flx_t(109)=+1.00;flx_SS(109)=0; flx_states(109)=ppZ5c
    flx_ostates(109)=ppZ6c
    flx_t(110)=-1.00;flx_SS(110)=1; flx_states(110)=ppZ6c
    flx_ostates(110)=ppR1c
    flx_t(111)=+1.00;flx_SS(111)=0; flx_states(111)=ppZ6c
    flx_ostates(111)=ppB1c
    flx_t(112)=+1.00;flx_SS(112)=0; flx_states(112)=ppZ6c
    flx_ostates(112)=ppP1c
    flx_t(113)=+1.00;flx_SS(113)=0; flx_states(113)=ppZ6c
    flx_ostates(113)=ppP2c
    flx_t(114)=+1.00;flx_SS(114)=0; flx_states(114)=ppZ6c
    flx_ostates(114)=ppP3c
    flx_t(115)=+1.00;flx_SS(115)=0; flx_states(115)=ppZ6c
    flx_ostates(115)=ppP4c
    flx_t(116)=+1.00;flx_SS(116)=0; flx_states(116)=ppZ6c
    flx_ostates(116)=ppP5c
    flx_t(117)=+1.00;flx_SS(117)=0; flx_states(117)=ppZ6c
    flx_ostates(117)=ppP6c
    flx_t(118)=+1.00;flx_SS(118)=0; flx_states(118)=ppZ6c
    flx_ostates(118)=ppPcc
    flx_t(119)=+1.00;flx_SS(119)=0; flx_states(119)=ppZ6c
    flx_ostates(119)=ppZ3c
    flx_t(120)=+1.00;flx_SS(120)=0; flx_states(120)=ppZ6c
    flx_ostates(120)=ppZ4c
    flx_t(121)=+1.00;flx_SS(121)=0; flx_states(121)=ppZ6c
    flx_ostates(121)=ppZ2c
    flx_t(122)=+1.00;flx_SS(122)=0; flx_states(122)=ppZ6c
    flx_ostates(122)=ppZ5c
    flx_t(123)=+1.00;flx_SS(123)=0; flx_states(123)=ppZ6c
    flx_ostates(123)=ppZ6c

    ! rrPTo=(O2o->*)        (flux):
    flx_calc_nr(7)= 124; flx_CalcIn(7)=iiPel; flx_option(7)=0
    flx_t(124)=+1.00;flx_SS(124)=1; flx_states(124)=ppO2o
    flx_ostates(124)=ppO2o

    ! flR2B1c=R2c->B1c        (flux):
    flx_calc_nr(8)= 125; flx_CalcIn(8)=iiPel; flx_option(8)=0
    flx_t(125)=+1.00;flx_SS(125)=1; flx_states(125)=ppR2c
    flx_ostates(125)=ppB1c

    ! flP1R2c=P1c->R2c        (flux):
    flx_calc_nr(9)= 126; flx_CalcIn(9)=iiPel; flx_option(9)=0
    flx_t(126)=+1.00;flx_SS(126)=1; flx_states(126)=ppP1c
    flx_ostates(126)=ppR2c

    ! flP6R3c=P6c->R3c        (flux):
    flx_calc_nr(10)= 127; flx_CalcIn(10)=iiPel; flx_option(10)=0
    flx_t(127)=+1.00;flx_SS(127)=1; flx_states(127)=ppP6c
    flx_ostates(127)=ppR3c

    ! qqutZIc=Z.c        (u):
    flx_calc_nr(11)= 132; flx_CalcIn(11)=iiPel; flx_option(11)=30
    flx_t(128)=+1.00;flx_SS(128)=1; flx_states(128)=ppZ3c;flx_ostates(128)=0
    flx_t(129)=+1.00;flx_SS(129)=1; flx_states(129)=ppZ4c;flx_ostates(129)=0
    flx_t(130)=+1.00;flx_SS(130)=1; flx_states(130)=ppZ2c;flx_ostates(130)=0
    flx_t(131)=+1.00;flx_SS(131)=1; flx_states(131)=ppZ5c;flx_ostates(131)=0
    flx_t(132)=+1.00;flx_SS(132)=1; flx_states(132)=ppZ6c;flx_ostates(132)=0

    ! qqutPIc=P.c        (u):
    flx_calc_nr(12)= 139; flx_CalcIn(12)=iiPel; flx_option(12)=30
    flx_t(133)=+1.00;flx_SS(133)=1; flx_states(133)=ppP1c;flx_ostates(133)=0
    flx_t(134)=+1.00;flx_SS(134)=1; flx_states(134)=ppP2c;flx_ostates(134)=0
    flx_t(135)=+1.00;flx_SS(135)=1; flx_states(135)=ppP3c;flx_ostates(135)=0
    flx_t(136)=+1.00;flx_SS(136)=1; flx_states(136)=ppP4c;flx_ostates(136)=0
    flx_t(137)=+1.00;flx_SS(137)=1; flx_states(137)=ppP5c;flx_ostates(137)=0
    flx_t(138)=+1.00;flx_SS(138)=1; flx_states(138)=ppP6c;flx_ostates(138)=0
    flx_t(139)=+1.00;flx_SS(139)=1; flx_states(139)=ppPcc;flx_ostates(139)=0

    ! qqutRc=R1c+R2c+R3c+R6c+RZc        (u):
    flx_calc_nr(13)= 144; flx_CalcIn(13)=iiPel; flx_option(13)=30
    flx_t(140)=+1.00;flx_SS(140)=1; flx_states(140)=ppR1c;flx_ostates(140)=0
    flx_t(141)=+1.00;flx_SS(141)=1; flx_states(141)=ppR2c;flx_ostates(141)=0
    flx_t(142)=+1.00;flx_SS(142)=1; flx_states(142)=ppR3c;flx_ostates(142)=0
    flx_t(143)=+1.00;flx_SS(143)=1; flx_states(143)=ppR6c;flx_ostates(143)=0
    flx_t(144)=+1.00;flx_SS(144)=1; flx_states(144)=ppRZc;flx_ostates(144)=0

    ! qqutDIc=O3c        (u):
    flx_calc_nr(14)= 145; flx_CalcIn(14)=iiPel; flx_option(14)=30
    flx_t(145)=+1.00;flx_SS(145)=1; flx_states(145)=ppO3c;flx_ostates(145)=0

    ! qqutBc=B1c+Bac        (u):
    flx_calc_nr(15)= 147; flx_CalcIn(15)=iiPel; flx_option(15)=30
    flx_t(146)=+1.00;flx_SS(146)=1; flx_states(146)=ppB1c;flx_ostates(146)=0
    flx_t(147)=+1.00;flx_SS(147)=1; flx_states(147)=ppBac;flx_ostates(147)=0

    ! qqut_w=1        (u):
    flx_calc_nr(16)= 148; flx_CalcIn(16)=iiPel; flx_option(16)=30
    flx_t(148)=+1.00;flx_SS(148)=1; flx_states(148)=0;flx_ostates(148)=0

    ! qqvtZIc=Z.c        (v):
    flx_calc_nr(17)= 153; flx_CalcIn(17)=iiPel; flx_option(17)=40
    flx_t(149)=+1.00;flx_SS(149)=1; flx_states(149)=ppZ3c;flx_ostates(149)=0
    flx_t(150)=+1.00;flx_SS(150)=1; flx_states(150)=ppZ4c;flx_ostates(150)=0
    flx_t(151)=+1.00;flx_SS(151)=1; flx_states(151)=ppZ2c;flx_ostates(151)=0
    flx_t(152)=+1.00;flx_SS(152)=1; flx_states(152)=ppZ5c;flx_ostates(152)=0
    flx_t(153)=+1.00;flx_SS(153)=1; flx_states(153)=ppZ6c;flx_ostates(153)=0

    ! qqvtPIc=P.c        (v):
    flx_calc_nr(18)= 160; flx_CalcIn(18)=iiPel; flx_option(18)=40
    flx_t(154)=+1.00;flx_SS(154)=1; flx_states(154)=ppP1c;flx_ostates(154)=0
    flx_t(155)=+1.00;flx_SS(155)=1; flx_states(155)=ppP2c;flx_ostates(155)=0
    flx_t(156)=+1.00;flx_SS(156)=1; flx_states(156)=ppP3c;flx_ostates(156)=0
    flx_t(157)=+1.00;flx_SS(157)=1; flx_states(157)=ppP4c;flx_ostates(157)=0
    flx_t(158)=+1.00;flx_SS(158)=1; flx_states(158)=ppP5c;flx_ostates(158)=0
    flx_t(159)=+1.00;flx_SS(159)=1; flx_states(159)=ppP6c;flx_ostates(159)=0
    flx_t(160)=+1.00;flx_SS(160)=1; flx_states(160)=ppPcc;flx_ostates(160)=0

    ! qqvtRc=R1c+R2c+R3c+R6c+RZc        (v):
    flx_calc_nr(19)= 165; flx_CalcIn(19)=iiPel; flx_option(19)=40
    flx_t(161)=+1.00;flx_SS(161)=1; flx_states(161)=ppR1c;flx_ostates(161)=0
    flx_t(162)=+1.00;flx_SS(162)=1; flx_states(162)=ppR2c;flx_ostates(162)=0
    flx_t(163)=+1.00;flx_SS(163)=1; flx_states(163)=ppR3c;flx_ostates(163)=0
    flx_t(164)=+1.00;flx_SS(164)=1; flx_states(164)=ppR6c;flx_ostates(164)=0
    flx_t(165)=+1.00;flx_SS(165)=1; flx_states(165)=ppRZc;flx_ostates(165)=0

    ! qqvtDIc=O3c        (v):
    flx_calc_nr(20)= 166; flx_CalcIn(20)=iiPel; flx_option(20)=40
    flx_t(166)=+1.00;flx_SS(166)=1; flx_states(166)=ppO3c;flx_ostates(166)=0

    ! qqvtBc=B1c+Bac        (v):
    flx_calc_nr(21)= 168; flx_CalcIn(21)=iiPel; flx_option(21)=40
    flx_t(167)=+1.00;flx_SS(167)=1; flx_states(167)=ppB1c;flx_ostates(167)=0
    flx_t(168)=+1.00;flx_SS(168)=1; flx_states(168)=ppBac;flx_ostates(168)=0

    ! qqvt_w=1        (v):
    flx_calc_nr(22)= 169; flx_CalcIn(22)=iiPel; flx_option(22)=40
    flx_t(169)=+1.00;flx_SS(169)=1; flx_states(169)=0;flx_ostates(169)=0


    ! grsPPm2=(P.c<-O3c)        (flux perm2):
    flx_calc_nr(23)= 176; flx_CalcIn(23)=iiPel; flx_option(23)=2
    flx_t(170)=+1.00;flx_SS(170)=0; flx_states(170)=ppP1c
    flx_ostates(170)=ppO3c
    flx_t(171)=+1.00;flx_SS(171)=0; flx_states(171)=ppP2c
    flx_ostates(171)=ppO3c
    flx_t(172)=+1.00;flx_SS(172)=0; flx_states(172)=ppP3c
    flx_ostates(172)=ppO3c
    flx_t(173)=+1.00;flx_SS(173)=0; flx_states(173)=ppP4c
    flx_ostates(173)=ppO3c
    flx_t(174)=+1.00;flx_SS(174)=0; flx_states(174)=ppP5c
    flx_ostates(174)=ppO3c
    flx_t(175)=+1.00;flx_SS(175)=0; flx_states(175)=ppP6c
    flx_ostates(175)=ppO3c
    flx_t(176)=+1.00;flx_SS(176)=0; flx_states(176)=ppPcc
    flx_ostates(176)=ppO3c

    ! netPPm2=(P.c<-O3c)-(P.c->O3c)        (flux perm2):
    flx_calc_nr(24)= 190; flx_CalcIn(24)=iiPel; flx_option(24)=2
    flx_t(177)=+1.00;flx_SS(177)=0; flx_states(177)=ppP1c
    flx_ostates(177)=ppO3c
    flx_t(178)=-1.00;flx_SS(178)=1; flx_states(178)=ppP1c
    flx_ostates(178)=ppO3c
    flx_t(179)=+1.00;flx_SS(179)=0; flx_states(179)=ppP2c
    flx_ostates(179)=ppO3c
    flx_t(180)=-1.00;flx_SS(180)=1; flx_states(180)=ppP2c
    flx_ostates(180)=ppO3c
    flx_t(181)=+1.00;flx_SS(181)=0; flx_states(181)=ppP3c
    flx_ostates(181)=ppO3c
    flx_t(182)=-1.00;flx_SS(182)=1; flx_states(182)=ppP3c
    flx_ostates(182)=ppO3c
    flx_t(183)=+1.00;flx_SS(183)=0; flx_states(183)=ppP4c
    flx_ostates(183)=ppO3c
    flx_t(184)=-1.00;flx_SS(184)=1; flx_states(184)=ppP4c
    flx_ostates(184)=ppO3c
    flx_t(185)=+1.00;flx_SS(185)=0; flx_states(185)=ppP5c
    flx_ostates(185)=ppO3c
    flx_t(186)=-1.00;flx_SS(186)=1; flx_states(186)=ppP5c
    flx_ostates(186)=ppO3c
    flx_t(187)=+1.00;flx_SS(187)=0; flx_states(187)=ppP6c
    flx_ostates(187)=ppO3c
    flx_t(188)=-1.00;flx_SS(188)=1; flx_states(188)=ppP6c
    flx_ostates(188)=ppO3c
    flx_t(189)=+1.00;flx_SS(189)=0; flx_states(189)=ppPcc
    flx_ostates(189)=ppO3c
    flx_t(190)=-1.00;flx_SS(190)=1; flx_states(190)=ppPcc
    flx_ostates(190)=ppO3c

    ! jnetR6c=(R6c<-P.c+Z.c+B1c)-(R6c->B1c)        (flux perm2):
    flx_calc_nr(25)= 204; flx_CalcIn(25)=iiPel; flx_option(25)=2
    flx_t(191)=-1.00;flx_SS(191)=1; flx_states(191)=ppR6c
    flx_ostates(191)=ppB1c
    flx_t(192)=+1.00;flx_SS(192)=0; flx_states(192)=ppR6c
    flx_ostates(192)=ppB1c
    flx_t(193)=+1.00;flx_SS(193)=0; flx_states(193)=ppR6c
    flx_ostates(193)=ppP1c
    flx_t(194)=+1.00;flx_SS(194)=0; flx_states(194)=ppR6c
    flx_ostates(194)=ppP2c
    flx_t(195)=+1.00;flx_SS(195)=0; flx_states(195)=ppR6c
    flx_ostates(195)=ppP3c
    flx_t(196)=+1.00;flx_SS(196)=0; flx_states(196)=ppR6c
    flx_ostates(196)=ppP4c
    flx_t(197)=+1.00;flx_SS(197)=0; flx_states(197)=ppR6c
    flx_ostates(197)=ppP5c
    flx_t(198)=+1.00;flx_SS(198)=0; flx_states(198)=ppR6c
    flx_ostates(198)=ppP6c
    flx_t(199)=+1.00;flx_SS(199)=0; flx_states(199)=ppR6c
    flx_ostates(199)=ppPcc
    flx_t(200)=+1.00;flx_SS(200)=0; flx_states(200)=ppR6c
    flx_ostates(200)=ppZ3c
    flx_t(201)=+1.00;flx_SS(201)=0; flx_states(201)=ppR6c
    flx_ostates(201)=ppZ4c
    flx_t(202)=+1.00;flx_SS(202)=0; flx_states(202)=ppR6c
    flx_ostates(202)=ppZ2c
    flx_t(203)=+1.00;flx_SS(203)=0; flx_states(203)=ppR6c
    flx_ostates(203)=ppZ5c
    flx_t(204)=+1.00;flx_SS(204)=0; flx_states(204)=ppR6c
    flx_ostates(204)=ppZ6c

    ! jnetR6n=(R6n<-P.n+B1n)-(R6n->B1n)        (flux perm2):
    flx_calc_nr(26)= 212; flx_CalcIn(26)=iiPel; flx_option(26)=2
    flx_t(205)=-1.00;flx_SS(205)=1; flx_states(205)=ppR6n
    flx_ostates(205)=ppB1n
    flx_t(206)=+1.00;flx_SS(206)=0; flx_states(206)=ppR6n
    flx_ostates(206)=ppB1n
    flx_t(207)=+1.00;flx_SS(207)=0; flx_states(207)=ppR6n
    flx_ostates(207)=ppP1n
    flx_t(208)=+1.00;flx_SS(208)=0; flx_states(208)=ppR6n
    flx_ostates(208)=ppP2n
    flx_t(209)=+1.00;flx_SS(209)=0; flx_states(209)=ppR6n
    flx_ostates(209)=ppP3n
    flx_t(210)=+1.00;flx_SS(210)=0; flx_states(210)=ppR6n
    flx_ostates(210)=ppP4n
    flx_t(211)=+1.00;flx_SS(211)=0; flx_states(211)=ppR6n
    flx_ostates(211)=ppP5n
    flx_t(212)=+1.00;flx_SS(212)=0; flx_states(212)=ppR6n
    flx_ostates(212)=ppP6n

    ! jnetR6p=(R6p<-P.p+B1p)-(R6p->B1p)        (flux perm2):
    flx_calc_nr(27)= 220; flx_CalcIn(27)=iiPel; flx_option(27)=2
    flx_t(213)=-1.00;flx_SS(213)=1; flx_states(213)=ppR6p
    flx_ostates(213)=ppB1p
    flx_t(214)=+1.00;flx_SS(214)=0; flx_states(214)=ppR6p
    flx_ostates(214)=ppB1p
    flx_t(215)=+1.00;flx_SS(215)=0; flx_states(215)=ppR6p
    flx_ostates(215)=ppP1p
    flx_t(216)=+1.00;flx_SS(216)=0; flx_states(216)=ppR6p
    flx_ostates(216)=ppP2p
    flx_t(217)=+1.00;flx_SS(217)=0; flx_states(217)=ppR6p
    flx_ostates(217)=ppP3p
    flx_t(218)=+1.00;flx_SS(218)=0; flx_states(218)=ppR6p
    flx_ostates(218)=ppP4p
    flx_t(219)=+1.00;flx_SS(219)=0; flx_states(219)=ppR6p
    flx_ostates(219)=ppP5p
    flx_t(220)=+1.00;flx_SS(220)=0; flx_states(220)=ppR6p
    flx_ostates(220)=ppP6p

    ! jnetR2c=(R2c<-P.c+Z.c)-(R2c->B1c)        (flux perm2):
    flx_calc_nr(28)= 233; flx_CalcIn(28)=iiPel; flx_option(28)=2
    flx_t(221)=-1.00;flx_SS(221)=1; flx_states(221)=ppR2c
    flx_ostates(221)=ppB1c
    flx_t(222)=+1.00;flx_SS(222)=0; flx_states(222)=ppR2c
    flx_ostates(222)=ppP1c
    flx_t(223)=+1.00;flx_SS(223)=0; flx_states(223)=ppR2c
    flx_ostates(223)=ppP2c
    flx_t(224)=+1.00;flx_SS(224)=0; flx_states(224)=ppR2c
    flx_ostates(224)=ppP3c
    flx_t(225)=+1.00;flx_SS(225)=0; flx_states(225)=ppR2c
    flx_ostates(225)=ppP4c
    flx_t(226)=+1.00;flx_SS(226)=0; flx_states(226)=ppR2c
    flx_ostates(226)=ppP5c
    flx_t(227)=+1.00;flx_SS(227)=0; flx_states(227)=ppR2c
    flx_ostates(227)=ppP6c
    flx_t(228)=+1.00;flx_SS(228)=0; flx_states(228)=ppR2c
    flx_ostates(228)=ppPcc
    flx_t(229)=+1.00;flx_SS(229)=0; flx_states(229)=ppR2c
    flx_ostates(229)=ppZ3c
    flx_t(230)=+1.00;flx_SS(230)=0; flx_states(230)=ppR2c
    flx_ostates(230)=ppZ4c
    flx_t(231)=+1.00;flx_SS(231)=0; flx_states(231)=ppR2c
    flx_ostates(231)=ppZ2c
    flx_t(232)=+1.00;flx_SS(232)=0; flx_states(232)=ppR2c
    flx_ostates(232)=ppZ5c
    flx_t(233)=+1.00;flx_SS(233)=0; flx_states(233)=ppR2c
    flx_ostates(233)=ppZ6c

    ! jnetR1c=(R1c<-P.c+Z.c+B1c)-(R1c->B1c)        (flux perm2):
    flx_calc_nr(29)= 247; flx_CalcIn(29)=iiPel; flx_option(29)=2
    flx_t(234)=-1.00;flx_SS(234)=1; flx_states(234)=ppR1c
    flx_ostates(234)=ppB1c
    flx_t(235)=+1.00;flx_SS(235)=0; flx_states(235)=ppR1c
    flx_ostates(235)=ppB1c
    flx_t(236)=+1.00;flx_SS(236)=0; flx_states(236)=ppR1c
    flx_ostates(236)=ppP1c
    flx_t(237)=+1.00;flx_SS(237)=0; flx_states(237)=ppR1c
    flx_ostates(237)=ppP2c
    flx_t(238)=+1.00;flx_SS(238)=0; flx_states(238)=ppR1c
    flx_ostates(238)=ppP3c
    flx_t(239)=+1.00;flx_SS(239)=0; flx_states(239)=ppR1c
    flx_ostates(239)=ppP4c
    flx_t(240)=+1.00;flx_SS(240)=0; flx_states(240)=ppR1c
    flx_ostates(240)=ppP5c
    flx_t(241)=+1.00;flx_SS(241)=0; flx_states(241)=ppR1c
    flx_ostates(241)=ppP6c
    flx_t(242)=+1.00;flx_SS(242)=0; flx_states(242)=ppR1c
    flx_ostates(242)=ppPcc
    flx_t(243)=+1.00;flx_SS(243)=0; flx_states(243)=ppR1c
    flx_ostates(243)=ppZ3c
    flx_t(244)=+1.00;flx_SS(244)=0; flx_states(244)=ppR1c
    flx_ostates(244)=ppZ4c
    flx_t(245)=+1.00;flx_SS(245)=0; flx_states(245)=ppR1c
    flx_ostates(245)=ppZ2c
    flx_t(246)=+1.00;flx_SS(246)=0; flx_states(246)=ppR1c
    flx_ostates(246)=ppZ5c
    flx_t(247)=+1.00;flx_SS(247)=0; flx_states(247)=ppR1c
    flx_ostates(247)=ppZ6c

    ! jnetR1n=(R1n<-P.n+B1n)-(R1n->B1n)        (flux perm2):
    flx_calc_nr(30)= 255; flx_CalcIn(30)=iiPel; flx_option(30)=2
    flx_t(248)=-1.00;flx_SS(248)=1; flx_states(248)=ppR1n
    flx_ostates(248)=ppB1n
    flx_t(249)=+1.00;flx_SS(249)=0; flx_states(249)=ppR1n
    flx_ostates(249)=ppB1n
    flx_t(250)=+1.00;flx_SS(250)=0; flx_states(250)=ppR1n
    flx_ostates(250)=ppP1n
    flx_t(251)=+1.00;flx_SS(251)=0; flx_states(251)=ppR1n
    flx_ostates(251)=ppP2n
    flx_t(252)=+1.00;flx_SS(252)=0; flx_states(252)=ppR1n
    flx_ostates(252)=ppP3n
    flx_t(253)=+1.00;flx_SS(253)=0; flx_states(253)=ppR1n
    flx_ostates(253)=ppP4n
    flx_t(254)=+1.00;flx_SS(254)=0; flx_states(254)=ppR1n
    flx_ostates(254)=ppP5n
    flx_t(255)=+1.00;flx_SS(255)=0; flx_states(255)=ppR1n
    flx_ostates(255)=ppP6n

    ! jnetR1p=(R1p<-P.p+B1p)-(R1p->B1p)        (flux perm2):
    flx_calc_nr(31)= 263; flx_CalcIn(31)=iiPel; flx_option(31)=2
    flx_t(256)=-1.00;flx_SS(256)=1; flx_states(256)=ppR1p
    flx_ostates(256)=ppB1p
    flx_t(257)=+1.00;flx_SS(257)=0; flx_states(257)=ppR1p
    flx_ostates(257)=ppB1p
    flx_t(258)=+1.00;flx_SS(258)=0; flx_states(258)=ppR1p
    flx_ostates(258)=ppP1p
    flx_t(259)=+1.00;flx_SS(259)=0; flx_states(259)=ppR1p
    flx_ostates(259)=ppP2p
    flx_t(260)=+1.00;flx_SS(260)=0; flx_states(260)=ppR1p
    flx_ostates(260)=ppP3p
    flx_t(261)=+1.00;flx_SS(261)=0; flx_states(261)=ppR1p
    flx_ostates(261)=ppP4p
    flx_t(262)=+1.00;flx_SS(262)=0; flx_states(262)=ppR1p
    flx_ostates(262)=ppP5p
    flx_t(263)=+1.00;flx_SS(263)=0; flx_states(263)=ppR1p
    flx_ostates(263)=ppP6p

    ! jDIPTn=(P.n<-N.n)-(P.n->N.n)        (flux perm2):
    flx_calc_nr(32)= 287; flx_CalcIn(32)=iiPel; flx_option(32)=2
    flx_t(264)=+1.00;flx_SS(264)=0; flx_states(264)=ppP1n
    flx_ostates(264)=ppN3n
    flx_t(265)=+1.00;flx_SS(265)=0; flx_states(265)=ppP1n
    flx_ostates(265)=ppN4n
    flx_t(266)=-1.00;flx_SS(266)=1; flx_states(266)=ppP1n
    flx_ostates(266)=ppN3n
    flx_t(267)=-1.00;flx_SS(267)=1; flx_states(267)=ppP1n
    flx_ostates(267)=ppN4n
    flx_t(268)=+1.00;flx_SS(268)=0; flx_states(268)=ppP2n
    flx_ostates(268)=ppN3n
    flx_t(269)=+1.00;flx_SS(269)=0; flx_states(269)=ppP2n
    flx_ostates(269)=ppN4n
    flx_t(270)=-1.00;flx_SS(270)=1; flx_states(270)=ppP2n
    flx_ostates(270)=ppN3n
    flx_t(271)=-1.00;flx_SS(271)=1; flx_states(271)=ppP2n
    flx_ostates(271)=ppN4n
    flx_t(272)=+1.00;flx_SS(272)=0; flx_states(272)=ppP3n
    flx_ostates(272)=ppN3n
    flx_t(273)=+1.00;flx_SS(273)=0; flx_states(273)=ppP3n
    flx_ostates(273)=ppN4n
    flx_t(274)=-1.00;flx_SS(274)=1; flx_states(274)=ppP3n
    flx_ostates(274)=ppN3n
    flx_t(275)=-1.00;flx_SS(275)=1; flx_states(275)=ppP3n
    flx_ostates(275)=ppN4n
    flx_t(276)=+1.00;flx_SS(276)=0; flx_states(276)=ppP4n
    flx_ostates(276)=ppN3n
    flx_t(277)=+1.00;flx_SS(277)=0; flx_states(277)=ppP4n
    flx_ostates(277)=ppN4n
    flx_t(278)=-1.00;flx_SS(278)=1; flx_states(278)=ppP4n
    flx_ostates(278)=ppN3n
    flx_t(279)=-1.00;flx_SS(279)=1; flx_states(279)=ppP4n
    flx_ostates(279)=ppN4n
    flx_t(280)=+1.00;flx_SS(280)=0; flx_states(280)=ppP5n
    flx_ostates(280)=ppN3n
    flx_t(281)=+1.00;flx_SS(281)=0; flx_states(281)=ppP5n
    flx_ostates(281)=ppN4n
    flx_t(282)=-1.00;flx_SS(282)=1; flx_states(282)=ppP5n
    flx_ostates(282)=ppN3n
    flx_t(283)=-1.00;flx_SS(283)=1; flx_states(283)=ppP5n
    flx_ostates(283)=ppN4n
    flx_t(284)=+1.00;flx_SS(284)=0; flx_states(284)=ppP6n
    flx_ostates(284)=ppN3n
    flx_t(285)=+1.00;flx_SS(285)=0; flx_states(285)=ppP6n
    flx_ostates(285)=ppN4n
    flx_t(286)=-1.00;flx_SS(286)=1; flx_states(286)=ppP6n
    flx_ostates(286)=ppN3n
    flx_t(287)=-1.00;flx_SS(287)=1; flx_states(287)=ppP6n
    flx_ostates(287)=ppN4n

    ! jDIPTp=(P.p<-N1p)-(P.p->N1p)        (flux perm2):
    flx_calc_nr(33)= 299; flx_CalcIn(33)=iiPel; flx_option(33)=2
    flx_t(288)=+1.00;flx_SS(288)=0; flx_states(288)=ppP1p
    flx_ostates(288)=ppN1p
    flx_t(289)=-1.00;flx_SS(289)=1; flx_states(289)=ppP1p
    flx_ostates(289)=ppN1p
    flx_t(290)=+1.00;flx_SS(290)=0; flx_states(290)=ppP2p
    flx_ostates(290)=ppN1p
    flx_t(291)=-1.00;flx_SS(291)=1; flx_states(291)=ppP2p
    flx_ostates(291)=ppN1p
    flx_t(292)=+1.00;flx_SS(292)=0; flx_states(292)=ppP3p
    flx_ostates(292)=ppN1p
    flx_t(293)=-1.00;flx_SS(293)=1; flx_states(293)=ppP3p
    flx_ostates(293)=ppN1p
    flx_t(294)=+1.00;flx_SS(294)=0; flx_states(294)=ppP4p
    flx_ostates(294)=ppN1p
    flx_t(295)=-1.00;flx_SS(295)=1; flx_states(295)=ppP4p
    flx_ostates(295)=ppN1p
    flx_t(296)=+1.00;flx_SS(296)=0; flx_states(296)=ppP5p
    flx_ostates(296)=ppN1p
    flx_t(297)=-1.00;flx_SS(297)=1; flx_states(297)=ppP5p
    flx_ostates(297)=ppN1p
    flx_t(298)=+1.00;flx_SS(298)=0; flx_states(298)=ppP6p
    flx_ostates(298)=ppN1p
    flx_t(299)=-1.00;flx_SS(299)=1; flx_states(299)=ppP6p
    flx_ostates(299)=ppN1p

    ! jB1DIn=(B1n->N4n)-(B1n<-N.n)        (flux perm2):
    flx_calc_nr(34)= 302; flx_CalcIn(34)=iiPel; flx_option(34)=2
    flx_t(300)=+1.00;flx_SS(300)=1; flx_states(300)=ppB1n
    flx_ostates(300)=ppN4n
    flx_t(301)=-1.00;flx_SS(301)=0; flx_states(301)=ppB1n
    flx_ostates(301)=ppN3n
    flx_t(302)=-1.00;flx_SS(302)=0; flx_states(302)=ppB1n
    flx_ostates(302)=ppN4n

    ! jB1DIp=(B1p->N1p)-(B1p<-N1p)        (flux perm2):
    flx_calc_nr(35)= 304; flx_CalcIn(35)=iiPel; flx_option(35)=2
    flx_t(303)=+1.00;flx_SS(303)=1; flx_states(303)=ppB1p
    flx_ostates(303)=ppN1p
    flx_t(304)=-1.00;flx_SS(304)=0; flx_states(304)=ppB1p
    flx_ostates(304)=ppN1p

    ! jPTMec=P.c->Z4c+Z3c+Z2c        (flux perm2):
    flx_calc_nr(36)= 325; flx_CalcIn(36)=iiPel; flx_option(36)=2
    flx_t(305)=+1.00;flx_SS(305)=1; flx_states(305)=ppP1c
    flx_ostates(305)=ppZ4c
    flx_t(306)=+1.00;flx_SS(306)=1; flx_states(306)=ppP1c
    flx_ostates(306)=ppZ3c
    flx_t(307)=+1.00;flx_SS(307)=1; flx_states(307)=ppP1c
    flx_ostates(307)=ppZ2c
    flx_t(308)=+1.00;flx_SS(308)=1; flx_states(308)=ppP2c
    flx_ostates(308)=ppZ4c
    flx_t(309)=+1.00;flx_SS(309)=1; flx_states(309)=ppP2c
    flx_ostates(309)=ppZ3c
    flx_t(310)=+1.00;flx_SS(310)=1; flx_states(310)=ppP2c
    flx_ostates(310)=ppZ2c
    flx_t(311)=+1.00;flx_SS(311)=1; flx_states(311)=ppP3c
    flx_ostates(311)=ppZ4c
    flx_t(312)=+1.00;flx_SS(312)=1; flx_states(312)=ppP3c
    flx_ostates(312)=ppZ3c
    flx_t(313)=+1.00;flx_SS(313)=1; flx_states(313)=ppP3c
    flx_ostates(313)=ppZ2c
    flx_t(314)=+1.00;flx_SS(314)=1; flx_states(314)=ppP4c
    flx_ostates(314)=ppZ4c
    flx_t(315)=+1.00;flx_SS(315)=1; flx_states(315)=ppP4c
    flx_ostates(315)=ppZ3c
    flx_t(316)=+1.00;flx_SS(316)=1; flx_states(316)=ppP4c
    flx_ostates(316)=ppZ2c
    flx_t(317)=+1.00;flx_SS(317)=1; flx_states(317)=ppP5c
    flx_ostates(317)=ppZ4c
    flx_t(318)=+1.00;flx_SS(318)=1; flx_states(318)=ppP5c
    flx_ostates(318)=ppZ3c
    flx_t(319)=+1.00;flx_SS(319)=1; flx_states(319)=ppP5c
    flx_ostates(319)=ppZ2c
    flx_t(320)=+1.00;flx_SS(320)=1; flx_states(320)=ppP6c
    flx_ostates(320)=ppZ4c
    flx_t(321)=+1.00;flx_SS(321)=1; flx_states(321)=ppP6c
    flx_ostates(321)=ppZ3c
    flx_t(322)=+1.00;flx_SS(322)=1; flx_states(322)=ppP6c
    flx_ostates(322)=ppZ2c
    flx_t(323)=+1.00;flx_SS(323)=1; flx_states(323)=ppPcc
    flx_ostates(323)=ppZ4c
    flx_t(324)=+1.00;flx_SS(324)=1; flx_states(324)=ppPcc
    flx_ostates(324)=ppZ3c
    flx_t(325)=+1.00;flx_SS(325)=1; flx_states(325)=ppPcc
    flx_ostates(325)=ppZ2c

    ! jPTMic=P.c->Z5c+Z6c        (flux perm2):
    flx_calc_nr(37)= 339; flx_CalcIn(37)=iiPel; flx_option(37)=2
    flx_t(326)=+1.00;flx_SS(326)=1; flx_states(326)=ppP1c
    flx_ostates(326)=ppZ5c
    flx_t(327)=+1.00;flx_SS(327)=1; flx_states(327)=ppP1c
    flx_ostates(327)=ppZ6c
    flx_t(328)=+1.00;flx_SS(328)=1; flx_states(328)=ppP2c
    flx_ostates(328)=ppZ5c
    flx_t(329)=+1.00;flx_SS(329)=1; flx_states(329)=ppP2c
    flx_ostates(329)=ppZ6c
    flx_t(330)=+1.00;flx_SS(330)=1; flx_states(330)=ppP3c
    flx_ostates(330)=ppZ5c
    flx_t(331)=+1.00;flx_SS(331)=1; flx_states(331)=ppP3c
    flx_ostates(331)=ppZ6c
    flx_t(332)=+1.00;flx_SS(332)=1; flx_states(332)=ppP4c
    flx_ostates(332)=ppZ5c
    flx_t(333)=+1.00;flx_SS(333)=1; flx_states(333)=ppP4c
    flx_ostates(333)=ppZ6c
    flx_t(334)=+1.00;flx_SS(334)=1; flx_states(334)=ppP5c
    flx_ostates(334)=ppZ5c
    flx_t(335)=+1.00;flx_SS(335)=1; flx_states(335)=ppP5c
    flx_ostates(335)=ppZ6c
    flx_t(336)=+1.00;flx_SS(336)=1; flx_states(336)=ppP6c
    flx_ostates(336)=ppZ5c
    flx_t(337)=+1.00;flx_SS(337)=1; flx_states(337)=ppP6c
    flx_ostates(337)=ppZ6c
    flx_t(338)=+1.00;flx_SS(338)=1; flx_states(338)=ppPcc
    flx_ostates(338)=ppZ5c
    flx_t(339)=+1.00;flx_SS(339)=1; flx_states(339)=ppPcc
    flx_ostates(339)=ppZ6c

    ! jPTRTc=P.c->R1c+R2c+R3c+R6c        (flux perm2):
    flx_calc_nr(38)= 367; flx_CalcIn(38)=iiPel; flx_option(38)=2
    flx_t(340)=+1.00;flx_SS(340)=1; flx_states(340)=ppP1c
    flx_ostates(340)=ppR1c
    flx_t(341)=+1.00;flx_SS(341)=1; flx_states(341)=ppP1c
    flx_ostates(341)=ppR2c
    flx_t(342)=+1.00;flx_SS(342)=1; flx_states(342)=ppP1c
    flx_ostates(342)=ppR3c
    flx_t(343)=+1.00;flx_SS(343)=1; flx_states(343)=ppP1c
    flx_ostates(343)=ppR6c
    flx_t(344)=+1.00;flx_SS(344)=1; flx_states(344)=ppP2c
    flx_ostates(344)=ppR1c
    flx_t(345)=+1.00;flx_SS(345)=1; flx_states(345)=ppP2c
    flx_ostates(345)=ppR2c
    flx_t(346)=+1.00;flx_SS(346)=1; flx_states(346)=ppP2c
    flx_ostates(346)=ppR3c
    flx_t(347)=+1.00;flx_SS(347)=1; flx_states(347)=ppP2c
    flx_ostates(347)=ppR6c
    flx_t(348)=+1.00;flx_SS(348)=1; flx_states(348)=ppP3c
    flx_ostates(348)=ppR1c
    flx_t(349)=+1.00;flx_SS(349)=1; flx_states(349)=ppP3c
    flx_ostates(349)=ppR2c
    flx_t(350)=+1.00;flx_SS(350)=1; flx_states(350)=ppP3c
    flx_ostates(350)=ppR3c
    flx_t(351)=+1.00;flx_SS(351)=1; flx_states(351)=ppP3c
    flx_ostates(351)=ppR6c
    flx_t(352)=+1.00;flx_SS(352)=1; flx_states(352)=ppP4c
    flx_ostates(352)=ppR1c
    flx_t(353)=+1.00;flx_SS(353)=1; flx_states(353)=ppP4c
    flx_ostates(353)=ppR2c
    flx_t(354)=+1.00;flx_SS(354)=1; flx_states(354)=ppP4c
    flx_ostates(354)=ppR3c
    flx_t(355)=+1.00;flx_SS(355)=1; flx_states(355)=ppP4c
    flx_ostates(355)=ppR6c
    flx_t(356)=+1.00;flx_SS(356)=1; flx_states(356)=ppP5c
    flx_ostates(356)=ppR1c
    flx_t(357)=+1.00;flx_SS(357)=1; flx_states(357)=ppP5c
    flx_ostates(357)=ppR2c
    flx_t(358)=+1.00;flx_SS(358)=1; flx_states(358)=ppP5c
    flx_ostates(358)=ppR3c
    flx_t(359)=+1.00;flx_SS(359)=1; flx_states(359)=ppP5c
    flx_ostates(359)=ppR6c
    flx_t(360)=+1.00;flx_SS(360)=1; flx_states(360)=ppP6c
    flx_ostates(360)=ppR1c
    flx_t(361)=+1.00;flx_SS(361)=1; flx_states(361)=ppP6c
    flx_ostates(361)=ppR2c
    flx_t(362)=+1.00;flx_SS(362)=1; flx_states(362)=ppP6c
    flx_ostates(362)=ppR3c
    flx_t(363)=+1.00;flx_SS(363)=1; flx_states(363)=ppP6c
    flx_ostates(363)=ppR6c
    flx_t(364)=+1.00;flx_SS(364)=1; flx_states(364)=ppPcc
    flx_ostates(364)=ppR1c
    flx_t(365)=+1.00;flx_SS(365)=1; flx_states(365)=ppPcc
    flx_ostates(365)=ppR2c
    flx_t(366)=+1.00;flx_SS(366)=1; flx_states(366)=ppPcc
    flx_ostates(366)=ppR3c
    flx_t(367)=+1.00;flx_SS(367)=1; flx_states(367)=ppPcc
    flx_ostates(367)=ppR6c

    ! jPTRTn=P.n->R1n+R6n        (flux perm2):
    flx_calc_nr(39)= 379; flx_CalcIn(39)=iiPel; flx_option(39)=2
    flx_t(368)=+1.00;flx_SS(368)=1; flx_states(368)=ppP1n
    flx_ostates(368)=ppR1n
    flx_t(369)=+1.00;flx_SS(369)=1; flx_states(369)=ppP1n
    flx_ostates(369)=ppR6n
    flx_t(370)=+1.00;flx_SS(370)=1; flx_states(370)=ppP2n
    flx_ostates(370)=ppR1n
    flx_t(371)=+1.00;flx_SS(371)=1; flx_states(371)=ppP2n
    flx_ostates(371)=ppR6n
    flx_t(372)=+1.00;flx_SS(372)=1; flx_states(372)=ppP3n
    flx_ostates(372)=ppR1n
    flx_t(373)=+1.00;flx_SS(373)=1; flx_states(373)=ppP3n
    flx_ostates(373)=ppR6n
    flx_t(374)=+1.00;flx_SS(374)=1; flx_states(374)=ppP4n
    flx_ostates(374)=ppR1n
    flx_t(375)=+1.00;flx_SS(375)=1; flx_states(375)=ppP4n
    flx_ostates(375)=ppR6n
    flx_t(376)=+1.00;flx_SS(376)=1; flx_states(376)=ppP5n
    flx_ostates(376)=ppR1n
    flx_t(377)=+1.00;flx_SS(377)=1; flx_states(377)=ppP5n
    flx_ostates(377)=ppR6n
    flx_t(378)=+1.00;flx_SS(378)=1; flx_states(378)=ppP6n
    flx_ostates(378)=ppR1n
    flx_t(379)=+1.00;flx_SS(379)=1; flx_states(379)=ppP6n
    flx_ostates(379)=ppR6n

    ! jPTRTp=P.n->R1p+R6p        (flux perm2):
    flx_calc_nr(40)= 391; flx_CalcIn(40)=iiPel; flx_option(40)=2
    flx_t(380)=+1.00;flx_SS(380)=1; flx_states(380)=ppP1n
    flx_ostates(380)=ppR1p
    flx_t(381)=+1.00;flx_SS(381)=1; flx_states(381)=ppP1n
    flx_ostates(381)=ppR6p
    flx_t(382)=+1.00;flx_SS(382)=1; flx_states(382)=ppP2n
    flx_ostates(382)=ppR1p
    flx_t(383)=+1.00;flx_SS(383)=1; flx_states(383)=ppP2n
    flx_ostates(383)=ppR6p
    flx_t(384)=+1.00;flx_SS(384)=1; flx_states(384)=ppP3n
    flx_ostates(384)=ppR1p
    flx_t(385)=+1.00;flx_SS(385)=1; flx_states(385)=ppP3n
    flx_ostates(385)=ppR6p
    flx_t(386)=+1.00;flx_SS(386)=1; flx_states(386)=ppP4n
    flx_ostates(386)=ppR1p
    flx_t(387)=+1.00;flx_SS(387)=1; flx_states(387)=ppP4n
    flx_ostates(387)=ppR6p
    flx_t(388)=+1.00;flx_SS(388)=1; flx_states(388)=ppP5n
    flx_ostates(388)=ppR1p
    flx_t(389)=+1.00;flx_SS(389)=1; flx_states(389)=ppP5n
    flx_ostates(389)=ppR6p
    flx_t(390)=+1.00;flx_SS(390)=1; flx_states(390)=ppP6n
    flx_ostates(390)=ppR1p
    flx_t(391)=+1.00;flx_SS(391)=1; flx_states(391)=ppP6n
    flx_ostates(391)=ppR6p

    ! jZ4Z3c=Z3c<-Z2c+Z4c        (flux perm2):
    flx_calc_nr(41)= 393; flx_CalcIn(41)=iiPel; flx_option(41)=2
    flx_t(392)=+1.00;flx_SS(392)=0; flx_states(392)=ppZ3c
    flx_ostates(392)=ppZ2c
    flx_t(393)=+1.00;flx_SS(393)=0; flx_states(393)=ppZ3c
    flx_ostates(393)=ppZ4c

    ! jZIR6c=Z.c->R6c        (flux perm2):
    flx_calc_nr(42)= 398; flx_CalcIn(42)=iiPel; flx_option(42)=2
    flx_t(394)=+1.00;flx_SS(394)=1; flx_states(394)=ppZ3c
    flx_ostates(394)=ppR6c
    flx_t(395)=+1.00;flx_SS(395)=1; flx_states(395)=ppZ4c
    flx_ostates(395)=ppR6c
    flx_t(396)=+1.00;flx_SS(396)=1; flx_states(396)=ppZ2c
    flx_ostates(396)=ppR6c
    flx_t(397)=+1.00;flx_SS(397)=1; flx_states(397)=ppZ5c
    flx_ostates(397)=ppR6c
    flx_t(398)=+1.00;flx_SS(398)=1; flx_states(398)=ppZ6c
    flx_ostates(398)=ppR6c

    ! jZ6Z5c=Z5c<-Z6c        (flux perm2):
    flx_calc_nr(43)= 399; flx_CalcIn(43)=iiPel; flx_option(43)=2
    flx_t(399)=+1.00;flx_SS(399)=0; flx_states(399)=ppZ5c
    flx_ostates(399)=ppZ6c

    ! jPTR2c=R2c<-P.c        (flux perm2):
    flx_calc_nr(44)= 406; flx_CalcIn(44)=iiPel; flx_option(44)=2
    flx_t(400)=+1.00;flx_SS(400)=0; flx_states(400)=ppR2c
    flx_ostates(400)=ppP1c
    flx_t(401)=+1.00;flx_SS(401)=0; flx_states(401)=ppR2c
    flx_ostates(401)=ppP2c
    flx_t(402)=+1.00;flx_SS(402)=0; flx_states(402)=ppR2c
    flx_ostates(402)=ppP3c
    flx_t(403)=+1.00;flx_SS(403)=0; flx_states(403)=ppR2c
    flx_ostates(403)=ppP4c
    flx_t(404)=+1.00;flx_SS(404)=0; flx_states(404)=ppR2c
    flx_ostates(404)=ppP5c
    flx_t(405)=+1.00;flx_SS(405)=0; flx_states(405)=ppR2c
    flx_ostates(405)=ppP6c
    flx_t(406)=+1.00;flx_SS(406)=0; flx_states(406)=ppR2c
    flx_ostates(406)=ppPcc

    ! J1_PTc=P.c        (sedimentation 1):
    flx_calc_nr(45)= 413; flx_CalcIn(45)=iiBen; flx_option(45)=20
    flx_t(407)=+1.00;flx_SS(407)=1; flx_states(407)=ppP1c;flx_ostates(407)=1
    flx_t(408)=+1.00;flx_SS(408)=1; flx_states(408)=ppP2c;flx_ostates(408)=1
    flx_t(409)=+1.00;flx_SS(409)=1; flx_states(409)=ppP3c;flx_ostates(409)=1
    flx_t(410)=+1.00;flx_SS(410)=1; flx_states(410)=ppP4c;flx_ostates(410)=1
    flx_t(411)=+1.00;flx_SS(411)=1; flx_states(411)=ppP5c;flx_ostates(411)=1
    flx_t(412)=+1.00;flx_SS(412)=1; flx_states(412)=ppP6c;flx_ostates(412)=1
    flx_t(413)=+1.00;flx_SS(413)=1; flx_states(413)=ppPcc;flx_ostates(413)=1

    ! J1_R6c=R6c        (sedimentation 1):
    flx_calc_nr(46)= 414; flx_CalcIn(46)=iiBen; flx_option(46)=20
    flx_t(414)=+1.00;flx_SS(414)=1; flx_states(414)=ppR6c;flx_ostates(414)=1

    ! J1_R2c=R2c        (sedimentation 1):
    flx_calc_nr(47)= 415; flx_CalcIn(47)=iiBen; flx_option(47)=20
    flx_t(415)=+1.00;flx_SS(415)=1; flx_states(415)=ppR2c;flx_ostates(415)=1

    ! J0_Ndn=N.n        (sum bot):
    flx_calc_nr(48)= 417; flx_CalcIn(48)=iiBen; flx_option(48)=12
    flx_t(416)=+1.00;flx_SS(416)=1; flx_states(416)=ppN3n
    flx_ostates(416)=ppN3n
    flx_t(417)=+1.00;flx_SS(417)=1; flx_states(417)=ppN4n
    flx_ostates(417)=ppN4n

    ! J0_Ntn=N.n+P.n+R6n        (sum bot):
    flx_calc_nr(49)= 426; flx_CalcIn(49)=iiBen; flx_option(49)=12
    flx_t(418)=+1.00;flx_SS(418)=1; flx_states(418)=ppN3n
    flx_ostates(418)=ppN3n
    flx_t(419)=+1.00;flx_SS(419)=1; flx_states(419)=ppN4n
    flx_ostates(419)=ppN4n
    flx_t(420)=+1.00;flx_SS(420)=1; flx_states(420)=ppP1n
    flx_ostates(420)=ppP1n
    flx_t(421)=+1.00;flx_SS(421)=1; flx_states(421)=ppP2n
    flx_ostates(421)=ppP2n
    flx_t(422)=+1.00;flx_SS(422)=1; flx_states(422)=ppP3n
    flx_ostates(422)=ppP3n
    flx_t(423)=+1.00;flx_SS(423)=1; flx_states(423)=ppP4n
    flx_ostates(423)=ppP4n
    flx_t(424)=+1.00;flx_SS(424)=1; flx_states(424)=ppP5n
    flx_ostates(424)=ppP5n
    flx_t(425)=+1.00;flx_SS(425)=1; flx_states(425)=ppP6n
    flx_ostates(425)=ppP6n
    flx_t(426)=+1.00;flx_SS(426)=1; flx_states(426)=ppR6n
    flx_ostates(426)=ppR6n

    ! qqrivO3c=O3c        (riv):
    flx_calc_nr(50)= 427; flx_CalcIn(50)=iiBen; flx_option(50)=50
    flx_t(427)=+1.00;flx_SS(427)=1; flx_states(427)=ppO3c;flx_ostates(427)=0

    ! qqrivO3h=O3h        (riv):
    flx_calc_nr(51)= 428; flx_CalcIn(51)=iiBen; flx_option(51)=50
    flx_t(428)=+1.00;flx_SS(428)=1; flx_states(428)=ppO3h;flx_ostates(428)=0

    ! qqriv_w= 1        (riv):
    flx_calc_nr(52)= 429; flx_CalcIn(52)=iiBen; flx_option(52)=50
    flx_t(429)=+1.00;flx_SS(429)=1; flx_states(429)=0;flx_ostates(429)=0

    ! qqrivNtn=N3n+N4n+R1n+R6n        (riv):
    flx_calc_nr(53)= 433; flx_CalcIn(53)=iiBen; flx_option(53)=50
    flx_t(430)=+1.00;flx_SS(430)=1; flx_states(430)=ppN3n;flx_ostates(430)=0
    flx_t(431)=+1.00;flx_SS(431)=1; flx_states(431)=ppN4n;flx_ostates(431)=0
    flx_t(432)=+1.00;flx_SS(432)=1; flx_states(432)=ppR1n;flx_ostates(432)=0
    flx_t(433)=+1.00;flx_SS(433)=1; flx_states(433)=ppR6n;flx_ostates(433)=0

    ! qqrivNtp=N1p+R1p+R6p        (riv):
    flx_calc_nr(54)= 436; flx_CalcIn(54)=iiBen; flx_option(54)=50
    flx_t(434)=+1.00;flx_SS(434)=1; flx_states(434)=ppN1p;flx_ostates(434)=0
    flx_t(435)=+1.00;flx_SS(435)=1; flx_states(435)=ppR1p;flx_ostates(435)=0
    flx_t(436)=+1.00;flx_SS(436)=1; flx_states(436)=ppR6p;flx_ostates(436)=0

    ! qqrivR6n=R6n        (riv):
    flx_calc_nr(55)= 437; flx_CalcIn(55)=iiBen; flx_option(55)=50
    flx_t(437)=+1.00;flx_SS(437)=1; flx_states(437)=ppR6n;flx_ostates(437)=0

    ! qqrivR6p=R6p        (riv):
    flx_calc_nr(56)= 438; flx_CalcIn(56)=iiBen; flx_option(56)=50
    flx_t(438)=+1.00;flx_SS(438)=1; flx_states(438)=ppR6p;flx_ostates(438)=0

    ! qqrivR6c=R6c        (riv):
    flx_calc_nr(57)= 439; flx_CalcIn(57)=iiBen; flx_option(57)=50
    flx_t(439)=+1.00;flx_SS(439)=1; flx_states(439)=ppR6c;flx_ostates(439)=0

    ! jBPQ2c=Q2c<-BP.c        (flux):
    flx_calc_nr(58)= 440; flx_CalcIn(58)=iiBen; flx_option(58)=0
    flx_t(440)=+1.00;flx_SS(440)=0; flx_states(440)=ppQ2c
    flx_ostates(440)=ppBP1c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of vars for  track of constituents
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! In this Setup is tracking not active


  end subroutine

