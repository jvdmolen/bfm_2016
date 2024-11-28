!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: InitTransportStateTypes
!
! DESCRIPTION
!   Defining way of transport/integration of for Statevariables

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine InitTransportStateTypes
!
! USES:
  use global_mem
  use mem_param,only:CalcPhytoPlankton,CalcMesoZooPlankton,CalcMicroZooPlankton
  use mem_PelBac,only:p_version_PelBac =>p_version
  use mem_phyto,only:p_qnR2c
  use mem

!  
!
! !AUTHORS
!   mfstep/ERSEM team
!
! !REVISION_HISTORY
!   ---
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
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
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Setting of type for transport/integration  Pelagic state variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  implicit none
  integer               ::i,j,k
  
  !Reset array D3STATETYPE
  D3STATETYPE(:)=NOTRANSPORT
  !Set elements
  D3STATETYPE(ppR9x)=SILTTRANSPORT
  D3STATETYPE(ppO2o)=ALLTRANSPORT
  D3STATETYPE(ppN1p)=ALLTRANSPORT
  D3STATETYPE(ppN3n)=ALLTRANSPORT
  D3STATETYPE(ppN4n)=ALLTRANSPORT
  D3STATETYPE(ppN5s)=ALLTRANSPORT
  D3STATETYPE(ppN6r)=ALLTRANSPORT
  D3STATETYPE(ppB1c)=ALLTRANSPORT
  D3STATETYPE(ppB1n)=ALLTRANSPORT+1
  D3STATETYPE(ppB1p)=ALLTRANSPORT+1
  if (p_version_PelBac== 4 .or. p_version_PelBac== 5) &
      D3STATETYPE(ppBac)=ALLTRANSPORT
  write(LOGUNIT,*) 'ppBac Transported?',p_version_PelBac,ppBac,D3STATETYPE(ppBac)

  do i=1,iiPhytoPlankton
    if (CalcPhytoPlankton(i)) then
      k=ALLTRANSPORT;!if ( i==iiP5) k=SILTTRANSPORT
      D3STATETYPE(ppPhytoPlankton(i,iiC))=k
      D3STATETYPE(ppPhytoPlankton(i,iiN))=k+1
      D3STATETYPE(ppPhytoPlankton(i,iiP))=k+1
      D3STATETYPE(ppPhytoPlankton(i,iiL))=k+1
      j=ppPhytoPlankton(i,iiS);if (j.gt.0) D3STATETYPE(j)=k+1
    endif
  enddo
  if (CalcPhytoPlankton(iiP6))D3STATETYPE(ppPcc)=ALLTRANSPORT+1
  do i=1,iiMesoZooPlankton
    if (CalcMesoZooPlankton(i)) then
        D3STATETYPE(ppMesoZooPlankton(i,iiC))=ALLTRANSPORT
        j=ppMesoZooPlankton(i,iiN);if(j.gt.0) D3STATETYPE(j)=ALLTRANSPORT+1
        j=ppMesoZooPlankton(i,iiP);if(j.gt.0) D3STATETYPE(j)=ALLTRANSPORT+1
     endif
  enddo
  do i=1,iiMicroZooPlankton
    if (CalcMesoZooPlankton(i)) then
        D3STATETYPE(ppMicroZooPlankton(i,iiC))=ALLTRANSPORT
        j=ppMicroZooPlankton(i,iiN);if(j.gt. 0) D3STATETYPE(j)=ALLTRANSPORT+1
        j=ppMicroZooPlankton(i,iiP);if(j.gt. 0) D3STATETYPE(j)=ALLTRANSPORT+1
     endif
  enddo
  D3STATETYPE(ppR1c)=ALLTRANSPORT
  D3STATETYPE(ppR1n)=ALLTRANSPORT
  D3STATETYPE(ppR1p)=ALLTRANSPORT
  D3STATETYPE(ppR2c)=ALLTRANSPORT
  if (p_qnR2c>NZERO) D3STATETYPE(ppR2n)=ALLTRANSPORT
  D3STATETYPE(ppR3c)=ALLTRANSPORT
  D3STATETYPE(ppRZc)=ALLTRANSPORT
  D3STATETYPE(ppR6c)=ALLTRANSPORT
  D3STATETYPE(ppR6n)=ALLTRANSPORT
  D3STATETYPE(ppR6p)=ALLTRANSPORT
  D3STATETYPE(ppR6s)=ALLTRANSPORT
! D3STATETYPE(ppR7c)=NOTRANSPORT
  if ( ppO3c > 0 ) D3STATETYPE(ppO3c)=ALLTRANSPORT
  if ( ppO3h > 0 ) D3STATETYPE(ppO3h)=ALLTRANSPORT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Setting of type for transport/integration  Benthic state variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  D2STATETYPE(ppY1c)=NOTRANSPORT
  D2STATETYPE(ppY2c)=NOTRANSPORT
  D2STATETYPE(ppY3c)=NOTRANSPORT
  D2STATETYPE(ppYy3c)=NOTRANSPORT
  D2STATETYPE(ppY4c)=NOTRANSPORT
  D2STATETYPE(ppY5c)=NOTRANSPORT
  D2STATETYPE(ppY1n)=NOTRANSPORT
  D2STATETYPE(ppY2n)=NOTRANSPORT
  D2STATETYPE(ppY3n)=NOTRANSPORT
  D2STATETYPE(ppY4n)=NOTRANSPORT
  D2STATETYPE(ppY5n)=NOTRANSPORT
  D2STATETYPE(ppY1p)=NOTRANSPORT
  D2STATETYPE(ppY2p)=NOTRANSPORT
  D2STATETYPE(ppY3p)=NOTRANSPORT
  D2STATETYPE(ppY4p)=NOTRANSPORT
  D2STATETYPE(ppY5p)=NOTRANSPORT
  D2STATETYPE(ppQ6c)=NOTRANSPORT
  D2STATETYPE(ppQ6n)=NOTRANSPORT
  D2STATETYPE(ppQ6p)=NOTRANSPORT
  D2STATETYPE(ppQ6s)=NOTRANSPORT
  D2STATETYPE(ppQ1c)=NOTRANSPORT
  D2STATETYPE(ppQ11c)=NOTRANSPORT
  D2STATETYPE(ppQ1n)=NOTRANSPORT
  D2STATETYPE(ppQ11n)=NOTRANSPORT
  D2STATETYPE(ppQ1p)=NOTRANSPORT
  D2STATETYPE(ppQ11p)=NOTRANSPORT
  D2STATETYPE(ppH1c)=NOTRANSPORT
  D2STATETYPE(ppH2c)=NOTRANSPORT
  D2STATETYPE(ppH3c)=NOTRANSPORT
  D2STATETYPE(ppH1n)=NOTRANSPORT
  D2STATETYPE(ppH2n)=NOTRANSPORT
  D2STATETYPE(ppH3n)=NOTRANSPORT
  D2STATETYPE(ppH1p)=NOTRANSPORT
  D2STATETYPE(ppH2p)=NOTRANSPORT
  D2STATETYPE(ppH3p)=NOTRANSPORT
  D2STATETYPE(ppK1p)=NOTRANSPORT
  D2STATETYPE(ppK11p)=NOTRANSPORT
  D2STATETYPE(ppK21p)=NOTRANSPORT
  D2STATETYPE(ppK4n)=NOTRANSPORT
  D2STATETYPE(ppK14n)=NOTRANSPORT
  D2STATETYPE(ppK24n)=NOTRANSPORT
  D2STATETYPE(ppK3n)=NOTRANSPORT
  D2STATETYPE(ppK5s)=NOTRANSPORT
  D2STATETYPE(ppK6r)=NOTRANSPORT
  D2STATETYPE(ppG2o)=NOTRANSPORT
  D2STATETYPE(ppG4n)=NOTRANSPORT
  D2STATETYPE(ppD1m)=NOTRANSPORT
  D2STATETYPE(ppD2m)=NOTRANSPORT
  D2STATETYPE(ppD6m)=NOTRANSPORT
  D2STATETYPE(ppD7m)=NOTRANSPORT
  D2STATETYPE(ppD8m)=NOTRANSPORT
  D2STATETYPE(ppD9m)=NOTRANSPORT
  if ( ppG3c > 0 ) D2STATETYPE(ppG3c)=NOTRANSPORT
  if ( ppG3h > 0 ) D2STATETYPE(ppG3h)=NOTRANSPORT
  if ( ppG13c > 0 ) D2STATETYPE(ppG13c)=NOTRANSPORT
  if ( ppG13h > 0 ) D2STATETYPE(ppG13h)=NOTRANSPORT
  end subroutine
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
