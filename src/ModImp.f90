! Copyright (C) 2002 Regents of the University of Michigan,
! portions used with permission
! For more information, see http://csem.engin.umich.edu/tools/swmf

module ModImp

  ! Inner Magnetosphere Precipitation model

  use ModIonosphere, ONLY: IONO_nTheta, IONO_nPsi

  use ModUtilities, ONLY: CON_stop, CON_set_do_test

  implicit none
  save

  real, parameter :: ImEfluxFloor = 1E-3, ImAveEFloor = 1E-6

contains
  !============================================================================
  subroutine imp_gen_fluxes(NameHemiIn, AvgEDiffe_II, AvgEDiffi_II, &
       AvgEMono_II, AvgEBbnd_II, EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II,&
       EfluxBbnd_II, LatIn_II)

    use ModMagnit, ONLY: monoenergetic_flux, broadband_flux, smooth_polar_cap,&
                  ConeEfluxDifp, ConeNfluxDifp, ConeEfluxDife, ConeNfluxDife
    use ModIonosphere, ONLY: IONO_NORTH_JR, IONO_SOUTH_JR, &
         IONO_NORTH_invB, IONO_SOUTH_invB, IONO_NORTH_Poynting, &
         IONO_SOUTH_Poynting, IONO_north_im_boundary, &
         IONO_south_im_boundary, DoPolarCapSmoothing, &
         iono_north_im_aveeElec, iono_south_im_aveeElec, &
         iono_north_im_efluxElec, iono_south_im_eFluxElec, &
         iono_north_im_aveeHydr, iono_south_im_aveeHydr, &
         iono_north_im_efluxHydr, iono_south_im_eFluxHydr, &
         DoUseIMSpectrum, iono_north_im_nElecPrec, iono_south_im_nElecPrec, &
         IONO_North_p, IONO_North_rho, IONO_South_p, IONO_South_rho, &
         DoUseGmPe, IONO_North_Pe, IONO_South_Pe
    use ModConst, ONLY: cKEV, cProtonMass, cElectronMass, cPi

    real, intent(out), dimension(IONO_nTheta, IONO_nPsi) :: &
         AvgEDiffe_II, AvgEDiffi_II, AvgEMono_II, AvgEBbnd_II, &
         EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II, EfluxBbnd_II

    real, intent(in), dimension(IONO_nTheta, IONO_nPsi) :: LatIn_II

    real, dimension(IONO_nTheta, IONO_nPsi) :: &
         FAC_II, OCFL_II, NfluxDiffe_II, ElectronTemp_II, Potential_II, &
         Poynting_II, ImBoundary_II, MagP_II, MagNp_II, MagPe_II, MagNe_II, &
         NfluxDiffi_II

    character(len=*), intent(in) :: NameHemiIn

    character(len=*), parameter:: NameSub = 'imp_gen_fluxes'
    !--------------------------------------------------------------------------
    if (trim(NameHemiIn) == 'south') then
      AvgEDiffe_II = iono_south_im_aveeElec ! CIMI gives in keV
      EfluxDiffe_II = iono_south_im_efluxElec / 1000.0 ! mW/m^2 to W/m^2
      AvgEDiffi_II = iono_south_im_aveeHydr
      EfluxDiffi_II = iono_south_im_efluxHydr / 1000.0
    else if (trim(NameHemiIn) == 'north') then
      AvgEDiffe_II = iono_north_im_aveeElec
      EfluxDiffe_II = iono_north_im_efluxElec / 1000.0
      AvgEDiffi_II = iono_north_im_aveeHydr
      EfluxDiffi_II = iono_north_im_efluxHydr / 1000.0
    else
       call CON_stop(NameSub//' : unrecognized hemisphere - '//&
                NameHemiIn)
    end if

    ! Calculate monoenergetic flux (same as MAGNIT)
    NfluxDiffe_II = EfluxDiffe_II / AvgEDiffe_II / cKEV
    where(NfluxDiffe_II < 1E11) AvgEDiffe_II = ImAveEFloor

    ElectronTemp_II = 2.0 * AvgEDiffe_II * cKEV ! kEV to J

    where(Poynting_II < 0) Poynting_II = 0

    if(NameHemiIn == 'north')then
      MagP_II = iono_north_p
      MagNp_II = iono_north_rho / cProtonMass
      FAC_II = IONO_NORTH_JR
      OCFL_II = IONO_NORTH_invB
      Poynting_II = -IONO_NORTH_Poynting
      if (DoUseGMPe) MagPe_II = iono_north_pe
      ImBoundary_II = IONO_north_im_boundary
    else if (NameHemiIn == 'south')then
      MagP_II = iono_south_p
      MagNp_II = iono_south_rho / cProtonMass
      FAC_II = IONO_SOUTH_JR
      OCFL_II = IONO_SOUTH_invB
      Poynting_II = -IONO_SOUTH_Poynting
      if (DoUseGMPe) MagPe_II = iono_south_pe
      ImBoundary_II = IONO_south_im_boundary
    end if

    if(.not. DoUseGMPe) then
      where(OCFL_II < 0)
        MagPe_II = MagP_II
        MagNe_II = MagNe_II
      end where
    else
      where(OCFL_II < 0)
        MagNe_II = MagNe_II
      end where
    end if

    ! Calculate diffuse precipitation: protons.
    AvgEDiffi_II  = MagP_II / MagNp_II  ! Temp = P/nk in Joules
    NfluxDiffi_II = ConeNfluxDifp * MagNp_II * AvgEDiffi_II**0.5 / &
                    sqrt(2 * cPi * cProtonMass)  ! units of #/m2/s
    ! units of W/m2
    EfluxDiffi_II = 2 * NfluxDiffi_II * AvgEDiffi_II * &
            ConeEfluxDifp/ConeNfluxDifp
    ! Recalc to make consistent with ConeFactors (and get units of keV)
    AvgEDiffi_II = EfluxDiffi_II / (NfluxDiffi_II * cKEV)

    where(OCFL_II < 0) 
      ElectronTemp_II  = MagPe_II / MagNe_II  ! T = P/nk in Joules
      NfluxDiffe_II = ConeNfluxDife * MagNe_II * ElectronTemp_II**0.5 / &
                  sqrt(2 * cPi * cElectronMass)  ! units of #/m2/s
      ! units of W/m2
      EfluxDiffe_II = 2 * NfluxDiffe_II * ElectronTemp_II * &
              ConeEfluxDife/ConeNfluxDife
      ! Recalc to make consistent with ConeFactors (and get units of keV)
      AvgEDiffe_II = EfluxDiffe_II / (NfluxDiffe_II * cKEV)
    end where

    ! Smooth area between closed and open field lines
    if (DoPolarCapSmoothing) then
      call smooth_polar_cap(AvgEDiffe_II, OCFL_II, NameHemiIn)
      call smooth_polar_cap(AvgEDiffi_II, OCFL_II, NameHemiIn)
      call smooth_polar_cap(EfluxDiffe_II, OCFL_II, NameHemiIn, ImBoundary_II)
      call smooth_polar_cap(EfluxDiffi_II, OCFL_II, NameHemiIn, ImBoundary_II)
    end if

    call monoenergetic_flux(FAC_II, OCFL_II, NfluxDiffe_II, ElectronTemp_II, &
            AvgEDiffe_II, LatIn_II, EfluxMono_II, AvgEMono_II, &
            PotOut_II=Potential_II, ImIn_II=ImBoundary_II)

    call broadband_flux(Poynting_II, EfluxBbnd_II, AvgEBbnd_II)

    if (DoUseIMSpectrum) then
      call imp_spectral_to_UA(NameHemiIn, Potential_II)
    end if

    contains
    !==========================================================================
      subroutine imp_spectral_to_UA(NameHemiIn, PotIn_II)
          use ModIonosphere, ONLY: &
                iono_north_im_nElecPrec, iono_south_im_nElecPrec, &
                iono_north_im_nHydrPrec, iono_south_im_nHydrPrec, &
                nEngIM, EngIM, nEngUA, EngUA, IONO_HYDR_NFlux, IONO_ELEC_NFlux

          real, intent(in), dimension(IONO_nTheta, IONO_nPsi) :: PotIn_II
          real, dimension(IONO_nTheta, IONO_nPsi, nEngIM) :: NewNflux_II

          character(len=*), intent(in) :: NameHemiIn
          real :: new_eGrid_I(nEngIM), hyd_weight, ele_weight, &
               engUAwidth(nEngUA), engIMwidth(2,nEngIM), engIMeV(2, nEngIM)
          integer :: i, j, k, l, hyd_index, ele_index

      character(len=*), parameter:: NameSub = 'imp_spectral_to_UA'
      !------------------------------------------------------------------------
      ! Should replace this later with spline interpolation
      IONO_HYDR_NFlux = 0.
      IONO_ELEC_NFlux = 0.

      engIMeV = EngIM * 1000

      engUAwidth(1) = EngUA(2) - EngUA(1)
      engUAwidth(nEngUA) = EngUA(nEngUA) - EngUA(nEngUA - 1)
      engUAwidth(2:nEngUA-1) = (EngUA(3:nEngUA) - EngUA(2:nEngUA-1))/2 &
           + (EngUA(2:nEngUA-1) - EngUA(1:nEngUA-2))/2

      do i = 1, 2
         engIMwidth(i,1) = engIMeV(i,2) - engIMeV(i,1)
         engIMwidth(i,nEngIM) = engIMeV(i,nEngIM) - engIMeV(i,nEngIM - 1)
         engIMwidth(i,2:nEngIM-1) =	(engIMeV(i,3:nEngIM) - &
                                    engIMeV(i,2:nEngIM-1))/2&
              + (engIMeV(i,2:nEngIM-1) - engIMeV(i,1:nEngIM-2))/2
         engIMwidth(i,:) = engIMwidth(i,:)
      end do
      ! write(*,*) 'engUAwidths: ', engUAwidth
      ! write(*,*) 'hydIMwidths: ', engIMwidth(1,:)
      ! write(*,*) 'eleIMwidths: ', engIMwidth(2,:)

      ! Now do the same for electrons, except the energy grid is different at
      ! every point
      do i = 1, IONO_nTheta; do j = 1, IONO_nPsi
        ! Create new electron energy grid
        new_eGrid_I = engIMeV(2, :) + PotIn_II(i, j)
        ! Create interp indices to GITM Energy grid
        UAs: do k = 1, nEngUA
          ! Less than smallest will be set to 0
          if (EngUA(k) < new_eGrid_I(1)) then
            ele_index = -1
            ele_weight = -1
          ! Greater than largest will be set to 0
          else if (EngUA(k) > new_eGrid_I(nEngIM)) then
            ele_index = -1
            ele_weight = -1
        ! Interpolate between bounds
          else
            eleIMs: do l = 1, nEngIM - 1
              if (EngUA(k) > new_eGrid_I(l) .and. &
                  EngUA(k) < new_eGrid_I(l+1)) then
                ele_index = l
                ele_weight = (log(EngUA(k)) - log(new_eGrid_I(l))) / &
                                (log(new_eGrid_I(l+1)) - log(new_eGrid_I(l)))
                EXIT eleIMs
              end if
            end do eleIMs
          end if
          ! Do for ions
          if (EngUA(k) < engIMeV(1,1)) then
            hyd_index = -1
            hyd_weight = -1
          ! Greater than largest will be set to 0
          else if (EngUA(k) > engIMeV(1,nEngIM)) then
            hyd_index = -1
            hyd_weight = -1
        ! Interpolate between bounds
          else
            hydIMs: do l = 1, nEngIM - 1
              if (EngUA(k) > engIMeV(1,l) .and. EngUA(k) < engIMeV(1,l+1)) then
                hyd_index = l
                hyd_weight = (log(EngUA(k)) - log(engIMeV(1,l))) / &
                                (log(engIMeV(1,l+1)) - log(engIMeV(1,l)))
                EXIT hydIMs
              end if
            end do hydIMs
          end if
          if (NameHemiIn == 'north') then
            if(hyd_index >= 0) &
              IONO_HYDR_NFlux(i,j,k) = (1 - hyd_weight) * &
                                    iono_north_im_nHydrPrec(i,j,hyd_index) &
                                    + hyd_weight * &
                                    iono_north_im_nHydrPrec(i,j,hyd_index + 1)
            if(ele_index >= 0) &
              IONO_ELEC_NFlux(i,j,k) = (1 - ele_weight) * &
                                    iono_north_im_nElecPrec(i,j,ele_index) &
                                    + ele_weight * &
                                    iono_north_im_nElecPrec(i,j,ele_index + 1)
          else if (NameHemiIn == 'south') then
            if(hyd_index >= 0) &
              IONO_HYDR_NFlux(i+IONO_nTheta-1,j,k) = (1 - hyd_weight) * &
                                    iono_south_im_nHydrPrec(i,j,hyd_index) &
                                    + hyd_weight * &
                                    iono_south_im_nHydrPrec(i,j,hyd_index + 1)
            if(ele_index >= 0) &
              IONO_ELEC_NFlux(i+IONO_nTheta-1,j,k) = (1 - ele_weight) * &
                                    iono_south_im_nElecPrec(i,j,ele_index) &
                                    + ele_weight * &
                                    iono_south_im_nElecPrec(i,j,ele_index + 1)
          else
            call CON_stop(NameSub//' : unrecognized hemisphere - '//&
                      NameHemiIn)
          end if
        end do UAs
        if (NameHemiIn == 'north') then
            IONO_HYDR_NFlux(i,j,:) = IONO_HYDR_NFlux(i,j,:) * &
              SUM(IONO_HYDR_NFlux(i,j,:) * engUAwidth) / &
              SUM(iono_north_im_nHydrPrec(i,j,:) * engIMwidth(1,:))
            IONO_ELEC_Nflux(i,j,:) = IONO_ELEC_NFlux(i,j,:) * &
              SUM(IONO_ELEC_NFlux(i,j,:) * engUAwidth) / &
              SUM(iono_north_im_nElecPrec(i,j,:) * engIMwidth(2,:))
        else if (NameHemiIn == 'south') then
            IONO_HYDR_NFlux(i+IONO_nTheta-1,j,:) = &
              IONO_HYDR_NFlux(i+IONO_nTheta-1,j,:) &
              * SUM(IONO_HYDR_NFlux(i+IONO_nTheta-1,j,:) * engUAwidth) &
              / SUM(iono_south_im_nHydrPrec(i,j,:) * engIMwidth(1,:))
            IONO_ELEC_Nflux(i+IONO_nTheta-1,j,:) = &
              IONO_ELEC_NFlux(i+IONO_nTheta-1,j,:) &
              * SUM(IONO_ELEC_NFlux(i+IONO_nTheta-1,j,:) * engUAwidth) &
              / SUM(iono_south_im_nElecPrec(i,j,:) * engIMwidth(2,:))
        else
          call CON_stop(NameSub//' : unrecognized hemisphere - '//&
                      NameHemiIn)
        endif
      enddo; enddo

    end subroutine imp_spectral_to_UA
    !==========================================================================
  end subroutine imp_gen_fluxes
  !============================================================================
end module ModImp
!==============================================================================
