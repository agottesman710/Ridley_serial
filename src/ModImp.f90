! Copyright (C) 2002 Regents of the University of Michigan,
! portions used with permission
! For more information, see http://csem.engin.umich.edu/tools/swmf

! Inner Magnetosphere Precipitation model

module ModImp

  use ModIonosphere, ONLY: IONO_nTheta, IONO_nPsi

  use ModUtilities, ONLY: CON_stop, CON_set_do_test

  implicit none
  save

  contains
  !============================================================================
  subroutine imp_gen_fluxes(NameHemiIn, AvgEDiffe_II, AvgEDiffi_II, &
          AvgEMono_II, AvgEBbnd_II, EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II,&
          EfluxBbnd_II, LatIn_II)

      use ModIonosphere, ONLY: DoUseIMSpectrum
      use ModMagnit, ONLY: monoenergetic_flux
      use ModIonosphere, ONLY: IONO_NORTH_JR, IONO_SOUTH_JR, &
              IONO_NORTH_invB, IONO_SOUTH_invB, &
              iono_north_im_aveeElec, iono_south_im_aveeElec, &
              iono_north_im_efluxElec, iono_south_im_eFluxElec, &
              iono_north_im_aveeHydr, iono_south_im_aveeHydr, &
              iono_north_im_efluxHydr, iono_south_im_eFluxHydr, &
              iono_north_im_nElecPrec, iono_south_im_nElecPrec             
      use ModConst, ONLY: cKEV

      real, intent(out), dimension(IONO_nTheta, IONO_nPsi) :: &
              AvgEDiffe_II, AvgEDiffi_II, AvgEMono_II, AvgEBbnd_II, &
              EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II, EfluxBbnd_II

      real, intent(in), dimension(IONO_nTheta, IONO_nPsi) :: LatIn_II

      real, dimension(IONO_nTheta, IONO_nPsi) :: &
              FAC_II, OCFL_II, NfluxDiffe_II, ElectronTemp_II, Potential_II

      character(len=*), intent(in) :: NameHemiIn

    character(len=*), parameter:: NameSub = 'imp_gen_fluxes'
    !--------------------------------------------------------------------------
      if (trim(NameHemiIn) == 'south') then
          AvgEDiffe_II = iono_south_im_aveeElec / 1000.0 ! eV to keV
          EfluxDiffe_II = iono_south_im_efluxElec / 1000.0 ! mW/m^2 to W/m^2
          AvgEDiffi_II = iono_south_im_aveeHydr / 1000.0
          EfluxDiffi_II = iono_south_im_efluxHydr / 1000.0
      else if (trim(NameHemiIn) == 'north') then
          AvgEDiffe_II = iono_north_im_aveeElec / 1000.0
          EfluxDiffe_II = iono_north_im_efluxElec / 1000.0
          AvgEDiffi_II = iono_north_im_aveeHydr / 1000.0
          EfluxDiffi_II = iono_north_im_efluxHydr / 1000.0
      else
          call CON_stop(NameSub//' : unrecognized hemisphere - '//&
                  NameHemiIn)
      end if

      ! Limits AvgE (primarily for low latitudes)
      where(AvgEMono_II > 100) AvgEMono_II = 100
      where(AvgEDiffe_II > 100)AvgEDiffe_II = 100
      where(AvgEDiffi_II > 200)AvgEDiffi_II= 200

      if(NameHemiIn == 'north')then
          FAC_II = IONO_NORTH_JR
          OCFL_II = IONO_NORTH_invB
      else if (NameHemiIn == 'south')then
          FAC_II = IONO_SOUTH_JR
          OCFL_II = IONO_SOUTH_invB
       end if

      ! Calculate monoenergetic flux (same as MAGNIT)
      NfluxDiffe_II = EfluxDiffe_II / AvgEDiffe_II / cKEV
      ElectronTemp_II = 2.0 * AvgEDiffe_II / cKEV ! kEV to J(????)
      call monoenergetic_flux(FAC_II, OCFL_II, NfluxDiffe_II, ElectronTemp_II, &
              AvgEDiffe_II, LatIn_II, EfluxMono_II, AvgEMono_II, Potential_II)

      if (DoUseIMSpectrum) then
        call imp_spectral_add_potential(NameHemiIn, Potential_II)
      end if

      ! Add bband
      contains
    !==========================================================================
    subroutine imp_spectral_add_potential(NameHemiIn, PotIn_II)
          use ModIonosphere, ONLY: &
                iono_north_im_nElecPrec, iono_south_im_nElecPrec, &
                nEngIM, EngIM, EngUA

          real, intent(in), dimension(IONO_nTheta, IONO_nPsi) :: PotIn_II
          real, dimension(IONO_nTheta, IONO_nPsi, nEngIM) :: NewNflux_II

          character(len=*), intent(in) :: NameHemiIn
      character(len=*), parameter:: NameSub = 'imp_spectral_add_potential'
      !------------------------------------------------------------------------
      ! need energy grid here
      ! how do we linearly sum energies???.
        ! energy += potential 
      
    end subroutine imp_spectral_add_potential
      !==========================================================================
  end subroutine imp_gen_fluxes
  !============================================================================
end module ModImp
!==============================================================================
