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
        call imp_spectral_to_UA(NameHemiIn, Potential_II)
      end if

      ! Add bband
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
               engUAwidth(nEngUA), engIMwidth(2,nEngIM)
          integer :: i, j, k, l, hyd_index, ele_index

      character(len=*), parameter:: NameSub = 'imp_spectral_add_potential'
      !------------------------------------------------------------------------
      ! Should replace this later with spline interpolation
      IONO_HYDR_NFlux = 0.
      IONO_ELEC_NFlux = 0.

      engUAwidth(1) = EngUA(2) - EngUA(1)
      engUAwidth(nEngUA) = EngUA(nEngUA) - EngUA(nEngUA - 1)
      engUAwidth(2:nEngUA-1) = (EngUA(3:nEngUA) - EngUA(2:nEngUA-1))/2 &
           + (EngUA(2:nEngUA-1) - EngUA(1:nEngUA-2))/2
      
      do i = 1, 2
         engIMwidth(i,1) = EngIM(i,2) - EngIM(i,1)
         engIMwidth(i,nEngIM) = EngIM(i,nEngIM) - EngIM(i,nEngIM - 1)
         engIMwidth(i,2:nEngIM-1) =	(EngIM(i,3:nEngIM) - EngIM(i,2:nEngIM-1))/2 &
              + (EngIM(i,2:nEngIM-1) - EngIM(i,1:nEngIM-2))/2
         engIMwidth(i,:) = engIMwidth(i,:) * 1000 ! convert to eV
      end do
      !write(*,*) 'engUAwidths: ', engUAwidth
      !write(*,*) 'hydIMwidths: ', engIMwidth(1,:)
      !write(*,*) 'eleIMwidths: ', engIMwidth(2,:)
     
      ! Now do the same for electrons, except the energy grid is different at
      ! every point
      do i = 1, IONO_nTheta; do j = 1, IONO_nPsi
        ! Create new electron energy grid
        new_eGrid_I = EngIM(2, :) * 1000 + PotIn_II(i, j)
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
          if (EngUA(k) < EngIM(1,1)) then 
            hyd_index = -1
            hyd_weight = -1
          ! Greater than largest will be set to 0
          else if (EngUA(k) > EngIM(1,nEngIM)) then 
            hyd_index = -1
            hyd_weight = -1
        ! Interpolate between bounds
          else
            hydIMs: do l = 1, nEngIM - 1
              if (EngUA(k) > EngIM(1,l) .and. EngUA(k) < EngIM(1,l+1)) then
                hyd_index = l
                hyd_weight = (log(EngUA(k)) - log(EngIM(1,l))) / &
                                (log(EngIM(1,l+1)) - log(EngIM(1,l)))
                EXIT hydIMs
              end if
            end do hydIMs
          end if
          if (NameHemiIn == 'north') then
            if(hyd_index >= 0) &
              IONO_HYDR_NFlux(i,j,k) = (1 - hyd_weight) * engUAwidth(k) * &
                                    iono_north_im_nHydrPrec(i,j,hyd_index) &
                                    / engIMwidth(1,hyd_index) &
                                    + hyd_weight * engUAwidth(k) * &
                                    iono_north_im_nHydrPrec(i,j,hyd_index + 1)&
                                    / engIMwidth(1,hyd_index+1)     
            if(ele_index >= 0) &
              IONO_ELEC_NFlux(i,j,k) = (1 - ele_weight) * engUAwidth(k) * &
                                    iono_north_im_nHydrPrec(i,j,ele_index) &
                                    / engIMwidth(2,ele_index) &
                                    + ele_weight * engUAwidth(k) * &
                                    iono_north_im_nHydrPrec(i,j,ele_index + 1)&
                                    / engIMwidth(2,ele_index+1)
          else if (NameHemiIn == 'south') then
            if(hyd_index >= 0) &
              IONO_HYDR_NFlux(i+IONO_nTheta-1,j,k) = (1 - hyd_weight) * engUAwidth(k) *&
                                    iono_south_im_nHydrPrec(i,j,hyd_index) &
                                    / engIMwidth(1,hyd_index) &
                                    + hyd_weight * engUAwidth(k) * &
                                    iono_south_im_nHydrPrec(i,j,hyd_index + 1) &
                                    / engIMwidth(1,hyd_index+1)
            if(ele_index >= 0) &                        
              IONO_ELEC_NFlux(i+IONO_nTheta-1,j,k) = (1 - ele_weight) * engUAwidth(k) *&
                                    iono_south_im_nHydrPrec(i,j,ele_index) &
                                    / engIMwidth(2,ele_index) &
                                    + ele_weight * engUAwidth(k) *&
                                    iono_south_im_nHydrPrec(i,j,ele_index + 1) &
                                    / engIMwidth(2,ele_index+1)
          else
            call CON_stop(NameSub//' : unrecognized hemisphere - '//&
                      NameHemiIn)
          end if
        end do UAs
        if (NameHemiIn == 'north') then
            IONO_HYDR_NFlux(i,j,:) = IONO_HYDR_NFlux(i,j,:) * &
              SUM(IONO_HYDR_NFlux(i,j,:) * engUAwidth) / SUM(iono_north_im_nHydrPrec(i,j,:) * engIMwidth(1,:))
            IONO_ELEC_Nflux(i,j,:) = IONO_ELEC_NFlux(i,j,:) * &
              SUM(IONO_ELEC_NFlux(i,j,:) * engUAwidth) / SUM(iono_north_im_nElecPrec(i,j,:) * engIMwidth(2,:))
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
