!> Utility Routines
!!
!! This module contains several routines useful for common
!! model tasks. The initial routines included output configuration
!! information about the PDAF library, and configuration information
!! about the assimilation parameters.
!! 
module mod_util_pdaf

  use mod_kind_pdaf

  implicit none
  save

contains

  !> This routine performs a model-sided screen output about
  !! the configuration of the data assimilation system.
  !!
  !! **Calling Sequence**
  !!
  !! - Called from: `init_pdaf`
  !!
  subroutine init_info_pdaf()

    use mod_assimilation_pdaf, & ! Variables for assimilation
         only: filtertype, subtype, dim_ens, delt_obs, model_error, &
         model_err_amp, forget, rank_analysis_enkf

    ! *****************************
    ! *** Initial Screen output ***
    ! *****************************

    if (filtertype == 1) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: SEIK'
       if (subtype == 2) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed error-space basis'
       else if (subtype == 3) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed state covariance matrix'
       else if (subtype == 4) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- use ensemble transformation'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    else if (filtertype == 2) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: EnKF'
       if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
       if (rank_analysis_enkf > 0) then
          write (*, '(a,6x, a, i5)') 'NEMO-PDAF', &
               'analysis with pseudo-inverse of HPH, rank:', rank_analysis_enkf
       end if
    else if (filtertype == 3) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: LSEIK'
       if (subtype == 2) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed error-space basis'
       else if (subtype == 3) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed state covariance matrix'
       else if (subtype == 4) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- use ensemble transformation'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    else if (filtertype == 4) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: ETKF'
       if (subtype == 0) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant using T-matrix'
       else if (subtype == 1) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant following Hunt et al. (2007)'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    else if (filtertype == 5) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: LETKF'
       if (subtype == 0) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant using T-matrix'
       else if (subtype == 1) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant following Hunt et al. (2007)'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    else if (filtertype == 6) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: ESTKF'
       if (subtype == 0) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Standard mode'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    else if (filtertype == 7) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: LESTKF'
       if (subtype == 0) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Standard mode'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    end if

  end subroutine init_info_pdaf

!-------------------------------------------------------------------------------

  !> This routine reads the namelist file with parameters
  !! controlling data assimilation with PDAF and outputs to
  !! screen.
  !!
  !! **Calling Sequence**
  !!
  !! - Called from: `init_pdaf`
  !!
  subroutine read_config_pdaf()

    use mod_parallel_pdaf, &
         only: mype_ens
    use mod_assimilation_pdaf, &
         only: filtertype, subtype, type_trans, type_sqrt, &
         locweight  
    use mod_assimilation_pdaf, &
         only: screen, dim_ens, ensscale, delt_obs, &
         type_forget, forget, &
         type_ens_init, type_central_state
    use mod_io_pdaf, &
         only: verbose, path_inistate, path_ens, file_ens, &
         coupling_nemo, save_var_time, save_state, add_slash, &
         file_inistate_date1, file_inistate_date2, &
         file_ens_date1, file_ens_date2, &
         ens_datelist, datestype
    use mod_obs_ssh_mgrid_pdafomi, &
         only: assim_ssh_mgrid, rms_ssh_mgrid, file_ssh_mgrid, &
         lradius_ssh_mgrid, sradius_ssh_mgrid, varname_ssh_mgrid
    use mod_obs_sst_cmems_pdafomi, &
         only: assim_sst_cmems, path_sst_cmems, file_sst_cmems, rms_obs_sst_cmems, &
         lradius_sst_cmems, sradius_sst_cmems, mode_sst_cmems, dist_sst_cmems, &
         varname_sst_cmems

    !< Namelist file
    character(lc) :: nmlfile

    namelist /pdaf_nml/ &
         screen, filtertype, subtype, type_trans, type_sqrt, &
         type_forget, forget, locweight, delt_obs

    namelist /init_nml/ &
         type_ens_init, type_central_state, ensscale, &
         path_inistate, path_ens, file_ens, coupling_nemo, &
         file_inistate_date1, file_inistate_date2, &
         file_ens_date1, file_ens_date2, &
         ens_datelist, datestype, &
         save_var_time, save_state, verbose

    namelist /obs_ssh_mgrid_nml/ &
         assim_ssh_mgrid, rms_ssh_mgrid, file_ssh_mgrid, &
         lradius_ssh_mgrid,  sradius_ssh_mgrid, varname_ssh_mgrid

    namelist /obs_sst_cmems_nml/ &
         assim_sst_cmems, path_sst_cmems, file_sst_cmems, rms_obs_sst_cmems, &
         lradius_sst_cmems,  sradius_sst_cmems, mode_sst_cmems, dist_sst_cmems, &
         varname_sst_cmems


    ! ****************************************************
    ! ***   Initialize PDAF parameters from namelist   ***
    ! ****************************************************

    nmlfile = 'namelist_cfg.pdaf'

    open (20, file=nmlfile)
    read (20, NML=pdaf_nml)
    rewind(20)
    read (20, NML=init_nml)
    rewind(20)
    read (20, NML=obs_ssh_mgrid_nml)
    rewind(20)
    read (20, NML=obs_sst_cmems_nml)
    rewind(20)
    close (20)

    ! *** Add trailing slash to paths ***
    call add_slash(path_sst_cmems)
    call add_slash(path_inistate)
    call add_slash(path_ens)


    ! Print PDAF parameters to screen
    showconf: if (mype_ens == 0) then

       write (*, '(/a,1x,a)') 'NEMO-PDAF','-- Overview of PDAF configuration --'
       write (*, '(a,3x,a)') 'NEMO-PDAF','[pdaf_nml]:'
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','screen       ', screen
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','filtertype   ', filtertype
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','subtype      ', subtype
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_trans   ', type_trans
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_sqrt    ', type_sqrt
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_forget  ', type_forget
       write (*, '(a,5x,a,f10.3)') 'NEMO-PDAF','forget       ', forget
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','locweight    ', locweight
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','delt_obs     ', delt_obs
       write (*, *) ''
       write (*, '(a,3x,a)') 'NEMO-PDAF','[init_nml]:'
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_ens_init      ', type_ens_init
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_central_state ', type_central_state
       write (*, '(a,5x,a,5x,a)') 'NEMO-PDAF','coupling_nemo        ', coupling_nemo
       write (*, '(a,5x,a,5x,f10.2)') 'NEMO-PDAF','ensscale     ', ensscale
       write (*, *) ''
       write (*, '(a,3x,a)') 'NEMO-PDAF','[obs_ssh_mgrid_nml]:'
       write (*, '(a,5x,a,5x,l)') 'NEMO-PDAF','assim_ssh_mgrid      ', assim_ssh_mgrid
       if (assim_ssh_mgrid) then
          write (*, '(a,5x,a,f12.4)') 'NEMO-PDAF','rms_ssh_mgrid      ', rms_ssh_mgrid
          write (*, '(a,5x,a,a)') 'NEMO-PDAF','file_ssh_mgrid     ', trim(file_ssh_mgrid)
          write (*, '(a,5x,a,a)') 'NEMO-PDAF','varname_ssh_mgrid     ', trim(varname_ssh_mgrid)
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','lradius_ssh_mgrid      ', lradius_ssh_mgrid
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','sradius_ssh_mgrid      ', sradius_ssh_mgrid
       end if
       write (*, *) ''
       write (*, '(a,3x,a)') 'NEMO-PDAF','[obs_sst_cmems_nml]:'
       write (*, '(a,5x,a,5x,l)') 'NEMO-PDAF','assim_sst_cmems      ', assim_sst_cmems
       if (assim_sst_cmems) then
          write (*, '(a,5x,a,f12.4)') 'NEMO-PDAF','rms_obs_sst_cmems  ', rms_obs_sst_cmems
          write (*, '(a,5x,a,a)') 'NEMO-PDAF','file_sst_cmems     ', trim(file_sst_cmems)
          write (*, '(a,5x,a,a)') 'NEMO-PDAF','varname_sst_cmems     ', trim(varname_sst_cmems)
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','lradius_sst_cmems      ', lradius_sst_cmems
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','sradius_sst_cmems      ', sradius_sst_cmems
          write (*, '(a,5x,a,i10)') 'NEMO-PDAF','mode_sst_cmems      ', mode_sst_cmems
          write (*, '(a,5x,a,8x,a)') 'NEMO-PDAF','dist_sst_cmems      ', dist_sst_cmems
       end if

       write (*, '(a,1x,a/)') 'NEMO-PDAF','-- End of PDAF configuration overview --'

    end if showconf

  end subroutine read_config_pdaf

!-------------------------------------------------------------------------------

  subroutine finalize_pdaf()

    !>Timing and clean-up of PDAF
    !!
    !! **Calling Sequence**
    !!
    !! - Called from: `nemogcm`
    !!
    !! - Calls: `PDAF_deallocate`
    !!
    use mod_parallel_pdaf, &
         only: mype_ens, comm_ensemble, mpierr
    use mod_iau_pdaf, &
         only: ssh_iau_pdaf, t_iau_pdaf, s_iau_pdaf, u_iau_pdaf, v_iau_pdaf


    ! Show allocated memory for PDAF
    if (mype_ens==0) call PDAF_print_info(2)
    if (mype_ens==0) call PDAF_print_info(11)

    ! Print PDAF timings onto screen
    if (mype_ens==0) call PDAF_print_info(3)

    ! Deallocate PDAF arrays
    call PDAF_deallocate()

    ! Deallocaite IAU arrays
    deallocate (ssh_iau_pdaf)
    deallocate (t_iau_pdaf, s_iau_pdaf, u_iau_pdaf, v_iau_pdaf)

    call mpi_barrier(comm_ensemble, mpierr)

  end subroutine finalize_pdaf

end module mod_util_pdaf
