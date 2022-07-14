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
         only: filtertype, subtype, type_trans, &
               type_sqrt, incremental, covartype, &
               rank_analysis_enkf
    use mod_assimilation_pdaf, &
         only: lradius, locweight, sradius
    use mod_assimilation_pdaf, &
         only: screen, dim_ens, ensscale, delt_obs, &
               type_forget, forget, &
               type_ens_init, type_central_state
    use mod_io_pdaf, &
         only: path_inistate, path_ens, coupling_nemo, &
               file_ens, verbose
    ! use mod_obs_ssh_mgrid_pdafomi, &
    !      only: assim_ssh_mgrid, rms_ssh_mgrid, file_ssh_mgrid
    

    !< Namelist file
    character(lc) :: nmlfile

    namelist /pdaf_nml/ screen, dim_ens, ensscale
    namelist /filter_nml/ filtertype, subtype, type_trans, &
                          type_sqrt, incremental, covartype, &
                          rank_analysis_enkf
    namelist /infl_nml/ type_forget, forget
    namelist /local_nml/ lradius, locweight, sradius
    namelist /obs_nml/ delt_obs !, assim_ssh_mgrid, rms_ssh_mgrid, file_ssh_mgrid
    namelist /init_nml/ type_ens_init, type_central_state, &
                        path_inistate, path_ens, coupling_nemo, &
                        file_ens, verbose

    ! ****************************************************
    ! ***   Initialize PDAF parameters from namelist   ***
    ! ****************************************************

    nmlfile = 'namelist_cfg.pdaf'

    open (20, file=nmlfile)
    read (20, NML=pdaf_nml)
    rewind(20)
    read (20, NML=filter_nml)
    rewind(20)
    read (20, NML=infl_nml)
    rewind(20)
    read (20, NML=local_nml)
    rewind(20)
    read (20, NML=obs_nml)
    rewind(20)
    read (20, NML=init_nml)
    rewind(20)
    close (20)

    ! Print PDAF parameters to screen
    showconf: if (mype_ens == 0) then

       write (*, '(/a,1x,a)') 'NEMO-PDAF','-- Overview of PDAF configuration --'
       write (*, '(a,3x,a)') 'NEMO-PDAF','[pdaf_nml]:'
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','screen       ', screen
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','dim_ens      ', dim_ens
       write (*, '(a,3x,a,f10.2)') 'NEMO-PDAF','ensscale     ', ensscale
       write (*, *)
       write (*, '(a,3x,a)') 'NEMO-PDAF','[filter_nml]:'
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','filtertype   ', filtertype
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','subtype      ', subtype
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','type_trans   ', type_trans
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','type_sqrt    ', type_sqrt
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','incremental  ', incremental
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','covartype    ', covartype
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','rank_analysis_enkf ', rank_analysis_enkf
       write (*, *)
       write (*, '(a,3x,a)') 'NEMO-PDAF','[infl_nml]:'
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','type_forget  ', type_forget
       write (*, '(a,3x,a,f10.3)') 'NEMO-PDAF','forget       ', forget
       write (*, *)
       write (*, '(a,3x,a)') 'NEMO-PDAF','[local_nml]:'
       write (*, '(a,3x,a,es10.2)') 'NEMO-PDAF','lradius   ', lradius
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','locweight    ', locweight
       write (*, '(a,3x,a,es10.2)') 'NEMO-PDAF','sradius   ', sradius
       write (*, *)
       write (*, '(a,3x,a)') 'NEMO-PDAF','[obs_nml]:'
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','delt_obs     ', delt_obs
       write (*, *)
       write (*, '(a,3x,a)') 'NEMO-PDAF','[init_nml]:'
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','type_ens_init      ', type_ens_init
       write (*, '(a,3x,a,i10)') 'NEMO-PDAF','type_central_state ', type_central_state
       write (*, '(a,3x,a,a)') 'NEMO-PDAF','coupling_nemo        ', coupling_nemo
       

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
