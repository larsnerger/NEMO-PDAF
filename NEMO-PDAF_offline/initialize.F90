!$Id$
!> Initialize model dimensions
!!
!! Routine to perform initialization of the model information for
!! PDAF. Here, the global size of the model domain, the global size
!! of the model state vector and the sizes for decomposition of the 
!! state vector need to be initialized.
!! Generally, this could also be joined with the routine init_pdaf().
!!
!! !REVISION HISTORY:
!! 2022-02 - Lars Nerger - Initial code
!! Later revisions - see svn log
!!
subroutine initialize()
  use mod_kind_pdaf
  use mpi
  use netcdf
  use mod_parallel_pdaf, &
     only: mype=>mype_model, npes=>npes_model, MPIerr, &
     abort_parallel
  use mod_nemo_pdaf, &
       only: jpiglo, jpjglo, jpk, tmask, &
             glamt, gphit, gdept_1d, &            ! reading from domain file
             domain_decomposition, read_decomp, & ! doing domain decomposition
             nimpp, njmpp, nldi, nldj, nlei, nlej, 
             set_nemo_grid, use_wet_state, ndastp
  use mod_statevector_pdaf, only: sfields
  use mod_io_pdaf, &
       only: check, add_slash
  use mod_memcount_pdaf, &
       only: memcount

  implicit none

! *** local variables *** 
  integer :: screen=1         ! Verbosity flag
  integer :: ncid             ! nc file ID
  integer :: varid            ! variable IDs

  integer :: ios              ! ascii file status

  integer :: i, j             ! counters
  logical :: have_pdafnml     ! Flag whether namelist file is present

  integer, allocatable :: k_bot(:, :)       ! bottom level see vertical grid section in nemo manual
  integer, allocatable :: k_top(:, :)       ! top level see vertical grid section in nemo manual

  character(len=50)    :: file_decomp='decomp.txt' ! Name of decomposition file 
  character(len=200)   :: file_ndastp              ! file for getting ndastp, usually a restart file
  ! *** File name and path to read grid information
  character(len=200)   :: path_dims         ! Path for NEMO file holding dimensions
  character(len=80)    :: cn_domcfg         ! File name NEMO file holding dimensions
  namelist /nemo_nml/ screen, path_dims, cn_domcfg, use_wet_state, &
       read_decomp, file_decomp, file_ndastp


! **********************
! *** Initialization ***
! **********************
  path_dims = './'
  cn_domcfg = 'nemo_output.nc'
  file_ndastp = 'restart.nc'

! *** Read namelist file for PDAF if the list is there ***

  inquire( FILE='namelist_cfg.pdaf', EXIST=have_pdafnml ) 
  if (have_pdafnml) then
     open (500,file='namelist_cfg.pdaf')
     read (500,NML=nemo_nml)
     close (500)
  end if
  call add_slash(path_dims)

  ! *************************************
  ! *** Read dimensions of model grid ***
  ! *************************************

  if (mype==0) then
     write (*,'(a,2x,a)') 'NEMO-PDAF', '*** Reading NEMO dimensions ***'
     write (*,'(a,2x,a,a)') 'NEMO-PDAF', 'File: ',trim(path_dims)//trim(cn_domcfg)
  end if

  ! Open the file
  call check( nf90_open(trim(path_dims)//trim(cn_domcfg), nf90_nowrite, ncid) )

  ! Get the dimensions
  call check( nf90_inq_varid(ncid,  'jpiglo', varid) )  
  call check( nf90_get_var(ncid, varid, jpiglo) )
   
  call check( nf90_inq_varid(ncid,  'jpjglo', varid) )  
  call check( nf90_get_var(ncid, varid, jpjglo) )

  call check( nf90_inq_varid(ncid,  'jpkglo', varid) )  
  call check( nf90_get_var(ncid, varid, jpk) )

  ! *************************************************
  ! *** Read coordinates ***
  ! *************************************************
  allocate(glamt(jpiglo, jpjglo), gphit(jpiglo, jpjglo))
  allocate(gdept_1d(jpk))

  call check( nf90_inq_varid(ncid, 'gphit', varid) )
  call check( nf90_get_var(ncid, varid, gphit(:,:), (/1, 1, 1/), (/jpiglo, jpjglo, 1/) ) )
  call check( nf90_inq_varid(ncid, 'glamt', varid) )
  call check( nf90_get_var(ncid, varid, glamt(:,:), (/1, 1, 1/), (/jpiglo, jpjglo, 1/) ) )
  call check( nf90_inq_varid(ncid, 'nav_lev', varid) )
  call check( nf90_get_var(ncid, varid, gdept_1d(:), (/1, 1/), (/jpk, 1/) ) )

  call domain_decomposition(screen, mype, npes, trim(file_decomp))

  ! *************************************************
  ! *** Read field to define mask ***
  ! *************************************************
  allocate(k_bot(nlei, nlej), k_top(nlei, nlej))
  allocate(tmask(nlei, nlej, jpk))

  call check( nf90_inq_varid(ncid, 'bottom_level', varid) )
  call check( nf90_get_var(ncid, varid, k_bot, &
                           (/nimpp+nldi-1, njmpp+nldj-1, 1/), &
                           (/nlei, nlej, 1/) ) )
  call check( nf90_inq_varid(ncid, 'top_level', varid) )
  call check( nf90_get_var(ncid, varid, k_top, &
                           (/nimpp+nldi-1, njmpp+nldj-1, 1/), &
                           (/nlei, nlej, 1/) ) )
  
  ! Close the file. 
  call check( nf90_close(ncid) )

  ! follows dom_msk in dommask.F90 in NEMO
  tmask = 0._pwp
  DO j = 1, nlej
     DO i = 1, nlei
        IF( k_top(i,j) /= 0 ) THEN       ! water in the column
           tmask(i, j, k_top(i,j):k_bot(i,j)) = 1._pwp
        ENDIF
     END DO  
  END DO
  deallocate(k_top, k_bot)
   
  ! get ndastp
  ! Open the file
  call check( nf90_open(trim(file_ndastp), nf90_nowrite, ncid) )

  ! Get the dimensions
  call check( nf90_inq_varid(ncid,  'ndastp', varid) )  
  call check( nf90_get_var(ncid, varid, ndastp) )
  call check( nf90_close(ncid) )

  call MPI_Barrier(MPI_COMM_WORLD, MPIerr)
end subroutine initialize