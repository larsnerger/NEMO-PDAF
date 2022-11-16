#!/bin/bash
#SBATCH --job-name=EGMDA
#SBATCH -p mpp
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=96
#SBATCH --time=1:00:00
##SBATCH --partition=standard96:test
#SBATCH --partition=standard96
##SBATCH --mail-user=
#SBATCH --mail-type=END
##SBATCH -A 

# Run script for HLRN-EMMY using Intel compiler version 2022

# To prepare only the run directories without running
# set $prepare=1, $dorun=0, $postproc=0 
# and run interactively in the shell

module load intel/2022.2
module load openmpi/intel/4.1.4
#module load hdf5-parallel/ompi/intel/1.12.1
#module load netcdf-parallel/ompi/intel.22/4.9.1

ulimit -s unlimited

export LD_LIBRARY_PATH=/sw/dataformats/netcdf-parallel/ompi/intel.22/4.9.1/skl/lib:/sw/dataformats/hdf5-parallel/ompi/intel.22/1.12.1/skl/lib:/sw/tools/oneapi/2022.2/compiler/2022.1.0/linux/compiler/lib/intel64_lin:/sw/tools/oneapi/2022.2/mkl/2022.1.0/lib/intel64

# ---------------------------------------------------------------------------------------------------

yy=2015 #Start year
mm=01   #Start month
dd=01   #Start day

tstr2='20150101' #End date yyyymmdd

np_nemo=186  #number of PE's for nemo
np_xios=6    #number of PE's for xios
np=$np_nemo

NENS=4  # Ensemble size

# Whether we initialize with distributed restart files: 1=true, 0=global restart files
restart_dis=1

# Whether the script prepares the run directories, runs the experiment, does posptprocessing
prepare=1
dorun=0
postproc=0

# ---------------------------------------------------------------------------------------------------

# Set initial date  $yy$mm$dd==$initial_date we use the restart files from the input directory

initial_date=20150101

# ---------------------------------------------------------------------------------------------------

#set -xv #debugging

# ---------------------------------------------------------------------------------------------------
#nemo_exe_dir='/home/hbkycsun/nemo4_dev-nemo4_ergom_allEuler/cfgs/NEMO_ERGOM_PDAF/BLD/bin'
#nemo_exe_dir='/home/hzfblner/SEAMLESS/nemo4_dev_allEuler/cfgs/NEMO-ERGOM-PDAF/BLD/bin'
nemo_exe_dir='/home/hzfblner/SEAMLESS/nemo4_dev_allEuler/cfgs/NEMO-ERGOM-PDAF_intel22/BLD/bin'
xios_exe_dir='/home/hzfblner/SEAMLESS/xios-2.0_par_intel22/bin'
#rebuild='/home_ad/bm1405/balmfc_git/nemo4_dev/tools/REBUILD_NEMO/'
restart_out='output/restarts'
archive='/scratch/usr/hzfblner/SEAMLESS/forcing_ergom_allEuler'
archive_ln='/scratch/usr/hzfblner/SEAMLESS/forcings'
inputs_nc='/scratch/usr/hzfblner/SEAMLESS/run/inputs_allEuler_nc4'
#setup_store='/scratch/usr/hbkycsun/data/setup_store'
setup_store='/scratch/usr/hzfblner/SEAMLESS/run/config_ERGOM_allEuler'
initialdir='/scratch/usr/hzfblner/SEAMLESS/restart'
disrestartdir='/scratch/usr/hzfblner/SEAMLESS/run/restart_dist_20150101'
# ---------------------------------------------------------------------------------------------------

export nemo_exe_dir
export restart_out
export archive
#export rebuild

# ---------------------------------------------------------------------------------------------------

# Prepare run directories and files
if [ $prepare -eq 1 ]; then

    # Prepare PDAF namelist
    cp $setup_store/namelist_cfg.pdaf_template ./
    cat namelist_cfg.pdaf_template     \
	| sed -e "s:_DIMENS_:$NENS:"     \
	> namelist_cfg.pdaf

    echo 'Prepare XML files ...'
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      ENSstr2=`printf %02d $i`
     #echo 'preparing context_nemo.xml...'
      cat $setup_store/context_nemoXXX.xml_template   \
	  | sed -e "s:_DIMENS3_:${ENSstr}:g"   \
	  > context_nemo${ENSstr}.xml
     #echo 'preparing file_def_nemo-ice.xml...'
      cat $setup_store/file_def_nemo-iceXXX.xml_template   \
	  | sed -e "{s:_DIMENS3_:${ENSstr}:g;s:_DIMENS2_:${ENSstr2}:g}"   \
	  > file_def_nemo-ice${ENSstr}.xml
     #echo 'preparing file_def_nemo-oce.xml...'
      cat $setup_store/file_def_nemo-oceXXX.xml_template   \
	  | sed -e "{s:_DIMENS3_:${ENSstr}:g;s:_DIMENS2_:${ENSstr2}:g}"   \
	  > file_def_nemo-oce${ENSstr}.xml
     #echo 'preparing file_def_nemo-ergom.xml...'
      cat $setup_store/file_def_nemo-ergomXXX.xml_template   \
	  | sed -e "{s:_DIMENS3_:${ENSstr}:g;s:_DIMENS2_:${ENSstr2}:g}"   \
	  > file_def_nemo-ergom${ENSstr}.xml
    done

    cp  $setup_store/iodef.xml_template iodef.xml
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      echo "  <context id=\"nemo_${ENSstr}\" src=\"./context_nemo${ENSstr}.xml\"/> " >> iodef.xml
    done
    echo "</simulation>" >> iodef.xml

    # set working directories for different ensemble members
    echo ' '
    echo 'Creating ensemble working directories...'
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      if [ ! -d ${ENSstr} ]; then
	  mkdir -p ${ENSstr}
      else
	  rm -rf ${ENSstr}    # delete output from previous test runs
      fi
      wdir=`pwd`/${ENSstr}
      export wdir
      mkdir -p $wdir/output/restarts
      mkdir -p $wdir/output/log
      mkdir -p $wdir/output/data
      mkdir -p $wdir/output/DA
      mkdir -p $wdir/initialstate
      mkdir -p $wdir/forcing

      echo 'Run directory: ' $wdir
    done
    echo ' '

# ---------------------------------------------------------------------------------------------------
    echo 'linking executables...'
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}
      export wdir
      if [ ! -f ${ENSstr}/xios_server.exe ]; then
	  ln -s $xios_exe_dir/xios_server.exe $wdir/
      fi
      if [ ! -f ${ENSstr}/nemo.exe ]; then
	  ln -s $nemo_exe_dir/nemo.exe $wdir/
      fi
    done
# ---------------------------------------------------------------------------------------------------
    tstr0=`date +%Y%m%d`

    sdte=`date --date "$dte0 0 day" +'%Y-%m-%d %H:%M:%S'`
    yyp1=`date +'%Y' -d"$sdte"`                     #  formating
    mmp1=`date +'%m' -d"$sdte"`
    ddp1=`date +'%d' -d"$sdte"`
    tstr="$yy$mm$dd"
    tstr_ini=$tstr    # Store date at start of run

    echo 'Simulation with nemo executable from ' $nemo_exe_dir

# ---------------------------------------------------------------------------------------------------

    #linking restart

    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}
      export wdir
      echo 'linking restart...'
      if [ $tstr -eq $initial_date ]; then
	  if [ $restart_dis -eq 0 ]; then
	      echo "Link initial restart files for full domain"
	      if [ ! -f $wdir/initialstate/restart_in.nc ]; then
		  ln -s $initialdir/NORDIC_2015010100_restart.nc $wdir/initialstate/restart_in.nc
	      fi
	      if [ ! -f $wdir/initialstate/restart_ice_in.nc ]; then
		  ln -s $initialdir/NORDIC_2015010100_restart_ice.nc $wdir/initialstate/restart_ice_in.nc
	      fi
	      if [ ! -f $wdir/initialstate/restart_trc_in.nc ]; then
		  ln -s $initialdir/NORDIC_2015010100_restart_trc.nc $wdir/initialstate/restart_trc_in.nc
	      fi
	      if [ ! -f $wdir/initialstate/NORDIC-NS1_restart_ptrc_in_iowfix.nc ]; then
		  ln -s $inputs_nc/initialstate/NORDIC-NS1_restart_ptrc_in_iowfix.nc $wdir/initialstate/NORDIC-NS1_restart_ptrc_in_iowfix.nc
		  ln -s $inputs_nc/initialstate/NORDIC-NS1_restart_ptrc_in_iowfix.nc $wdir/NORDIC-NS1_restart_ptrc_in_iowfix.nc
	      fi
	      if [ ! -f $wdir/initialstate/NORDIC-NS1_restart_ptrc_in.nc ]; then
		  ln -s $inputs_nc/initialstate/NORDIC-NS1_restart_ptrc_in.nc $wdir/initialstate/NORDIC-NS1_restart_ptrc_in.nc
		  ln -s $inputs_nc/initialstate/NORDIC-NS1_restart_ptrc_in.nc $wdir/NORDIC-NS1_restart_ptrc_in.nc
	      fi
	  else
	      echo "Link distributed initial restart files"
	      for((i=1;i<=$NENS;i++))
		do
		ENSstr=`printf %03d $i`
		wdir=`pwd`/${ENSstr}
		export wdir
		ln -s ${disrestartdir}/* ${wdir}/initialstate
	      done
	  fi
      else
	  echo "Use distributed restart files from previous run"
      fi
    done

# ---------------------------------------------------------------------------------------------------
    echo '----------------------------------------'
    echo 'Preparing forcing'
    echo 'Time period from ' $tstr 'until' $tstr2

    echo 'Linking files with fixed values ...'
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}
      export wdir
    #echo 'linking input...'
      if [ ! -f $wdir/forcing/river_data.nc ]; then
	  ln -s $inputs_nc/river_data.nc $wdir/forcing
      fi
      if [ ! -f $wdir/forcing/dum12_y2015.nc ]; then
	  ln -s $inputs_nc/dum12_y2015.nc $wdir/forcing
      fi
      if [ ! -f $wdir/forcing/weights_bilin.nc ]; then
	  ln -s $inputs_nc/weights_bilin.nc $wdir/forcing
	  ln -s $inputs_nc/weights_bilin.nc $wdir/
      fi
      if [ ! -f $wdir/forcing/weights_bicubic.nc ]; then
	  ln -s $inputs_nc/weights_bicubic.nc $wdir/forcing
	  ln -s $inputs_nc/weights_bicubic.nc $wdir/
      fi

      ln -s $inputs_nc/bdytide*.nc            $wdir/
      ln -s $inputs_nc/coordinates.bdy.nc     $wdir/
      ln -s $inputs_nc/chlorophyll.nc         $wdir/
      ln -s $inputs_nc/benthos_null.nc        $wdir/
      ln -s $inputs_nc/domain_cfg.nc          $wdir/
      ln -s $inputs_nc/eddy_viscosity_3D.nc   $wdir/
      ln -s $inputs_nc/eddy_diffusivity_3D.nc $wdir/
    #ln -s $inputs_nc/np_ergom.nc            $wdir/

    # for ERGOM_allEuler
      ln -s $inputs_nc/bfr_roughness.nc       $wdir/
      ln -s $inputs_nc/carbon.nc              $wdir/
      ln -s $inputs_nc/iron_dummy.nc          $wdir/
      ln -s $inputs_nc/np_ergom_c16.nc        $wdir/
      ln -s $inputs_nc/sed_init_1k.nc         $wdir/
      ln -s $inputs_nc/z2d_ben201401.nc       $wdir/

    # for PDAF
      ln -s $setup_store/*.txt       $wdir/

    done

    #Loop over dates

    ndays=0
    while [ "$(date -d "$tstr" +%Y%m%d)" -le "$(date -d "$tstr2" +%Y%m%d)" ]; do
	date_nemo=$(date -d "$tstr" +y%Ym%md%d)
	date_force=$(date -d "$tstr" +%Y%m%d)
	date_NS01=$(date -d "$tstr" +y%Ym%m)
	tstr=$(date -d "$tstr" +%Y%m%d)

	echo 'Link time-varying forcing files ... date ' $tstr
	for((i=1;i<=$NENS;i++))
	  do
	  ENSstr=`printf %03d $i`
	  wdir=`pwd`/${ENSstr}
	  export wdir
#    echo 'Preparing input... member ' $i

          #1. get river forcing
	  if [ -f $wdir/forcing/EHYPE_$date_nemo.nc ]; then
	      echo 'Use file ' $wdir/forcing/EHYPE_$date_nemo'+024H.nc'
	  else
	      ln -s $archive_ln/EHYPE_$date_nemo'+024H.nc' $wdir/forcing/EHYPE_$date_nemo.nc || { echo '1. failed' ; exit 1; }
#      echo '************************************'
#      echo '1. river forcing done'
	  fi

          #2. get atm. force
	  if [ -f $wdir/forcing/FORCE_$date_nemo.nc ]; then
	      echo 'Use file ' $wdir/forcing/FORCE_$tstr'+24.nc'
	  else
	      ln -s $archive_ln/FORCE_$tstr'+24.nc' $wdir/forcing/FORCE_$date_nemo.nc  #|| { echo '2. failed' ; exit 1; }
#      echo '************************************'
#      echo '2. atm. forcing done'
	  fi

          #3a. get uhv bdy
	  if [ -f $wdir/forcing/bdy_uvh_$date_nemo.nc ]; then
	      echo 'Use file ' $wdir/forcing/bdy_uvh_$tstr.nc
	  else
	      ln -s $archive_ln/bdy_uvh_$tstr.nc  $wdir/forcing/bdy_uvh_$date_nemo.nc || { echo '3. failed' ; exit 1; }
#      echo '************************************'
#      echo '3. bdy done'
	  fi

          #4. get ts bdy
	  if [ -f $wdir/forcing/bdy_ts_$date_nemo.nc ]; then
	      echo 'Use file ' $wdir/forcing/bdy_ts_$tstr.nc
	  else
	      ln -s $archive_ln/bdy_ts_$tstr.nc  $wdir/forcing/bdy_ts_$date_nemo.nc || { echo '4. failed' ; exit 1; }
#      echo '************************************'
#      echo '4. obc done'
	  fi

          #5. ERGOM river loads
	  if [ -f $wdir/forcing/ERGOM_CBC_$date_nemo.nc ]; then
	      echo 'Use file ' $wdir/forcing/ERGOM_CBC_$date_nemo.nc
	  else
	      ln -s $archive/ERGOM_CBC_$date_nemo.nc  $wdir/forcing/ERGOM_CBC_$date_nemo.nc || { echo '5. failed' ; exit 1; }
	    #$wdir/do_ergom_cbc-rivers $date_nemo || { echo '5. failed' ; exit 1; }
#      echo '************************************'
#      echo '5. ergom river loads done'
	  fi

          #6. ergom surface bdy conditions
	  if [ -f $wdir/forcing/ERGOM_SBC_$date_nemo.nc ]; then
	      echo 'Use file ' $wdir/forcing/ERGOM_SBC_$date_nemo.nc
	  else
	      ln -s $archive/ERGOM_SBC_$date_nemo.nc  $wdir/forcing/ERGOM_SBC_$date_nemo.nc || { echo '6. failed' ; exit 1; }
      #$wdir/do_ergom_sbc-depo $date_nemo || { echo '5. failed' ; exit 1; }
#      echo '************************************'
#      echo '6. ergom sbc done'
	  fi

          #7. open ergom boundary conditions
	  if [ -f $wdir/forcing/ERGOM_OBC_$date_nemo.nc ]; then
	      echo 'Use file ' $wdir/forcing/ERGOM_OBC_$date_nemo.nc
	  else
	      ln -s $archive/ERGOM_OBC_$date_nemo.nc  $wdir/forcing/ERGOM_OBC_$date_nemo.nc || { echo '7. failed' ; exit 1; }
      #$wdir/do_ergom_obc-bounds $date_nemo || { echo '6. failed' ; exit 1; }
#      echo '************************************'
#      echo '7. ergom obc done '
	  fi
	done

	tstr=$(date -I -d "$tstr + 1 day ")
	
        # Count days
	ndays=$(($ndays + 1))
    done

    echo 'Number of days in this run: ' $ndays

    rn_rdt=90
    rnl=$(($ndays * 24 * 3600 / $rn_rdt ))
    nnstep=`printf %08d $rnl`
    echo 'run length rnl: ' $rnl


    echo 'FORCING READY FOR SIMULATION '
    echo '************************************'
    echo '************************************'
# ---------------------------------------------------------------------------------------------------

    # Prepare namelists
    echo 'Prepare namelist files'

    cp $setup_store/ACCESS ./

    ACCESS=`cat ACCESS | head`
    echo 'ACCESS: ' $ACCESS

    if [ "$ACCESS" == ".true." ];then
	echo 'hot start'
	nn_itend=$rnl
	ln_rstart=.true.
	ln_rsttrc=.true.
	ln_trcdta=.false.
	ln_apr_dyn=.true.
	ln_tsd_init=.false.
    else
	echo 'cold start'
	nn_itend=$rnl
	ln_rstart=.false.
	ln_rsttrc=.false.
	ln_trcdta=.false.
#       ln_trcdta=.true.
	ln_apr_dyn=.true.
	ln_tsd_init=.true.
    fi

   #define NEMO_001 output
    FILE101flag=.true. #oce SURF_grid_T
    FILE201flag=.true. #oce grid_T
    FILE301flag=.true. #oce grid_U
    FILE401flag=.true. #oce grid_V
    FILE501flag=.true. #oce grid_W
    FILE601flag=.true. #ice ice_grid_T
    FILE701flag=.true. #ergom _ERGOM_T
   #define NEMO_00X output
    FILE10Xflag=.false. #oce SURF_grid_T
    FILE20Xflag=.false. #oce grid_T
    FILE30Xflag=.false. #oce grid_U
    FILE40Xflag=.false. #oce grid_V
    FILE50Xflag=.false. #oce grid_W
    FILE60Xflag=.false. #ice ice_grid_T
    FILE70Xflag=.false. #.false. #ergom _ERGOM_T

    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      if [ $i -eq 1 ]; then
	  sed -i "s:_FILE1flag_:${FILE101flag}:g" file_def_nemo-oce001.xml
	  sed -i "s:_FILE2flag_:${FILE201flag}:g" file_def_nemo-oce001.xml
	  sed -i "s:_FILE3flag_:${FILE301flag}:g" file_def_nemo-oce001.xml
	  sed -i "s:_FILE4flag_:${FILE401flag}:g" file_def_nemo-oce001.xml
	  sed -i "s:_FILE5flag_:${FILE501flag}:g" file_def_nemo-oce001.xml
	  sed -i "s:_FILE6flag_:${FILE601flag}:g" file_def_nemo-ice001.xml
	  sed -i "s:_FILE7flag_:${FILE701flag}:g" file_def_nemo-ergom001.xml
      else
	  sed -i "s:_FILE1flag_:${FILE10Xflag}:g" file_def_nemo-oce${ENSstr}.xml
	  sed -i "s:_FILE2flag_:${FILE20Xflag}:g" file_def_nemo-oce${ENSstr}.xml
	  sed -i "s:_FILE3flag_:${FILE30Xflag}:g" file_def_nemo-oce${ENSstr}.xml
	  sed -i "s:_FILE4flag_:${FILE40Xflag}:g" file_def_nemo-oce${ENSstr}.xml
	  sed -i "s:_FILE5flag_:${FILE50Xflag}:g" file_def_nemo-oce${ENSstr}.xml
	  sed -i "s:_FILE6flag_:${FILE60Xflag}:g" file_def_nemo-ice${ENSstr}.xml
	  sed -i "s:_FILE7flag_:${FILE70Xflag}:g" file_def_nemo-ergom${ENSstr}.xml
      fi
    done

    echo "Copy configuration files"

    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}
      export wdir
      cp $setup_store/*.nml $wdir/
      cp $setup_store/*.xml $wdir/
      cp iodef.xml $wdir/
      cp $setup_store/namelist* $wdir/
      cp $setup_store/ACCESS $wdir/
      cp `pwd`/namelist_cfg.pdaf $wdir/
      cp `pwd`/*.xml $wdir/

      cat $wdir/namelist_cfg_template		       \
          | sed -e "s:_nn_itend_:$nn_itend:"    \
          | sed -e "s:_ln_rstart_:$ln_rstart:"  \
          | sed -e "s:_nn_date0_:$tstr_ini:"        \
          | sed -e "s:_rn_rdt_:$rn_rdt:"        \
          | sed -e "s:_ln_tsd_init_:$ln_tsd_init:"    \
          | sed -e "s:_cn_ocerst_outdir_:$restart_out:"  \
          | sed -e "s:_ln_apr_dyn_:$ln_apr_dyn:" \
          > $wdir/namelist_ref

      cat $wdir/namelist_top_cfg_template		         \
          | sed -e "s:_ln_rstart_:$ln_rsttrc:"  \
          | sed -e "s:_ln_trcdta_:$ln_trcdta:"       \
          > $wdir/namelist_top_ref
    done

   # Clean up xml prepared xml files
    rm `pwd`/*.xml 


# ---------------------------------------------------------------------------------------------------

    # Create MPMD configuration file

    echo 'Prepare mpmd.conf ...'

    if [ -f mpmd.conf ]; then
	rm mpmd.conf
    fi

    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      #echo 'preparing mpmd.conf...'${ENSstr}

      echo '#!/bin/sh' > nemo${ENSstr}.sh
      echo 'cd '${ENSstr} >> nemo${ENSstr}.sh
      echo './nemo.exe' >> nemo${ENSstr}.sh
      chmod +x nemo${ENSstr}.sh

      echo '#!/bin/sh' > xios${ENSstr}.sh
      echo 'cd '${ENSstr} >> xios${ENSstr}.sh
      echo './xios_server.exe' >> xios${ENSstr}.sh
      chmod +x xios${ENSstr}.sh

      echo $(((i-1)*($np_nemo+$np_xios)))'-'$(((i-1)*($np_nemo+$np_xios)+$np_nemo-1))' ./nemo'${ENSstr}.sh >> mpmd.conf
      #echo $(((i-1)*($np_nemo+$np_xios)+$np_nemo))'-'$(((i)*($np_nemo+$np_xios)-1))' ./xios'${ENSstr}.sh >> mpmd.conf
      echo $(((i-1)*($np_nemo+$np_xios)+$np_nemo))'-'$(((i)*($np_nemo+$np_xios)-1))' ./xios001.sh' >> mpmd.conf
    done

    echo 'JOB PREPARATIONS COMPLETED'

fi  # if $prepare==1

cat mpmd.conf

# Execute the run
if [ $dorun -eq 1 ]; then
    srun -l --cpu_bind=cores --multi-prog mpmd.conf

# ---------------------------------------------------------------------------------------------------
# error check in ocean.output
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}
      export wdir
  
      if grep -q 'E R R O R' $wdir/ocean.output
        then
          { echo 'ERROR in ocean.output' ;  }
        exit 1
      fi
    done

fi # if $dorun==1

#echo "EXIT SCRIPT FOR DEBUGGING"
#exit

if [ $postproc -eq 1 ]; then

    echo "START POSTPROCESSING"
    date

# ---------------------------------------------------------------------------------------------------
# Move output files

    echo "Move NORDIC and DA output files in task 001"
    wdir=`pwd`
    mv $wdir/001/???_NORDIC_* $wdir/001/output/data/
    mv $wdir/001/state_*.nc $wdir/001/output/DA/
    mv $wdir/001/variance_*.nc $wdir/001/output/DA/


# ---------------------------------------------------------------------------------------------------
#Remove old restart files

    echo "clean up directories initialstate/ and forcing/"
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}

      rm $wdir/initialstate/restart_in_*
      rm $wdir/initialstate/restart_ice_in_*
      rm $wdir/initialstate/restart_trc_in_*

    # remove forcing files
      rm $wdir/forcing/bdy_ts_*.nc
      rm $wdir/forcing/bdy_uvh_*.nc
      rm $wdir/forcing/EHYPE_*.nc
      rm $wdir/forcing/FORCE_*.nc
      rm $wdir/forcing/ERGOM_OBC_*.nc
      rm $wdir/forcing/ERGOM_SBC_*.nc
      rm $wdir/forcing/ERGOM_CBC_*.nc
    done


# ---------------------------------------------------------------------------------------------------
#save restarts

    echo "Save restart files"

    np=$(( $np_nemo - 1 ))
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}

      # Save log file
      cp $wdir/ocean.output $wdir/output/log/ocean.output_${tstr_ini}-${tstr}

      # Rename restart files to include date
      for n in `seq -f "%04g" 0 $np`;do
	  mv $wdir/$restart_out/'NORDIC_'$nnstep'_restart_out_'$n'.nc'     $wdir/$restart_out/'/restart_in_'$n'_'$date_nemo'.nc'
	  mv $wdir/$restart_out/'NORDIC_'$nnstep'_restart_ice_out_'$n'.nc' $wdir/$restart_out/'/restart_ice_in_'$n'_'$date_nemo'.nc'
	  mv $wdir/$restart_out/'NORDIC_'$nnstep'_restart_trc_out_'$n'.nc' $wdir/$restart_out/'/restart_trc_in_'$n'_'$date_nemo'.nc'
      done
    done

    # Link restart files for restarting
    echo "Link new restart files into initialstate"
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}

      for n in `seq -f "%04g" 0 $np`;do
	  ln -s $wdir/$restart_out/'/restart_in_'$n'_'$date_nemo'.nc'     $wdir'/initialstate/restart_in_'$n'.nc'
	  ln -s $wdir/$restart_out/'/restart_ice_in_'$n'_'$date_nemo'.nc' $wdir'/initialstate/restart_ice_in_'$n'.nc'
	  ln -s $wdir/$restart_out/'/restart_trc_in_'$n'_'$date_nemo'.nc' $wdir'/initialstate/restart_trc_in_'$n'.nc'
      done
    done

    echo "END POSTPROCESSING"
    date

fi # if postproc=1

# ---------------------------------------------------------------------------------------------------

echo 'JOB DONE'
exit
