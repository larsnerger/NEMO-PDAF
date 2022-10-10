#!/bin/bash
#SBATCH --job-name=NEMOEGM
###SBATCH -p smp
#SBATCH -p mpp
##SBATCH --mem=100G
#SBATCH --ntasks=3840
#SBATCH --ntasks-per-node=96
#SBATCH --time=12:00:00
##SBATCH --partition=standard96:test
#SBATCH --partition=standard96
#SBATCH --mail-user=yuchen.sun@awi.de
#SBATCH --mail-type=END
##SBATCH -A zzz0002

module load intel/18.0.6
module load openmpi/intel/3.1.6
module load netcdf-parallel/ompi/intel.18/4.7.4

ulimit -s unlimited

export LD_LIBRARY_PATH=/sw/dataformats/netcdf-parallel/ompi/intel.18/4.7.4/skl/lib:/sw/dataformats/hdf5-parallel/ompi/intel.18/1.10.6/skl/lib:/sw/compiler/intel/compilers_and_libraries_2018.6.288/linux/compiler/lib/intel64_lin:/sw/compiler/intel/compilers_and_libraries_2018.6.288/linux/mkl/lib/intel64_lin

# ---------------------------------------------------------------------------------------------------

yy=2015 #Start year
mm=01   #Start month
dd=01   #Start day

tstr2='20150101' #End date yyyymmdd

np_nemo=186 #186 102    #number of PE's for nemo
np_xios=6    #number of PE's for xios
np=$np_nemo

# ---------------------------------------------------------------------------------------------------

# Set initial date  $yy$mm$dd==$initial_date we use the restart files from the input directory

initial_date=20150101


# ---------------------------------------------------------------------------------------------------

NENS=20  # number of ensembles

# ---------------------------------------------------------------------------------------------------

#set -xv #debugging

# ---------------------------------------------------------------------------------------------------
nemo_exe_dir='/home/hbkycsun/nemo4_dev-nemo4_ergom_allEuler/cfgs/NEMO_ERGOM_PDAF/BLD/bin'
#nemo_exe_dir='/home/hzfblner/SEAMLESS/nemo4_dev_allEuler/cfgs/NEMO-ERGOM-PDAF/BLD/bin'
#nemo_exe_dir='/home/hzfblner/SEAMLESS/nemo4_dev_allEuler/cfgs/NEMO_ERGOM_run/BLD/bin'
#nemo_exe_dir='/home/hbkycsun/nemo4_dev/cfgs/NEMO_ERGOM_ens/BLD/bin'
#nemo_exe_dir='/home/hbkycsun/nemo4_ergom_allEuler/cfgs/NEMO_ERGOM_ens/BLD/bin'
xios_exe_dir='/home/hzfblner/SEAMLESS/xios-2.0_par/bin'
#rebuild='/home_ad/bm1405/balmfc_git/nemo4_dev/tools/REBUILD_NEMO/'
restart_out='output/restarts'
#archive='/work/ollie/fdaryabo/ergom_data'
archive='/scratch/usr/hzfblner/SEAMLESS/forcing_ergom_allEuler'
archive_ln='/scratch/usr/hzfblner/SEAMLESS/forcings'
#inputs_ln='/scratch/usr/hzfblner/SEAMLESS/inputs'
inputs_ln='/scratch/usr/hzfblner/SEAMLESS/inputs_allEuler'
inputs_nc='/scratch/usr/hzfblner/SEAMLESS/run/Eens/inputs_nc'
setup_store='/scratch/usr/hbkycsun/data/setup_store'
setup_store_allEular='/scratch/usr/hbkycsun/data/setup_ERGOM_allEuler'
initialdir='/scratch/usr/hzfblner/SEAMLESS/restart'
disrestartdir='/scratch/usr/hzfblner/SEAMLESS/run/Eens/restart_20150101'
# ---------------------------------------------------------------------------------------------------

# Prepare PDAF namelist
cp $setup_store_allEular/namelist_cfg.pdaf_template ./
cat namelist_cfg.pdaf_template     \
   | sed -e "s:_DIMENS_:$NENS:"     \
   > namelist_cfg.pdaf


for((i=1;i<=$NENS;i++))
   do
     ENSstr=`printf %03d $i`
     ENSstr2=`printf %02d $i`
     echo 'preparing context_nemo.xml...'
     cat $setup_store_allEular/context_nemoXXX.xml_template   \
	    | sed -e "s:_DIMENS3_:${ENSstr}:g"   \
	    > context_nemo${ENSstr}.xml
     echo 'preparing file_def_nemo-ice.xml...'
	   cat $setup_store_allEular/file_def_nemo-iceXXX.xml_template   \
	    | sed -e "{s:_DIMENS3_:${ENSstr}:g;s:_DIMENS2_:${ENSstr2}:g}"   \
	    > file_def_nemo-ice${ENSstr}.xml
     echo 'preparing file_def_nemo-oce.xml...'
	   cat $setup_store_allEular/file_def_nemo-oceXXX.xml_template   \
	    | sed -e "{s:_DIMENS3_:${ENSstr}:g;s:_DIMENS2_:${ENSstr2}:g}"   \
	    > file_def_nemo-oce${ENSstr}.xml
     echo 'preparing file_def_nemo-ergom.xml...'
     cat $setup_store_allEular/file_def_nemo-ergomXXX.xml_template   \
 	    | sed -e "{s:_DIMENS3_:${ENSstr}:g;s:_DIMENS2_:${ENSstr2}:g}"   \
 	    > file_def_nemo-ergom${ENSstr}.xml
done

# set working directories for different ensemble members
for((i=1;i<=$NENS;i++))
  do
    ENSstr=`printf %03d $i`
    echo 'creating ensemble working directories...'
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
    mkdir -p $wdir/initialstate
    mkdir -p $wdir/forcing

    echo ' '
    echo 'Run directory: ' $wdir
done

# ---------------------------------------------------------------------------------------------------
#export wdir
export nemo_exe_dir
#export rebuild
export restart_out
export archive
# ---------------------------------------------------------------------------------------------------
#linking executables
#if [ ! -f xios_server.exe ]; then
#  ln -s $xios_exe_dir/xios_server.exe
#fi
#if [ ! -f nemo.exe ]; then
#  ln -s $nemo_exe_dir/nemo.exe
#fi
# ---------------------------------------------------------------------------------------------------
for((i=1;i<=$NENS;i++))
  do
    ENSstr=`printf %03d $i`
    wdir=`pwd`/${ENSstr}
    export wdir
    echo 'linking executables...'
    if [ ! -f ${ENSstr}/xios_server.exe ]; then
      ln -s $xios_exe_dir/xios_server.exe $wdir/
    fi
    if [ ! -f ${ENSstr}/nemo.exe ]; then
      ln -s $nemo_exe_dir/nemo.exe $wdir/
    fi
done
# ---------------------------------------------------------------------------------------------------
tstr0=`date +%Y%m%d`

rn_rdt=90
rnl=$((24 * 3600 / $rn_rdt ))
if [ $rn_rdt -ge 90 ]; then
   nnstep=00000$rnl #00000 depends on rnl  -- for 90
else
   nnstep=0000$rnl #00000 depends on rnl   --- for 60
fi
echo 'run length rnl: ' $rnl
sdte=`date --date "$dte0 0 day" +'%Y-%m-%d %H:%M:%S'`
yyp1=`date +'%Y' -d"$sdte"`                     #  formating
mmp1=`date +'%m' -d"$sdte"`
ddp1=`date +'%d' -d"$sdte"`
tstr="$yy$mm$dd"
#echo 'fast: -fp-model precise -Dmpi -O3 -r8'               > log.$tstr0
echo 'Simulation with nemo executable from ' $nemo_exe_dir
# ---------------------------------------------------------------------------------------------------

#linking restart

for((i=1;i<=$NENS;i++))
  do
    ENSstr=`printf %03d $i`
    wdir=`pwd`/${ENSstr}
    export wdir
    echo 'linking restart...'
    if [ $yy$mm$dd -eq $initial_date ]; then
      echo "Link initial restart files"
      if [ $restart_dis -eq 0 ]; then
        echo "Link integrated restart files"
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
           ln -s $inputs_ln/initialstate/NORDIC-NS1_restart_ptrc_in_iowfix.nc $wdir/initialstate/NORDIC-NS1_restart_ptrc_in_iowfix.nc
        fi
        if [ ! -f $wdir/initialstate/NORDIC-NS1_restart_ptrc_in.nc ]; then
           ln -s $inputs_ln/initialstate/NORDIC-NS1_restart_ptrc_in.nc $wdir/initialstate/NORDIC-NS1_restart_ptrc_in.nc
        fi
      else
        echo "Link distributed restart files"
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

#if [ ! -f $wdir/forcing/river_data.nc ]; then
#    ln -s $inputs_ln/river_data.nc $wdir/forcing
#fi
#if [ ! -f $wdir/forcing/dum12_y2015.nc ]; then
#ln -s $inputs_ln/dum12_y2015.nc $wdir/forcing
#fi
#if [ ! -f $wdir/forcing/weights_bilin.nc ]; then
#ln -s $inputs_ln/weights_bilin.nc $wdir/forcing
#fi
#if [ ! -f $wdir/forcing/weights_bicubic.nc ]; then
#ln -s $inputs_ln/weights_bicubic.nc $wdir/forcing
#fi

for((i=1;i<=$NENS;i++))
  do
    ENSstr=`printf %03d $i`
    wdir=`pwd`/${ENSstr}
    export wdir
    echo 'linking input...'
    if [ ! -f $wdir/forcing/river_data.nc ]; then
        ln -s $inputs_nc/river_data.nc $wdir/forcing
        ln -s $inputs_nc/river_data.nc $wdir/
    fi
    if [ ! -f $wdir/forcing/dum12_y2015.nc ]; then
        ln -s $inputs_nc/dum12_y2015.nc $wdir/forcing
        ln -s $inputs_nc/dum12_y2015.nc $wdir/
    fi
    if [ ! -f $wdir/forcing/weights_bilin.nc ]; then
        ln -s $inputs_nc/weights_bilin.nc $wdir/forcing
        ln -s $inputs_nc/weights_bilin.nc $wdir/
    fi
    if [ ! -f $wdir/forcing/weights_bicubic.nc ]; then
        ln -s $inputs_nc/weights_bicubic.nc $wdir/forcing
        ln -s $inputs_nc/weights_bicubic.nc $wdir/
    fi

    ln -s $inputs_nc/dum12_y2015.nc         $wdir/
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
    ln -s $setup_store_allEular/*.txt       $wdir/

done

#Loop over dates

while [ "$(date -d "$tstr" +%Y%m%d)" -le "$(date -d "$tstr2" +%Y%m%d)" ]; do
   date_nemo=$(date -d "$tstr" +y%Ym%md%d)
   date_force=$(date -d "$tstr" +%Y%m%d)
   date_NS01=$(date -d "$tstr" +y%Ym%m)
   tstr=$(date -d "$tstr" +%Y%m%d)
#   echo 'tstr: ' $tstr
#   echo 'tstr2: ' $tstr2
#   echo 'date_nemo: ' $date_nemo
#   echo 'date_NS01: ' $date_NS01

for((i=1;i<=$NENS;i++))
  do
    ENSstr=`printf %03d $i`
    wdir=`pwd`/${ENSstr}
    export wdir
    echo 'Preparing input...'

   #1. get river forcing
   if [ -f $wdir/forcing/EHYPE_$date_nemo.nc ]; then
      echo 'Use file ' $wdir/forcing/EHYPE_$date_nemo'+024H.nc'
   else
      ln -s $archive_ln/EHYPE_$date_nemo'+024H.nc' $wdir/forcing/EHYPE_$date_nemo.nc || { echo '1. failed' ; exit 1; }
      echo '************************************'
      echo '1. river forcing done'
   fi

   #2. get atm. force
   if [ -f $wdir/forcing/FORCE_$date_nemo.nc ]; then
      echo 'Use file ' $wdir/forcing/FORCE_$tstr'+24.nc'
   else
      ln -s $archive_ln/FORCE_$tstr'+24.nc' $wdir/forcing/FORCE_$date_nemo.nc  #|| { echo '2. failed' ; exit 1; }
      echo '************************************'
      echo '2. atm. forcing done'
   fi

   #3a. get bdy
   if [ -f $wdir/forcing/bdy_uvh_$date_nemo.nc ]; then
      echo 'Use file ' $wdir/forcing/bdy_uvh_$tstr.nc
   else
      ln -s $archive_ln/bdy_uvh_$tstr.nc  $wdir/forcing/bdy_uvh_$date_nemo.nc || { echo '3. failed' ; exit 1; }
      echo '************************************'
      echo '3. bdy done'
   fi

   #4. get bdy
   if [ -f $wdir/forcing/bdy_ts_$date_nemo.nc ]; then
      echo 'Use file ' $wdir/forcing/bdy_ts_$tstr.nc
   else
      ln -s $archive_ln/bdy_ts_$tstr.nc  $wdir/forcing/bdy_ts_$date_nemo.nc || { echo '4. failed' ; exit 1; }
#	ln -s /cmems_archive/bm1405/force_nemo4_braxen_2nm/bdy_test/bdy_ts_$tstr.nc  $wdir/forcing/bdy_ts_$date_nemo.nc || { echo '4. failed' ; exit 1; }
      echo '************************************'
      echo '4. obc done'
   fi

   #5. ERGOM river loads
   if [ -f $wdir/forcing/ERGOM_CBC_$date_nemo.nc ]; then
      echo 'Use file ' $wdir/forcing/ERGOM_CBC_$date_nemo.nc
   else
      ln -s $archive/ERGOM_CBC_$date_nemo.nc  $wdir/forcing/ERGOM_CBC_$date_nemo.nc || { echo '5. failed' ; exit 1; }
	    #$wdir/do_ergom_cbc-rivers $date_nemo || { echo '5. failed' ; exit 1; }
      echo '************************************'
      echo '5. ergom river loads done'
   fi

   #6. ergom surface bdy conditions
   if [ -f $wdir/forcing/ERGOM_SBC_$date_nemo.nc ]; then
      echo 'Use file ' $wdir/forcing/ERGOM_SBC_$date_nemo.nc
   else
      ln -s $archive/ERGOM_SBC_$date_nemo.nc  $wdir/forcing/ERGOM_SBC_$date_nemo.nc || { echo '6. failed' ; exit 1; }
      #$wdir/do_ergom_sbc-depo $date_nemo || { echo '5. failed' ; exit 1; }
      echo '************************************'
      echo '6. ergom sbc done'
   fi

   #7. open ergom boundary conditions
   if [ -f $wdir/forcing/ERGOM_OBC_$date_nemo.nc ]; then
      echo 'Use file ' $wdir/forcing/ERGOM_OBC_$date_nemo.nc
   else
      ln -s $archive/ERGOM_OBC_$date_nemo.nc  $wdir/forcing/ERGOM_OBC_$date_nemo.nc || { echo '7. failed' ; exit 1; }
      #$wdir/do_ergom_obc-bounds $date_nemo || { echo '6. failed' ; exit 1; }
      echo '************************************'
      echo '7. ergom obc done '
   fi

    #8. ERGOM TCC
#    if [ -f $wdir/forcing/ERGOM_TCC_$date_nemo.nc ]; then
#	echo 'Use file ' $wdir/forcing/ERGOM_TCC_$date_nemo.nc
#    else
#	ln -s $archive/ERGOM_TCC_$date_nemo.nc  $wdir/forcing/ERGOM_TCC_$date_nemo.nc || { echo '8. failed' ; exit 1; }
#	echo '************************************'
#	echo '8. ergom tcc done '
#    fi

done

   echo 'FORCING READY FOR SIMULATION '
   echo '************************************'
   echo '************************************'
# ---------------------------------------------------------------------------------------------------
# Prepare namelists
   echo 'Prepare namelist files'

   cp $setup_store_allEular/ACCESS ./

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
   FILE70Xflag=.false. #ergom _ERGOM_T

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



   for((i=1;i<=$NENS;i++))
     do
       ENSstr=`printf %03d $i`
       wdir=`pwd`/${ENSstr}
       export wdir
       cp $setup_store_allEular/*.nml $wdir/
       cp $setup_store_allEular/*.xml $wdir/
       cp $setup_store_allEular/namelist* $wdir/
       cp $setup_store_allEular/ACCESS $wdir/
       cp $setup_store/ens_dates.txt $wdir/
       #cp $setup_store/files*.txt $wdir/
       cp `pwd`/namelist_cfg.pdaf $wdir/
       cp `pwd`/*.xml $wdir/
       ##sed -e 's/nemo_000/nemo_${ENSstr}/' $wdir/../context_nemo_000.xml >$wdir/context_nemo.xml
       #ln -s /scratch/usr/hzfblner/SEAMLESS/out_free/NORDIC_1d_SURF_grid_T_20150115-20150115.nc $wdir/my_nemo_ssh_file.nc
       cat $wdir/namelist_cfg_template		       \
          | sed -e "s:_nn_itend_:$nn_itend:"    \
          | sed -e "s:_ln_rstart_:$ln_rstart:"  \
          | sed -e "s:_nn_date0_:$tstr:"        \
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

   #cat namelist_cfg_template		         \
   #| sed -e "s:_nn_itend_:$nn_itend:"    \
   #| sed -e "s:_ln_rstart_:$ln_rstart:"  \
   #| sed -e "s:_nn_date0_:$tstr:"        \
   #| sed -e "s:_rn_rdt_:$rn_rdt:"        \
   #| sed -e "s:_ln_tsd_init_:$ln_tsd_init:"    \
   #| sed -e "s:_cn_ocerst_outdir_:$restart_out:"  \
   #| sed -e "s:_ln_apr_dyn_:$ln_apr_dyn:" \
   #> namelist_ref
#
   #cat namelist_top_cfg_template		         \
   #| sed -e "s:_ln_rstart_:$ln_rsttrc:"  \
   #| sed -e "s:_ln_trcdta_:$ln_trcdta:"       \
   #> namelist_top_ref

# ---------------------------------------------------------------------------------------------------
#start model run


# Create MPMD configuration file
if [ -f mpmd.conf ]; then
  rm mpmd.conf
fi

for((i=1;i<=$NENS;i++))
do
      ENSstr=`printf %03d $i`
      echo 'preparing mpmd.conf...'${ENSstr}

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

cat mpmd.conf

   #(( NCPU = np_nemo + np_xios ))
   #(( NEMO_LAST_CPU = np_nemo - 1 ))
   #(( XIOS_LAST_CPU = NCPU  - 1 ))
   #echo 0-${NEMO_LAST_CPU} ./nemo.exe > mpmd.conf
   #echo ${np_nemo}-${XIOS_LAST_CPU} ./xios_server.exe >> mpmd.conf
   #cat mpmd.conf
#
   #if [ -z "$np" ]; then np=1;fi
   #if [ -z "$xios" ]; then xios=1;fi
   #echo 'Start the run with ' $np_nemo ' PES for NEMO and ' $np_xios 'for xios #server'

   srun -l --cpu_bind=cores --multi-prog mpmd.conf
#exit
# ---------------------------------------------------------------------------------------------------
# error check in ocean.output
for((i=1;i<=$NENS;i++))
  do
    ENSstr=`printf %03d $i`
    wdir=`pwd`/${ENSstr}
    export wdir

    if grep -q 'E R R O R' $wdir/ocean.output
      then
      #{ echo 'ERROR in ocean.output' ;./test_err.sh; exit 1; }  >> log.$tstr0
        { echo 'ERROR in ocean.output' ;  }
        exit 1
    fi
done

#echo "EXIT SCRIPT FOR DEBUGGING"
#exit

# ---------------------------------------------------------------------------------------------------
# Move output files
#   mv $wdir/NORDIC_* $wdir/output/data/
#   mv $wdir/station_* $wdir/output/data/

# call compression backround (and rebuild)
#./compress.sh $date_nemo $wdir/output/data/  &>> log.$tstr0
for((i=1;i<=$NENS;i++))
  do
    ENSstr=`printf %03d $i`
    wdir=`pwd`/${ENSstr}
    export wdir
    mv $wdir/../001/${ENSstr}_NORDIC_* $wdir/output/data/
    #mv $wdir/NORDIC_1h* $wdir/output/data/
    #mv $wdir/NORDIC_6h* $wdir/output/data/
    #mv $wdir/NORDIC_1d* $wdir/output/data/
    #mv $wdir/station_* $wdir/output/data/
    #mv $wdir/NORDIC_1ts_SURF_grid_T_2* $wdir/output/data/
    # call compression backround (and rebuild)
    #./compress.sh $date_nemo $wdir/output/data/  &>> log.$tstr0
done

# ---------------------------------------------------------------------------------------------------

#cp $wdir'/initialstate/restart_in.nc'     $wdir'/initialstate/restart_in'$date_nemo'.nc'                              >> log.$tstr0
#cp $wdir'/initialstate/restart_ice_in.nc' $wdir'/initialstate/restart_ice_in'$date_nemo'.nc'                          >> log.$tstr0
#cp $wdir'/initialstate/restart_trc_in.nc' $wdir'/initialstate/restart_trc_in'$date_nemo'.nc'                          >> log.$tstr0
#gzip $wdir'/initialstate/restart_'$tstr0'.nc'
# ---------------------------------------------------------------------------------------------------
#cd $wdir/$restart_out/
#ln -s $rebuild/rebuild_nemo.exe
#for n in `seq 1 3`;do
#	cat $wdir/$restart_out/nam_rebuild$n'_template'		         \
#	| sed -e "s:_ndomain_:$np:"  \
#	> $wdir/$restart_out/nam_rebuild$n
#	./rebuild_nemo.exe nam_rebuild$n  >> log.$tstr0
#done
##

# ---------------------------------------------------------------------------------------------------
#Remove old restart files
   #echo "Remove previous restart files from initialstate"
   #rm -f ${wdir}/initialstate/restart_in*
   #rm -f ${wdir}/initialstate/restart_ice_in*
   #rm -f ${wdir}/initialstate/restart_trc_in*
   for((i=1;i<=$NENS;i++))
     do
       ENSstr=`printf %03d $i`
       wdir=`pwd`/${ENSstr}
       export wdir
       echo "Remove old restart files from initialstate"
       rm $wdir'/initialstate/restart_in_*'
       rm $wdir'/initialstate/restart_ice_in_*'
       rm $wdir'/initialstate/restart_trc_in_*'
   done


# ---------------------------------------------------------------------------------------------------
#save restarts
   np=$(( $np_nemo - 1 ))
   for((i=1;i<=$NENS;i++))
     do
       ENSstr=`printf %03d $i`
       wdir=`pwd`/${ENSstr}
       export wdir

       #save restarts


   echo "Save restart files"

   for n in `seq -f "%04g" 0 $np`;do
      mv $wdir/$restart_out/'NORDIC_'$nnstep'_restart_out_'$n'.nc'     $wdir'/initialstate/restart_in_'$n'.nc'
      mv $wdir/$restart_out/'NORDIC_'$nnstep'_restart_ice_out_'$n'.nc' $wdir'/initialstate/restart_ice_in_'$n'.nc'
      mv $wdir/$restart_out/'NORDIC_'$nnstep'_restart_trc_out_'$n'.nc' $wdir'/initialstate/restart_trc_in_'$n'.nc'
   done

   for n in `seq -f "%04g" 0 $np`;do
      cp $wdir'/initialstate/restart_in_'$n'.nc'     $wdir/$restart_out/'/restart_in_'$n'_'$date_nemo'.nc'
      cp $wdir'/initialstate/restart_ice_in_'$n'.nc' $wdir/$restart_out/'/restart_ice_in_'$n'_'$date_nemo'.nc'
      cp $wdir'/initialstate/restart_trc_in_'$n'.nc' $wdir/$restart_out/'/restart_trc_in_'$n'_'$date_nemo'.nc'
   done

#   for n in `seq -f "%04g" 0 $np`;do
#      cp $wdir/$restart_out/'NORDIC_'$nnstep'_restart_out_'$n'.nc'     $wdir'/initialstate/restart_in_'$n'.nc'
#      cp $wdir/$restart_out/'NORDIC_'$nnstep'_restart_ice_out_'$n'.nc' $wdir'/initialstate/restart_ice_in_'$n'.nc'
#      cp $wdir/$restart_out/'NORDIC_'$nnstep'_restart_trc_out_'$n'.nc' $wdir'/initialstate/restart_trc_in_'$n'.nc'
#   done

#   for n in `seq -f "%04g" 0 $np`;do
#      if [ -f $wdir'/initialstate/restart_in_'$n'.nc' ] && [ $wdir'/initialstate/restart_ice_in_'$n'.nc' ] && [ $wdir'/initialstate/restart_trc_in_'$n'.nc' ]; then
#	rm $wdir/$restart_out/'NORDIC_'$nnstep'_restart_out_'$n'.nc'
#	rm $wdir/$restart_out/'NORDIC_'$nnstep'_restart_ice_out_'$n'.nc'
#	rm $wdir/$restart_out/'NORDIC_'$nnstep'_restart_trc_out_'$n'.nc'
#      else
#	{ echo ' failed' ; exit 1; }
#      fi
#   done


#for n in `seq -f "%04g" 0 $np`;do
#rm $wdir/$restart_out/'NORDIC_'$nnstep'_restart_out_'$n'.nc'     >> log.$tstr0
#rm $wdir/$restart_out/'NORDIC_'$nnstep'_restart_ice_out_'$n'.nc' >> log.$tstr0
#rm $wdir/$restart_out/'NORDIC_'$nnstep'_restart_trc_out_'$n'.nc' >> log.$tstr0
#done


   cp $wdir/ocean.output $wdir/output/log/ocean.output$tstr
   cp $wdir/timing.output $wdir/output/log/timing.output$tstr
#mv $wdir/output/ /cmems_shared_work/bm1405/nemo_output/output_$date_nemo

# ---------------------------------------------------------------------------------------------------
# remove forcing files

   rm $wdir/forcing/bdy_ts_$date_nemo'.nc'
   rm $wdir/forcing/bdy_uvh_$date_nemo'.nc'
   rm $wdir/forcing/EHYPE_$date_nemo'.nc'
   rm $wdir/forcing/FORCE_$date_nemo'.nc'
#rm $wdir/forcing/ERGOM_TCC_$date_nemo'.nc'
   rm $wdir/forcing/ERGOM_OBC_$date_nemo'.nc'
   rm $wdir/forcing/ERGOM_SBC_$date_nemo'.nc'
   rm $wdir/forcing/ERGOM_CBC_$date_nemo'.nc'
done

# ---------------------------------------------------------------------------------------------------
# add one day
   tstr=$(date -I -d "$tstr + 1 day ")
   echo 'next date: '$tstr

   echo '.true.' > ACCESS
done

echo 'JOB DONE'
exit
