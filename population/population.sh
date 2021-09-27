#!/bin/bash -l
#SBATCH --job-name=PoMCTL
#SBATCH --account=nn9249k
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --output=out.txt
#SBATCH --error=err.txt

source /cluster/bin/jobsetup

module load openmpi.intel/1.7
MPIRUN=mpirun

module load zlib

array=( 484 510 652 880 928 946 959 1044 1067 1086 1105 1162 1172 1195 1252 1306 1340 1391 1639 1799 1928 2029 2266 2409 2656 2674 2794 2797 2828 )

params_file=/work/users/marcia/rabbit/source_v6/PoM_2Hz_3/params_CTL.txt
ramm_file=/work/users/marcia/rabbit/source_v6/RAMM.out

for model in "${array[@]}"
do
  cp /work/users/marcia/rabbit/source_v6/PoM_2Hz_3/Parameters_4.txt ./Parameters_$model.txt
  cp /work/users/marcia/rabbit/source_v6/PoM_2Hz_3/run.sh ./
  p=(`awk -v line="$model" 'NR==line' $params_file`)    

  sed -i 's/ICaL_junc_scl=ICaL_junc_scl/ICaL_junc_scl='"${p[0]}"'/' Parameters_$model.txt
  sed -i 's/ICaL_sl_scl=ICaL_sl_scl/ICaL_sl_scl='"${p[0]}"'/' Parameters_$model.txt
  sed -i 's/ICaT_junc_scl=ICaT_junc_scl/ICaT_junc_scl='"${p[1]}"'/' Parameters_$model.txt
  sed -i 's/ICaT_sl_scl=ICaT_sl_scl/ICaT_sl_scl='"${p[1]}"'/' Parameters_$model.txt
  sed -i 's/INaK_junc_scl=INaK_junc_scl/INaK_junc_scl='"${p[2]}"'/' Parameters_$model.txt
  sed -i 's/INaK_sl_scl=INaK_sl_scl/INaK_sl_scl='"${p[2]}"'/' Parameters_$model.txt
  sed -i 's/INaCa_junc_scl=INaCa_junc_scl/INaCa_junc_scl='"${p[3]}"'/' Parameters_$model.txt
  sed -i 's/INaCa_sl_scl=INaCa_sl_scl/INaCa_sl_scl='"${p[3]}"'/' Parameters_$model.txt
  sed -i 's/INa_junc_scl=INa_junc_scl/INa_junc_scl='"${p[4]}"'/' Parameters_$model.txt
  sed -i 's/INa_sl_scl=INa_sl_scl/INa_sl_scl='"${p[4]}"'/' Parameters_$model.txt
  sed -i 's/ITo_scl=ITo_scl/ITo_scl='"${p[5]}"'/' Parameters_$model.txt
  sed -i 's/IKr_scl=IKr_scl/IKr_scl='"${p[6]}"'/' Parameters_$model.txt
  sed -i 's/IKs_scl=IKs_scl/IKs_scl='"${p[7]}"'/' Parameters_$model.txt
  sed -i 's/IK1_scl=IK1_scl/IK1_scl='"${p[8]}"'/' Parameters_$model.txt
  sed -i 's/ICaP_junc_scl=ICaP_junc_scl/ICaP_junc_scl='"${p[9]}"'/' Parameters_$model.txt
  sed -i 's/ICaP_sl_scl=ICaP_sl_scl/ICaP_sl_scl='"${p[9]}"'/' Parameters_$model.txt
  sed -i 's/ICab_junc_scl=ICab_junc_scl/ICab_junc_scl='"${p[10]}"'/' Parameters_$model.txt
  sed -i 's/ICab_sl_scl=ICab_sl_scl/ICab_sl_scl='"${p[10]}"'/' Parameters_$model.txt
  sed -i 's/INab_junc_scl=INab_junc_scl/INab_junc_scl='"${p[11]}"'/' Parameters_$model.txt
  sed -i 's/INab_sl_scl=INab_sl_scl/INab_sl_scl='"${p[11]}"'/' Parameters_$model.txt
  sed -i 's/IClb_scl=IClb_scl/IClb_scl='"${p[12]}"'/' Parameters_$model.txt
  sed -i 's/X0_File=states.txt/X0_File=states_'"$model"'.txt/' Parameters_$model.txt

  #sed -i 's/Parameters.txt/Parameters_'"$model"'.txt/' run.sh
  #sed -i 's/out.txt/out_'"$model"'.txt/' run.sh
  #sed -i 's/outputfile/output_'"$model"'.txt/' run.sh

  #sbatch run.sh
  $ramm_file Parameters_$model.txt output_$model.txt CELL

done

