#!/bin/bash -l

params_file=params_CTL.txt
#params_file=params_pop.txt

ramm_file=RAMM.out

# 1 - Choose models to simulate
# Run the population of 16 control models with <params_CTL.txt> file
indices=( $(seq 1 17 ) )

# 2 - Run the population of 16 control models with <params_pop.txt> file
#indices=( 1 171 268 452 518 759 946 959 1105 1172 1195 1391 1639 1650 1711 1822 2065 )

# 3 - Run the entire population of 3000 models with the <params_pop.txt> file
#indices=( $(seq 1 3000 ) )



for model in "${indices[@]}"
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

  $ramm_file Parameters_$model.txt output_$model.txt CELL

done

