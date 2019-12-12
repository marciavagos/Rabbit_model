#!/bin/bash -l


ramm=/Users/marciavagos/repos/Rabbit_model/source_v7/RAMM_v6.out

bcl=500
startsave=115000
ft=120000
Istim=0
Nai_initial=6.86

array=( 171 268 452 518 759 946 959 1105 1172  1195 1391 1639 1650 1711 1822 2065 )

for model in "${array[@]}"
do
  cp /Users/marciavagos/repos/Rabbit_model/population_2/Parameters_$model.txt ./

  sed -i '' 's/bcl=bcl/bcl='"$bcl"'/' Parameters_$model.txt
  sed -i '' 's/startsave=startsave/startsave='"$startsave"'/' Parameters_$model.txt
  sed -i '' 's/ft=ft/ft='"$ft"'/' Parameters_$model.txt
  sed -i '' 's/Istim=Istim/Istim='"$Istim"'/' Parameters_$model.txt
  sed -i '' 's/Nai_initial=Nai_initial/Nai_initial='"$Nai_initial"'/' Parameters_$model.txt

  $ramm Parameters_$model.txt output CELL
  mv output_"$bcl"_ms.txt quiescent_"$model".txt
  mv output_"$bcl"_ms_APD.txt quiescent_"$model"_APD.txt
  mv states.txt states_$model.txt
done
