#!/bin/bash -l

startsave=115000
ft=120000 # 2 min pacing
Istim=-20
Nai_initial=10

ramm=/Users/marciavagos/repos/Rabbit_model/source_v7/RAMM_v6.out

array=( 268 452 518 759 946 959 1105 1172  1195 1391 1639 1650 1711 1822 2065 )

for model in "${array[@]}"
do

  mkdir -p model_$model
  cd model_$model

  for bcl in 5000 2000 1000 500 333
  do
    cp /Users/marciavagos/repos/Rabbit_model/population_2/Parameters_$model.txt ./

    sed -i '' 's/bcl=bcl/bcl='"$bcl"'/' Parameters_$model.txt
    sed -i '' 's/startsave=startsave/startsave='"$startsave"'/' Parameters_$model.txt
    sed -i '' 's/ft=ft/ft='"$ft"'/' Parameters_$model.txt
    sed -i '' 's/Istim=Istim/Istim='"$Istim"'/' Parameters_$model.txt
    sed -i '' 's/Nai_initial=Nai_initial/Nai_initial='"$Nai_initial"'/' Parameters_$model.txt

    $ramm Parameters_$model.txt output CELL
    mv states.txt states_$rate.txt

  done
  cd ..
done
