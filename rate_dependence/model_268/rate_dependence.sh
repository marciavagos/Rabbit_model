#!/bin/bash -l

bcl=500
startsave=115000
ft=120000
Istim=-20

model=/Users/marciavagos/repos/Rabbit_model/source_v6/RAMM_v6.out


for rate in 5000 2000 1000 500 333
do
  cp /Users/marciavagos/repos/Rabbit_model/rate_dependence/Parameters_RD.txt ./

  sed -i '' 's/bcl=bcl/bcl='"$rate"'/' Parameters_RD.txt
  sed -i '' 's/startsave=startsave/startsave='"$startsave"'/' Parameters_RD.txt
  sed -i '' 's/ft=ft/ft='"$ft"'/' Parameters_RD.txt
  sed -i '' 's/Istim=Istim/Istim='"$Istim"'/' Parameters_RD.txt

  $model Parameters_RD.txt output CELL
  mv states.txt states_$rate.txt
  #i=$(($i+1))
done
