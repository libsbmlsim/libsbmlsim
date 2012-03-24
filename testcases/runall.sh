#!/bin/bash

### directories and files
# ~/Downloads/SBML-testcases/cases/semantic/00980/
#   00980-settings.txt
#   00980-sbml-l2v4.xml
#   00980-results.csv

BaseDir="./cases/semantic"

### 00001-settings.txt
# duration: 5
# steps: 50
# variables: S1, S2, S3, S4
# amount: S1, S2
# concentration:

# Simulate following models with small dt
fine_delta="00952 00953 00962 00963 00964 00965 00966 00967"

for i in {0000{1..9},000{10..99},00{100..980}}; do
#for i in 00957 00958 00959; do
  head="$BaseDir/$i/$i"
  # 00926 and 00927 contains '\r' before '\n' on each line (DOS format),
  # so we have to call tr -d '\r' before parsing the settings file...
  duration=`tr -d '\r' < "$head-settings.txt" | grep duration | cut -d" " -f 2`
  steps=`tr -d '\r' < "$head-settings.txt" | grep steps | cut -d" " -f 2`
  variables=`tr -d '\r' < "$head-settings.txt" | grep variables | cut -d":" -f 2 | sed -e "s/ //g"`
  amount=`tr -d '\r' < "$head-settings.txt" | grep amount | cut -d":" -f 2 | sed -e "s/ //g"`
  concentration=`tr -d '\r' < "$head-settings.txt" | grep amount | cut -d":" -f 2 | sed -e "s/ //g"`
  #sbml="$head-sbml-l2v4.xml"
  sbml="$head-sbml-l3v1.xml"
  result="$i-results.csv"
  opt_delta=""
  opt_amount=""
  echo "$i: $duration : $steps : [$variables] : [$amount] : [$concentration]"
  if [[ "$fine_delta" == *"$i"* ]]; then
    echo "  simulate with fine delta"
    opt_delta="-d 0.00001"
  fi
  if [ -n "$amount" ] ; then
	  echo "  print amount"
	  opt_amount="-a"
  else
	  echo "  print concentration"
  fi
  ./sim -t $duration -s $steps $opt_delta -m 1 -n $opt_amount $sbml && \
  ./genresult.pl out.csv $variables $steps > $result
  ./compare.pl $i
  
  unset head
  unset duration
  unset steps
  unset variables
  unset sbml
  unset result
  unset opt_delta
  unset opt_amount
# end
done
unset BaseDir
unset fine_delta
rm -f out.csv
#zip result.zip 00*-results.csv
#rm 00*-results.csv
