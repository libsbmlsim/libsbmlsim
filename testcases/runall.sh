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
# absolute: 1.000000e-007
# relative: 0.0001
# amount: S1, S2
# concentration:

# Simulate following models with small dt
fine_delta="00952 00953 00962 00963 00964 00965 00966 00967"

# Simulate following models with small atol
fine_atol="00430 00431 00863 00893 00945 00946 00947 00948"

# Simulate following models with small rtol
fine_rtol="00944 00945 00946 00947 00948"

# Simulate following models with small facmax
fine_facmax="00426 00430 00431 00863"

# for testcase 00001 .. 00980
for i in {0000{1..9},000{10..99},00{100..980}}; do
#for i in 00952 00953; do
#for i in 00266 00267 00268 00893; do # check for default atol
#for i in 00430 00431 00863 00893 00945 00946 00947 00948; do # check for fine atol
#for i in 00952 00962 00964 00965; do # check for fine delta

# for testcase 00001 .. 01123
#for i in {0000{1..9},000{10..99},00{100..999},0100{0..9},010{10..99},011{00..23}}; do
#for i in {00{965..999},0100{1..9},010{10..99},011{00..23}}; do

# for testcase 00981 .. 01123
#for i in {00{981..999},0100{0..9},010{10..99},0110{0..9},011{10..23}}; do

  head="$BaseDir/$i/$i"
  # 00926 and 00927 contains '\r' before '\n' on each line (DOS format),
  # so we have to call tr -d '\r' before parsing the settings file...
  duration=`tr -d '\r' < "$head-settings.txt" | grep duration | cut -d" " -f 2`
  steps=`tr -d '\r' < "$head-settings.txt" | grep steps | cut -d" " -f 2`
  variables=`tr -d '\r' < "$head-settings.txt" | grep variables | cut -d":" -f 2 | sed -e "s/ //g"`
  amount=`tr -d '\r' < "$head-settings.txt" | grep amount | cut -d":" -f 2 | sed -e "s/ //g"`
  concentration=`tr -d '\r' < "$head-settings.txt" | grep amount | cut -d":" -f 2 | sed -e "s/ //g"`
  atol='1e-16'
  rtol='1e-10'
  #sbml="$head-sbml-l2v4.xml"
  sbml="$head-sbml-l3v1.xml"
  result="$i-results.csv"
  opt_delta=""
  opt_amount=""
  opt_facmax=""
  print_msg=""
  dbl_sp="\040\040"
  if [[ "$fine_delta" == *"$i"* ]]; then
    tmp_msg="${dbl_sp}simulate with fine delta\n"
    print_msg="${print_msg}${tmp_msg}"
    opt_delta="-d 0.00001"
  fi
  if [[ "$fine_atol" == *"$i"* ]]; then
    tmp_msg="${dbl_sp}simulate with fine absolute tolerance\n"
    print_msg=${print_msg}${tmp_msg}
    atol="1e-22"
  fi
  if [[ "$fine_rtol" == *"$i"* ]]; then
    print_msg=${print_msg}"${dbl_sp}simulate with fine relative tolerance\n"
    rtol="1e-11"
  fi
  if [[ "$fine_facmax" == *"$i"* ]]; then
    print_msg=${print_msg}"${dbl_sp}simulate with very fine facmax\n"
    opt_facmax="-M 1.0000001"
  fi
  if [ -n "$amount" ] ; then
    print_msg=${print_msg}"${dbl_sp}print amount\n"
	  opt_amount="-a"
  else
    print_msg=${print_msg}"${dbl_sp}print concentration\n"
  fi
  echo "$i: $duration : $steps : [$variables] : [$amount] : [$concentration] : $atol : $rtol"
  echo -en $print_msg
  #./simulateSBML -t $duration -s $steps $opt_delta -m 1 -n $opt_amount $sbml && \
  ./simulateSBML -t $duration -s $steps $opt_delta -m 13 -A $atol -R $rtol $opt_facmax -n $opt_amount $sbml && \
  ./genresult.pl out.csv $variables $steps > $result 
  echo -en $dbl_sp
  ./compare.pl $i

  unset head
  unset duration
  unset steps
  unset variables
  unset sbml
  unset result
  unset opt_delta
  unset opt_amount
  unset opt_facmax
  unset atol
  unset rtol
  unset print_msg
  unset tmp_msg
# end
done
unset BaseDir
unset fine_delta
rm -f out.csv
zip result.zip 00*-results.csv
#rm 00*-results.csv
