#!/bin/bash

### directories and files
# ~/Downloads/SBML-testcases/cases/semantic/00980/
#   00980-settings.txt
#   00980-sbml-l2v4.xml
#   00980-results.csv

### parse args
while [[ $# -gt 1 ]]
do
  test_case_dir="$1"
  shift
  num="$1"
  shift
  out_dir="$1"
  shift
  level="$1"
  shift
  version="$1"
  shift
done
basedir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
sim="$basedir/simulateSBML"
genresult="$basedir/genresult.pl"

### 00001-settings.txt
# duration: 5
# steps: 50
# variables: S1, S2, S3, S4
# absolute: 1.000000e-007
# relative: 0.0001
# amount: S1, S2
# concentration:

# Simulate following models with small dt
fine_delta="00952 00953 00962 00963 00964 00965 00966 00967 00986 00987 00988 01000 01113 01121 01123"

# Simulate following models with small atol
fine_atol="00430 00431 00863 00893 00945 00946 00947 00948 01000 01104 01121 01123"

# Simulate following models with small rtol
fine_rtol="00944 00945 00946 00947 00948 00989 01000 01104 01121 01123"

# Simulate following models with small facmax
fine_facmax="00426 00430 00431 00863 01000 01104 01121 01123"

head="$test_case_dir/$num/$num"
#sbml="$head-sbml-l2v4.xml"
sbml="${head}-sbml-l${level}v${version}.xml"
if [[ ! -f "$sbml" ]]; then
  echo $sbml not found.
  continue
fi
# 00926 and 00927 contains '\r' before '\n' on each line (DOS format),
# so we have to call tr -d '\r' before parsing the settings file...
duration=`tr -d '\r' < "$head-settings.txt" | grep duration | cut -d" " -f 2`
steps=`tr -d '\r' < "$head-settings.txt" | grep steps | cut -d" " -f 2`
variables=`tr -d '\r' < "$head-settings.txt" | grep variables | cut -d":" -f 2 | sed -e "s/ //g"`
amount=`tr -d '\r' < "$head-settings.txt" | grep amount | cut -d":" -f 2 | sed -e "s/ //g"`
concentration=`tr -d '\r' < "$head-settings.txt" | grep amount | cut -d":" -f 2 | sed -e "s/ //g"`
atol='1e-16'
rtol='1e-10'
result="$num.csv"
opt_delta=""
opt_amount=""
opt_facmax=""
print_msg=""
dbl_sp="\040\040"
if [[ "$fine_delta" == *"$num"* ]]; then
  tmp_msg="${dbl_sp}simulate with fine delta\n"
  print_msg="${print_msg}${tmp_msg}"
  opt_delta="-d 0.00001"
fi
if [[ "$fine_atol" == *"$num"* ]]; then
  tmp_msg="${dbl_sp}simulate with fine absolute tolerance\n"
  print_msg=${print_msg}${tmp_msg}
  atol="1e-22"
fi
if [[ "$fine_rtol" == *"$num"* ]]; then
  print_msg=${print_msg}"${dbl_sp}simulate with fine relative tolerance\n"
  rtol="1e-11"
fi
if [[ "$fine_facmax" == *"$num"* ]]; then
  print_msg=${print_msg}"${dbl_sp}simulate with very fine facmax\n"
  opt_facmax="-M 1.0000001"
fi
if [ -n "$amount" ] ; then
  print_msg=${print_msg}"${dbl_sp}print amount\n"
  opt_amount="-a"
else
  print_msg=${print_msg}"${dbl_sp}print concentration\n"
fi
echo "$num: $duration : $steps : [$variables] : [$amount] : [$concentration] : $atol : $rtol"
echo -en $print_msg
#$sim -t $duration -s $steps $opt_delta -m 13 -A $atol -R $rtol $opt_facmax -n $opt_amount -o $result $sbml && \
$sim -t $duration -s $steps $opt_delta -m 1 -n $opt_amount -o $result $sbml && \
  $genresult $result $variables $steps > $out_dir/$result
rm -f $result
