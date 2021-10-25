!/bin/bash

#This script copies a list of temperatures (2nd column) in a numerated list (1st column) from a FILE_LIST
# having that information and creates new copies of gromacs input file.
#Each new copy differs in its temperature and name (listed in FILE_LIST) from the original one (FILE_INPUT).


FILE_LIST=$1
FILE_INPUT=$2
MDP="$(echo "$FILE_INPUT" | cut -d'.' -f1)"

while IFS= read -r line
do
  num=$(echo "$line" | awk '{print $1}')
  temp=$(echo "$line" | awk '{print $2}')
  cp $FILE_INPUT ${MDP}_${num}.mdp
  sed -i "s/ref_t                   = 300/ref_t                   = $temp/g" ${MDP}_${num}.mdp
  sed -i "s/gen_temp                = 300/gen_temp                = $temp/g" ${MDP}_${num}.mdp
done < $FILE_LIST

exit
