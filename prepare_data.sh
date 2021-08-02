## declare an array variable
#declare -a array=(amazon youtube eu2005 livejournal orkut patents uk2002 friendster twitter)
#declare -a skip_char=(¨#¨ "#" "#" "#" "#" "#" "#" "#" "%")

declare -a array=(amazon)
declare -a skip_char=("#")

dataset_path=sample_dataset

# get length of an array
arraylength=${#array[@]}

# use for loop to read all values and indexes
for (( i=1; i<${arraylength}+1; i++ ));
do
  echo ${array[$i-1]}
	./build/toolset/EdgeList2CSR.out $dataset_path/${array[$i-1]} $dataset_path/${array[$i-1]} ${skip_char[$i-1]}
  ./build/toolset/AssignWeight.out $dataset_path/${array[$i-1]} 1.0 5.0
  ./build/toolset/AssignEdgeLabel.out $dataset_path/${array[$i-1]} 5
  ./build/random_walk/alias_preprocess.out -f $dataset_path/${array[$i-1]}
  ./build/random_walk/its_preprocess.out -f $dataset_path/${array[$i-1]}
  ./build/random_walk/rj_preprocess.out -f $dataset_path/${array[$i-1]}
  echo "Done " $i " / " ${arraylength} " : " ${array[$i-1]}
done
