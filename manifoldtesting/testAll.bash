# get the desired input directory from the user
echo "please enter the directory of the input files"
inputDir="../handout_models"
outputFile="manifold_test_results.txt"
echo "please enter the file path to the built program"
programPath="../cmake-build-release/manifoldtesting/manifoldtesting"
# clear the output file
echo -n "" > $outputFile
# for each file in the input directory ending in .tri
for file in $inputDir/*.tri
do
    # get the filename without the path
    filename=$(basename -- "$file")
    # echo the filename to the console
    echo $filename
    # run the program with the file as input and store the output in a variable
    output=$(./$programPath $file 2>&1)
    # write the filename to the output file
      echo -n $filename >> $outputFile
      echo -n " " >> $outputFile
    # for each line in the output
    while IFS= read -r line
    do
      echo $line

      # if the li begins with "Mesh is not manifold"
      if [[ $line == "Mesh is not manifold"* ]]
      then
#        # write the line to the output file
        echo $line >> $outputFile
      fi
      # if the line begins with "The total genus of the mesh is"
      if [[ $line == "The total genus of the mesh is"* ]]
      then
        # get the genus from the line
        genus=$(echo $line | grep -o -E '[0-9]+')
        echo "Genus: $genus" >> $outputFile
      fi
    done <<< "$output"
done