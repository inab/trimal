#!/bin/bash

inputFolder=$1
outputFile=$2

echo -e "Preparing the MakeFile"
(cmake -DCMAKE_BUILD_TYPE=RELEASE .) >> /dev/null
echo -e "Compiling using MakeFile"
(make -j8) >> /dev/null  
echo -e "Copying executable"
(cp bin/trimal bin/trimalStable)
echo -e "Starting the timing tests"

# thresholds=(0.5)
# cons=(10)
automethod=('-fasta' '-gappyout' '-strict' '-strictplus' '-automated1' '-nogaps' '-noallgaps')
# # resoverlap=(0.5)
# # seqoverlap=(5)
# clusters=(2)
# # max_identities=(0.5)
# # blocks=(2)
# # formats=(-nbrf -mega -nexus -clustal -fasta -fasta_m10 -phylip -phylip_m10 -phylip_paml -phylip_paml_m10 -phylip3.2 -phylip3.2_m10)
# # external=(-sgc -sgt -ssc -sst -sfc -sft -sident -soverlap)
files=($(ls -Sr "$1"))

for method in ${automethod[@]}
do
    echo -e "Method ${method}" >> ${outputFile}
    sleep 1
    for file in ${files[@]} 
    do
        echo -e "\nFile ${file}" >> ${outputFile}
        stat --printf="%s" ${inputFolder}/${file} >> ${outputFile}
        echo -e "\n\t[new_trimAl]" >> ${outputFile}
        (time bin/trimalStable -in ${inputFolder}/${file} -out /tmp/null) &>> ${outputFile}

        echo -e "\n\t[trimAl]" >> ${outputFile}
        (time ../Strimal/source/trimal -in ${inputFolder}/${file} -out /tmp/null) &>> ${outputFile}
    done
done

rm bin/trimalStable





# for file in ${files[*]}
# do
# echo -e $file
# done
# 
# echo -e "File $1" > $2

# echo -e "\n\n[new_trimAl] vanilla execution" >> $2
# (time bin/trimal -in $inputFile -out /tmp/null) &>> $2
# 
# echo -e "\n\n[trimAl] vanilla execution" >> $2
# (time ../Strimal/source/trimal -in $inputFile -out /tmp/null) &>> $2
# 
# 
# 
# echo -e "\n\n[new_trimAl] automated1 execution" >> $2
# (time (bin/trimal -in $inputFile -automated1 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# echo -e "\n\n[trimAl] automated1 execution" >> $2
# (time (../Strimal/source/trimal -in $inputFile -automated1 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# 
# 
# echo -e "\n\n[new_trimAl] gappyout execution" >> $2
# (time (bin/trimal -in $inputFile -gappyout -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# echo -e "\n\n[trimAl] gappyout execution" >> $2
# (time (../Strimal/source/trimal -in $inputFile -gappyout -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# 
# 
# echo -e "\n\n[new_trimAl] strict execution" >> $2
# (time (bin/trimal -in $inputFile -strict -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# echo -e "\n\n[trimAl] strict execution" >> $2
# (time (../Strimal/source/trimal -in $inputFile -strict -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# 
# 
# echo -e "\n\n[new_trimAl] strictplus execution" >> $2
# (time (bin/trimal -in $inputFile -strictplus -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# echo -e "\n\n[trimAl] strictplus execution" >> $2
# (time (../Strimal/source/trimal -in $inputFile -strictplus -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# 
# 
# echo -e "\n\n[new_trimAl] nogaps execution" >> $2
# (time (bin/trimal -in $inputFile -nogaps -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# echo -e "\n\n[trimAl] nogaps execution" >> $2
# (time (../Strimal/source/trimal -in $inputFile -nogaps -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# 
# 
# echo -e "\n\n[new_trimAl] noallgaps execution" >> $2
# (time (bin/trimal -in $inputFile -noallgaps -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# echo -e "\n\n[trimAl] noallgaps execution" >> $2
# (time (../Strimal/source/trimal -in $inputFile -noallgaps -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# 
# 
# echo -e "\n\n[new_trimAl] cluster 4 execution" >> $2
# (time (bin/trimal -in $inputFile -cluster 4 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# echo -e "\n\n[trimAl] cluster 4 execution" >> $2
# (time (../Strimal/source/trimal -in $inputFile -cluster 4 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# 
# echo -e "\n\n[new_trimAl] max identity 0.5 " >> $2
# (time (bin/trimal -in $inputFile -maxidentity 0.5 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# echo -e "\n\n[trimAl] cluster 4 execution" >> $2
# (time (../Strimal/source/trimal -in $inputFile -maxidentity 0.5 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# 
# 
# echo -e "\n\n[readAlMS] to fasta" >> $2
# (time (bin/readalMS -in $inputFile -fasta -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# echo -e "\n\n[readal] to fasta" >> $2
# (time (../Strimal/source/readal -in $inputFile -fasta -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# 
# 
# echo -e "\n\n[readAlMS] to clustal" >> $2
# (time (bin/readalMS -in $inputFile -clustal -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# echo -e "\n\n[readal] to clustal" >> $2
# (time (../Strimal/source/readal -in $inputFile clustal -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# 
# 
# echo -e "\n\n[readAlMS] to phylip" >> $2
# (time (bin/readalMS -in $inputFile -phylip -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
# 
# echo -e "\n\n[readal] to phylip" >> $2
# (time (../Strimal/source/readal -in $inputFile -phylip -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
