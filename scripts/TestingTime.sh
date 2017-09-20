#!/bin/bash

inputFile=$1
outputFile=$2

echo -e "File $1" > $2

echo -e "\n\n[new_trimAl] vanilla execution" >> $2
(time bin/trimal -in $inputFile -out /tmp/null) &>> $2

echo -e "\n\n[trimAl] vanilla execution" >> $2
(time ../Strimal/source/trimal -in $inputFile -out /tmp/null) &>> $2



echo -e "\n\n[new_trimAl] automated1 execution" >> $2
(time (bin/trimal -in $inputFile -automated1 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2

echo -e "\n\n[trimAl] automated1 execution" >> $2
(time (../Strimal/source/trimal -in $inputFile -automated1 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2



echo -e "\n\n[new_trimAl] gappyout execution" >> $2
(time (bin/trimal -in $inputFile -gappyout -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2

echo -e "\n\n[trimAl] gappyout execution" >> $2
(time (../Strimal/source/trimal -in $inputFile -gappyout -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2



echo -e "\n\n[new_trimAl] strict execution" >> $2
(time (bin/trimal -in $inputFile -strict -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2

echo -e "\n\n[trimAl] strict execution" >> $2
(time (../Strimal/source/trimal -in $inputFile -strict -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2



echo -e "\n\n[new_trimAl] strictplus execution" >> $2
(time (bin/trimal -in $inputFile -strictplus -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2

echo -e "\n\n[trimAl] strictplus execution" >> $2
(time (../Strimal/source/trimal -in $inputFile -strictplus -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2



echo -e "\n\n[new_trimAl] nogaps execution" >> $2
(time (bin/trimal -in $inputFile -nogaps -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2

echo -e "\n\n[trimAl] nogaps execution" >> $2
(time (../Strimal/source/trimal -in $inputFile -nogaps -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2



echo -e "\n\n[new_trimAl] noallgaps execution" >> $2
(time (bin/trimal -in $inputFile -noallgaps -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2

echo -e "\n\n[trimAl] noallgaps execution" >> $2
(time (../Strimal/source/trimal -in $inputFile -noallgaps -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2



echo -e "\n\n[new_trimAl] cluster 4 execution" >> $2
(time (bin/trimal -in $inputFile -cluster 4 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2

echo -e "\n\n[trimAl] cluster 4 execution" >> $2
(time (../Strimal/source/trimal -in $inputFile -cluster 4 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2


echo -e "\n\n[new_trimAl] max identity 0.5 " >> $2
(time (bin/trimal -in $inputFile -maxidentity 0.5 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2

echo -e "\n\n[trimAl] cluster 4 execution" >> $2
(time (../Strimal/source/trimal -in $inputFile -maxidentity 0.5 -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2



echo -e "\n\n[readAlMS] to fasta" >> $2
(time (bin/readalMS -in $inputFile -fasta -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2

echo -e "\n\n[readal] to fasta" >> $2
(time (../Strimal/source/readal -in $inputFile -fasta -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2



echo -e "\n\n[readAlMS] to clustal" >> $2
(time (bin/readalMS -in $inputFile -clustal -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2

echo -e "\n\n[readal] to clustal" >> $2
(time (../Strimal/source/readal -in $inputFile clustal -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2



echo -e "\n\n[readAlMS] to phylip" >> $2
(time (bin/readalMS -in $inputFile -phylip -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2

echo -e "\n\n[readal] to phylip" >> $2
(time (../Strimal/source/readal -in $inputFile -phylip -out /tmp/null > /tmp/null 2>&1 ) ) &>> $2
