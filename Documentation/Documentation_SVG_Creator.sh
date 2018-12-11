#!/bin/bash

echo "|------------------------------------------------------------|"
echo "| Tool to create SVGs used in the documentation              |"
echo "|   along with other included htmls such as colnumbering     |"
echo "| Various transformations are done on the SVGs produced      |"
echo "|   to keep the interactability and selected/rejected colors |"
echo "|------------------------------------------------------------|"

I=../dataset/example.010.AA.fasta

function cleanSVG { 
	sed -i'' "s/#selected/#selected_$1/g" $2
	sed -i'' "s/#rejected/#rejected_$1/g" $2
	sed -i'' "s/id=\"selected/id=\"selected_$1/g" $2
	sed -i'' "s/id=\"rejected/id=\"rejected_$1/g" $2
	sed -i'' "s/Label/Label_$1/g" $2
	sed -i'' "s/Line/Line_$1/g" $2
}

function perform {
	echo "|-- > Performing \"$2\""
	../bin/trimal -out /dev/null -in $I -svgout $2 "${@:3}"
	cleanSVG $1 $2
}

rm -rf ./tmp
mkdir tmp

echo "|-> Creating files for 01_CleaningMethods"

mkdir tmp/01_CleaningMethods/

perform strict tmp/01_CleaningMethods/02_Strict.svg -strict

perform strictplus tmp/01_CleaningMethods/02_StrictPlus.svg -strict

perform gappyout tmp/01_CleaningMethods/02_Gappyout.svg -gappyout

perform gt05 tmp/01_CleaningMethods/03_gt05.svg -gt 0.5

perform st05 tmp/01_CleaningMethods/03_st05.svg -st 0.5

perform st05gt1 tmp/01_CleaningMethods/03_st05gt1.svg -st 0.5 -gt 1

perform s90r01 tmp/01_CleaningMethods/04_s90r01.svg -seqoverlap 90 -resoverlap 0.1

perform s955r0875 tmp/01_CleaningMethods/04_s955r0875.svg -seqoverlap 95.5 -resoverlap 0.875

perform rselect tmp/01_CleaningMethods/05_rselect0-4-6-12-16-18.svg -selectcols { 0,4-6,12,16-18 }

perform sselect tmp/01_CleaningMethods/05_sselect0-2.svg -selectseqs { 0-2 }

perform rsselect tmp/01_CleaningMethods/05_rselect0-4-6-12-16-18_sselect0-2.svg -selectcols { 0,4-6,12,16-18 } -selectseqs { 0-2 }

perform nogaps tmp/01_CleaningMethods/06_nogaps.svg -nogaps

# No all gaps need a special alignment, with allgaps in one column
#   Available alignments don't have this characteristic.
#   We found example 10 to have columns with all gaps but one residue.
#   Removing the sequence with the residue, we get the 'allgaps' column.
I=../dataset/example.010b.AA.fasta
perform noallgaps tmp/01_CleaningMethods/06_noallgaps.svg -noallgaps

# Back to normal
I=../dataset/example.010.AA.fasta

echo "|-> Creating files for 02_Trimming"

mkdir tmp/02_Trimming/

perform nocomplementary tmp/02_Trimming/02_gt05.svg -gt 0.5
perform complementary tmp/02_Trimming/02_gt05C.svg -gt 0.5 -complementary

perform noterminal tmp/02_Trimming/03_noterminal.svg -st 1 -gt 1
perform terminal tmp/02_Trimming/03_terminal.svg -st 1 -gt 1 -terminalonly

perform blo1 tmp/02_Trimming/04_b1.svg -strict -block 1
perform blo2 tmp/02_Trimming/04_b2.svg -strict -block 2
perform blo4 tmp/02_Trimming/04_b4.svg -strict -block 4

perform clus2 tmp/02_Trimming/05_clusters2.svg -clusters 2
perform clus4 tmp/02_Trimming/05_clusters4.svg -clusters 4
perform clus6 tmp/02_Trimming/05_clusters6.svg -clusters 6

perform identity2 tmp/02_Trimming/06_identity2.svg -maxidentity 0.2
perform identity8 tmp/02_Trimming/06_identity8.svg -maxidentity 0.8


# Compareset works differently
    echo "|-- > Performing \"tmp/02_Trimming/07_compare.svg\""
    ../bin/trimal -out /dev/null -compareset ../dataset/alignments_comparison.5 -ct 0.8 -svgout tmp/02_Trimming/07_compare.svg > tmp/02_Trimming/07_compare.html
    cleanSVG compare tmp/02_Trimming/07_compare.svg


perform cons0 tmp/02_Trimming/08_aligconserv0.svg -st 0.99 -cons 0
perform cons40 tmp/02_Trimming/08_aligconserv40.svg -st 0.99 -cons 40
perform cons60 tmp/02_Trimming/08_aligconserv60.svg -st 0.99 -cons 60
perform cons80 tmp/02_Trimming/08_aligconserv80.svg -st 0.99 -cons 80


echo "|-> Creating files for 03_OtherTweaks"

mkdir tmp/03_OtherTweaks/

# Colnumber works differently as we want to keep the colnumbering string
    echo "|-- > Performing \"tmp/03_OtherTweaks/01_colnumbering.svg\""
    ../bin/trimal -out /dev/null -svgout tmp/03_OtherTweaks/01_colnumbering.svg -colnumbering -strict -in $I | tr -d '\012\015' > tmp/03_OtherTweaks/01_colnumbering.html
    cleanSVG compare tmp/02_Trimming/07_compare.svg


perform svgout tmp/03_OtherTweaks/03_SVG_out.svg -htmlout tmp/03_OtherTweaks/03_HTML_out.html -strict