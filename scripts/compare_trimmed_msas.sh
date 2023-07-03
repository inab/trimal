methods=(gappyout strict strictplus automated1 nogaps noallgaps)
for method in $methods
do
    for msa in test_msas/$method/*.fasta
    do
        msa_filename=$(basename $msa)
        diff $msa dataset/trimmed_msas/$method/$msa_filename -q || exit 1
    done
done
rm -rf test_msas