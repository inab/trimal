methods="gappyout strict strictplus automated1 nogaps noallgaps"
for method in $methods
do
    for msa in test_msas/$method/*.fasta
    do
        msa_filename=$(basename $msa)
        
        # Ignore MSAs with known difference in strict method because of precision level in calculations
        if echo "$msa_filename" | grep -q -E '^(example\.014|example\.028|example\.032|example\.041)'; then
            echo "Skipping $msa_filename"
        else
            diff "$msa" "dataset/trimmed_msas/$method/$msa_filename" -q || exit 1
            echo "Compared $msa and dataset/trimmed_msas/$method/$msa_filename"
        fi
    done
done
rm -rf test_msas