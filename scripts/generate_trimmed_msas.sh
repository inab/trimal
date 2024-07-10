mkdir test_msas
mkdir test_msas/gappyout
mkdir test_msas/strict
mkdir test_msas/strictplus
mkdir test_msas/automated1
mkdir test_msas/nogaps
mkdir test_msas/noallgaps

./source/trimal -in dataset/example.004.AA.fasta -gappyout > test_msas/gappyout/example.004.AA.fasta
./source/trimal -in dataset/example.004.AA.fasta -strict > test_msas/strict/example.004.AA.fasta
./source/trimal -in dataset/example.004.AA.fasta -strictplus > test_msas/strictplus/example.004.AA.fasta
./source/trimal -in dataset/example.004.AA.fasta -automated1 > test_msas/automated1/example.004.AA.fasta
./source/trimal -in dataset/example.004.AA.fasta -nogaps > test_msas/nogaps/example.004.AA.fasta
./source/trimal -in dataset/example.004.AA.fasta -noallgaps > test_msas/noallgaps/example.004.AA.fasta

