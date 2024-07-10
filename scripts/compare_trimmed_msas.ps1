# Define methods
$methods = "gappyout", "strict", "strictplus", "automated1", "nogaps", "noallgaps"

# Compare files
foreach ($method in $methods) {
    $files = Get-ChildItem -Path "test_msas\$method" -Filter *.fasta
    foreach ($file in $files) {
        $msa_filename = $file.Name
        $result = Compare-Object (Get-Content "$file") (Get-Content "dataset\trimmed_msas\$method\$msa_filename") -SyncWindow 0
        if ($result) {
            Write-Output "Files $file and dataset\trimmed_msas\$method\$msa_filename differ."
            exit 1
        }
        else {
            Write-Output "Compared $file and dataset\trimmed_msas\$method\$msa_filename"
        }
    }
}

# Remove test_msas directory
Remove-Item -Recurse -Force test_msas
