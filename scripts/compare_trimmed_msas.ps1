# Define methods
$methods = "gappyout", "strict", "strictplus", "automated1", "nogaps", "noallgaps"

# Compare files
foreach ($method in $methods) {
    $files = Get-ChildItem -Path "test_msas\$method" -Filter *.fasta
    foreach ($file in $files) {
        $msa_filename = $file.Name
        $test_file = $file.FullName
        $reference_file = Join-Path -Path "dataset\trimmed_msas\$method" -ChildPath $msa_filename

        if (-Not (Test-Path $reference_file)) {
            Write-Output "Reference file $reference_file does not exist."
            exit 1
        }

        if (-Not (Test-Path $test_file)) {
            Write-Output "Test file $test_file does not exist."
            exit 1
        }

        $test_content = Get-Content $test_file
        $reference_content = Get-Content $reference_file

        # Check if both files are empty
        if ($test_content.Count -eq 0 -and $reference_content.Count -eq 0) {
            Write-Output "Both $test_file and $reference_file are empty and considered equal."
            continue
        }

        # Compare contents
        $result = Compare-Object $test_content $reference_content -SyncWindow 0
        if ($result) {
            Write-Output "Files $test_file and $reference_file differ."
            exit 1
        }
        else {
            Write-Output "Compared $test_file and $reference_file"
        }
    }
}

# Remove test_msas directory
Remove-Item -Recurse -Force test_msas
