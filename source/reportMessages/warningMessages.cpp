#include "reportsystem.h"

const std::map<WarningCode, const char *> reporting::reportManager::WarningMessages =
{
    {WarningCode::RemovingOnlyGapsSequence,
            "Removing sequence '[tag]' composed only by gaps"},

    {WarningCode::KeepingOnlyGapsSequence,
            "Keeping sequence '[tag]' composed only by gaps"},

    {WarningCode::SequenceWillBeCut,
            "Sequence \"[tag]\" will be cut at position [tag] (length:[tag])"},

    {WarningCode::IncludingIndeterminationSymbols,
            "Sequence \"[tag]\" has some indetermination symbols 'X' "
            "at the end of sequence. They will be included in the final Alignment."},

    {WarningCode::LessNucleotidesThanExpected,
            "Sequence \"[tag]\" has less nucleotides ([tag]) than expected ([tag])."
            " It will be added N's to complete the sequence"},

    {WarningCode::HeaderWillBeCut,
            "Original sequence header will be cut by 10 characters on format [tag]"},

    {WarningCode::DonorAlreadyAdded,
            "The donor \"[tag]\" is present on more than one VCF. "
            "Overlaping SNPs will be overwritten."},

    {WarningCode::SNPAlreadApplied,
            "SNP already applied to \"[tag]\":\"[tag]\" at position [tag] \"[tag]\"->\"[tag]\""},

    {WarningCode::OverwrittingFile,
            "[[tag]] Overwritting file [tag]."},

    {WarningCode::RenamingOutputPreventOverride,
            "[[tag]] -> To prevent overriding file [tag] a suffix has been added. Final filename: [tag]"}
};
