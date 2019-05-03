/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2009-2019
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@crg.es)

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

***************************************************************************** */

#include "reportsystem.h"

const std::map<WarningCode, const char *> reporting::reportManager::WarningMessages =
{
    {WarningCode::RemovingOnlyGapsSequence,
            "Removing sequence '[tag]' composed only by gaps after trimming"},

    {WarningCode::KeepingOnlyGapsSequence,
            "Keeping sequence '[tag]' composed only by gaps after trimming"},

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
