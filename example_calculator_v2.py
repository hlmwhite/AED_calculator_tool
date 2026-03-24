import argparse
import re

import gffpandas.gffpandas as gffpd


def extract_exon_ints(
    df, tx_id, fType, attr_col="attributes"
):  # extract exon hints into a list of intervals for a given transcript ID
    """
    Return rows for the transcript feature (ID=tx_id) and all child features (Parent=tx_id).
    df: a pandas DataFrame with an attributes column as a string.
    tx_id: transcript ID string, e.g. "transcript:XM_12345" or "tx1"
    """
    # match whole attribute fields to avoid partial matches
    if fType == "ANN":
        id_pat = re.compile(rf"(^|;)ID={re.escape(tx_id)}(;|$)")
        parent_pat = re.compile(rf"(^|;)Parent={re.escape(tx_id)}(;|$)")
        attrs = df[attr_col].astype(str)
        mask = attrs.str.contains(id_pat) | attrs.str.contains(parent_pat)
    elif fType == "STRINGTIE":
        tx_pat = re.compile(rf'(^|;\s*)transcript_id\s+"{re.escape(tx_id)}"\s*;')
        attrs = df[attr_col].astype(str)
        mask = attrs.str.contains(tx_pat)

    # return df.loc[mask].copy()
    subset_df = df.loc[mask].copy()
    exon_df = subset_df[subset_df["type"] == "exon"].copy()
    # exon_df = mask.filter_feature_of_type(["exon"])
    exon_df["start"] = exon_df["start"].astype(int)
    exon_df["end"] = exon_df["end"].astype(int)
    # intervals = exon_df[["start", "end"]].values.tolist()
    intervals = list(zip(exon_df["start"], exon_df["end"]))
    return intervals


def get_overlaps_sen(intervals_ann, intervals_ev, overlap_count, non_overlap_count):
    for int2 in intervals_ann:
        start2, end2 = int2
        # Loop through every integer value in the range int2 (inclusive)
        for value in range(start2, end2 + 1):
            found_overlap = False
            # Check if the current value falls within any interval in intervals1
            for int1 in intervals_ev:
                start1, end1 = int1
                if start1 <= value <= end1:  # inclusive boundaries
                    found_overlap = True
                    break
            if found_overlap:
                overlap_count += 1
            else:
                non_overlap_count += 1
    tp_ovlCount = overlap_count
    fn_nonOvlCount = non_overlap_count
    return tp_ovlCount, fn_nonOvlCount


def calculate_sensitivity(tp, fn):
    if (tp + fn) == 0:
        return 0.0
    sensitivity = tp / (tp + fn)
    return sensitivity


def get_overlaps_spec(intervals_ann, intervals_ev, overlap_count):
    for int2 in intervals_ev:
        start2, end2 = int2
        # Loop through every integer value in the range int2 (inclusive)
        for value in range(start2, end2 + 1):
            found_overlap = False
            # Check if the current value falls within any interval in intervals1
            for int1 in intervals_ann:
                start1, end1 = int1
                if start1 <= value <= end1:  # inclusive boundaries
                    found_overlap = True
                    break
            if found_overlap:
                overlap_count += 1
    return overlap_count


def calculate_specificity(tot, ovlCount):
    specificity = ovlCount / tot if tot > 0 else 0
    return specificity


def calculate_accuracy(sn, sp):
    accuracy = (sn + sp) / 2
    return accuracy


def calculate_distance(acc):
    distance = 1 - acc
    return distance


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="AED calculator with parallel processing",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("-A", "--annotation", required=True, help="Gene annotation/prediction in GFF format")
    parser.add_argument("-E", "--evidence", required=True, help="Evidence GFF file(s)")
    parser.add_argument(
        "-a",
        "--annotation_id",
        required=True,
        help="ID of the transcript in the annotation GFF to compare",
        default=None,
    )
    parser.add_argument(
        "-e",
        "--evidence_id",
        required=True,
        help="ID of the transcript in the evidence GFF to compare",
        default=None,
    )

    args = parser.parse_args()

    annotationFile = gffpd.read_gff3(args.annotation)
    evidenceFile = gffpd.read_gff3(args.evidence)  # Assuming only one evidence file for now

    intervals_ann = extract_exon_ints(annotationFile.df, args.annotation_id, "ANN")
    intervals_ev = extract_exon_ints(evidenceFile.df, args.evidence_id, "STRINGTIE")

    tp, fn = get_overlaps_sen(intervals_ann, intervals_ev, 0, 0)

    print("Number of overlapping values (true positives):", tp)
    print("Number of non-overlapping values (false negatives):", fn)

    SN = calculate_sensitivity(tp, fn)
    print("SN (or sensitivity):", SN)

    # SP can be thought of as the fraction of i overlapping j. where i is the evidence and j is the annotation/predicition

    # overlap_count = 0
    total_values = sum(end - start + 1 for start, end in intervals_ev)

    overlap_count = get_overlaps_spec(intervals_ann, intervals_ev, 0)

    print("Number of overlapping values:", overlap_count)

    # total_values = sum(end - start + 1 for start, end in intervals_ann)

    # SP can be thought of as the fraction of i overlapping j.

    SP = calculate_specificity(total_values, overlap_count)
    print("SP (or specificity):", SP)

    AC = calculate_accuracy(SN, SP)
    print("distance:", calculate_distance(AC))
