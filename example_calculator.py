import gffpandas.gffpandas as gffpd

annotation = gffpd.read_gff3("gene1_ann.gff")
evidence = gffpd.read_gff3("gene2_ev.gff")

# annotation_file = gffpd.read_gff3("copy.genomic.gff")
# evidence_file = gffpd.read_gff3("out.stringtie.gtf")

# gffcmp_loci_file = pd.read_csv("gffcmp.loci", sep="\t", header=None)

ann_exon_df = annotation.filter_feature_of_type(["exon"])
intervals_ann = ann_exon_df.df[["start", "end"]].values.tolist()

ev_exon_df = evidence.filter_feature_of_type(["exon"])
intervals_ev = ev_exon_df.df[["start", "end"]].values.tolist()

# overlap_count = 0
# non_overlap_count = 0


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


tp, fn = get_overlaps_sen(intervals_ann, intervals_ev, 0, 0)

print("Number of overlapping values (true positives):", tp)
print("Number of non-overlapping values (false negatives):", fn)


def calculate_sensitivity(tp, fn):
    if (tp + fn) == 0:
        return 0.0
    sensitivity = tp / (tp + fn)
    return sensitivity


SN = calculate_sensitivity(tp, fn)
print("SN (or sensitivity):", SN)

# SP can be thought of as the fraction of i overlapping j. where i is the evidence and j is the annotation/predicition

# overlap_count = 0
total_values = sum(end - start + 1 for start, end in intervals_ev)


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


overlap_count = get_overlaps_spec(intervals_ann, intervals_ev, 0)

print("Number of overlapping values:", overlap_count)

# total_values = sum(end - start + 1 for start, end in intervals_ann)

# SP can be thought of as the fraction of i overlapping j.


def calculate_specificity(tot, ovlCount):
    specificity = ovlCount / tot if tot > 0 else 0
    return specificity


SP = calculate_specificity(total_values, overlap_count)
print("SP (or specificity):", SP)


def calculate_accuracy(sn, sp):
    accuracy = (sn + sp) / 2
    return accuracy


def calculate_distance(acc):
    distance = 1 - acc
    return distance


AC = calculate_accuracy(SN, SP)

print("distance:", calculate_distance(AC))
