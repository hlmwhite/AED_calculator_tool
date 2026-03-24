#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.annotation_gff = ""
params.evidence_gff = ""


process GFFCOMPARE {
    input:
    path ann_gff
    path ev_gff

    output:
    path 'gffcmp.loci.test'

    shell:
    '''
    gffcompare -r !{ann_gff} !{ev_gff}
    head -n 40 gffcmp.loci > gffcmp.loci.test
    '''
}

process CALCULATE_AED {
    input:
    val gffcmp_string
    path ann_gff
    path ev_gff

    output:
    stdout

    shell:
    '''
    bash /Users/mark/Documents/work_stuff/python_course/aed_calc_bash.sh "!{gffcmp_string}" !{ann_gff} !{ev_gff} > aed_results.txt 2>&1
    '''
}

//params.input_list = 'files.txt'


workflow {
    annGff = file(params.annotation_gff)
    evGff = file(params.evidence_gff)

    input_list = GFFCOMPARE(annGff, evGff)

    //files_ch = Channel

    //input_list.splitText().map{it.trim()}

    CALCULATE_AED(input_list.splitText().map{it.trim()}, annGff, evGff )
}


    //annIDs=$(printf "!{gffcmp_string}" | awk '{print $4}' | tr ',' '\n')
    //evIDs=$(printf "!{gffcmp_string}" | awk '{print $5}'| tr ',' '\n')
    //for annID in $annIDs
    //do
    //    for evID in $evIDs
    //    do
    //        echo "Calculating AED for annotation ID: $annID and evidence ID: $evID"
    //        python3 example_calculator_v2.py -A !{ann_gff} -E !{ev_gff} -a $annID -e $evID > aed_results.txt 2>&1
    //    done
    //done 