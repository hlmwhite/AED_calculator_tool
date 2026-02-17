from example_calculator_v2 import (calculate_accuracy, calculate_distance,
                                   calculate_sensitivity,
                                   calculate_specificity, extract_exon_ints)


def test_sensitivity_basic():
    result = calculate_sensitivity(3, 2)
    assert result == 0.6


def test_specificity_basic_1():
    result = calculate_specificity(10, 5)
    assert result == 0.5


def test_specificity_basic_2():
    result = calculate_specificity(10, 0)
    assert result == 0.0


def test_calculate_distance():
    acc = 0.75
    result = calculate_distance(acc)
    assert result == 0.25


def test_calculate_accuracy():
    sn = 0.8
    sp = 0.6
    result = calculate_accuracy(sn, sp)
    assert result == 0.7


def test_extract_exon_ints():
    gff_file = "gene2_ev.gff"
    tx_id = "evm.model.1ctg.1717"
    exon_intervals = extract_exon_ints(gff_file, tx_id, "ANN")
    expected_intervals = [
        (17210569, 17210871),
        (17204220, 17204546),
        (17202665, 17202808),
        (17200809, 17202570),
        (17199444, 17200615),
        (17198724, 17199371),
        (17195834, 17198463),
    ]
    assert exon_intervals == expected_intervals


def test_placeholder():
    assert True


# This is a placeholder test file. Replace with actual tests.
