# Copyright (c) 2025 HaploExplore

from typing import Literal, List, Tuple, Optional
from genetic_structures import SNPInformations, Haploblock
import utils
import time
import psutil
import os
import logging
from multiprocessing import Manager
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import gc
import threading
from queue import Queue
import math

script_dir = os.path.dirname(os.path.realpath(__file__))
    
def calculate_ld_between_two_snps(snp_a : SNPInformations, snp_b : SNPInformations) -> tuple[float, float, float] :
    """
    Calculates r^2, D and D' for two SNPs.
    The SNPInformations have to come from the same vcf file and have the same count of genotypes to properly compute D.
    """
    # Make sure that data is phased for both SNPs a and b.
    if not snp_a.is_phased() and snp_b.is_phased() :
        raise ValueError(f"Cannot compute LD for unphased SNPs : {snp_a.get_snp().get_rs_id()} and {snp_b.get_snp().get_rs_id()} are both unphased")
    if not snp_a.is_phased() :
        raise ValueError(f"Cannot compute LD for unphased SNPs : {snp_a.get_snp().get_rs_id()} is unphased")
    if not snp_b.is_phased() :
        raise ValueError(f"Cannot compute LD for unphased SNPs : {snp_b.get_snp().get_rs_id()} is unphased")

    p_a, p_b = snp_a.get_maf(), snp_b.get_maf()
    q_a, q_b = 1 - p_a, 1 - p_b

    individuals_a, individuals_b = snp_a.get_individuals_genotyped(), snp_b.get_individuals_genotyped()

    if individuals_a.sort() != individuals_b.sort() :
        logging.warning(f"Warning : SNPs {snp_a.get_snp().get_rs_id()} and {snp_b.get_snp().get_rs_id()} don't have the same genotyped individuals : this can lead to LD miscalculation")

    # Count individuals that have both A and B alleles on chromosome #1.
    carriers_allele1_a, carriers_allele1_b = snp_a.get_alt_allele_carriers_chr_1(), snp_b.get_alt_allele_carriers_chr_1()
    common_carriers_allele1 = set(carriers_allele1_a).intersection(carriers_allele1_b)

    # Count individuals that have both A and B alleles on chromosome #2.
    carriers_allele2_a, carriers_allele2_b = snp_a.get_alt_allele_carriers_chr_2(), snp_b.get_alt_allele_carriers_chr_2()
    common_carriers_allele2 = set(carriers_allele2_a).intersection(carriers_allele2_b)

    # Calculate frequency of {A,B} genotype.
    common_carriers_count = len(common_carriers_allele1) + len(common_carriers_allele2)
    p_ab = common_carriers_count/(2*min(snp_a.get_individuals_count(), snp_b.get_individuals_count()))

    # Compute D.
    d = p_ab - p_a * p_b
    # Compute maximum value for |D|.
    if d <= 0 :
        d_max = min(p_a*p_b, q_a*q_b)
    else :
        d_max = min(p_a*q_b, q_a*p_b)

    # Compute D'.
    if d_max == 0:
        d_prime = 0 
    else:
        d_prime = d / d_max

    # Compute r^2;
    denominator = p_a * q_a * p_b * q_b
    if denominator == 0:
        r_square = 0  
    else:
        r_square = d**2 / denominator

    return r_square, d, d_prime

def calculate_percentage_carrier(ref_snp: SNPInformations, snp_test: SNPInformations) -> tuple[float, float, float, float]:
    """
    For a given reference SNP and SNP to test, returns the percentage of individuals carrying 
    the reference allele who also carry the tested allele.
    This handles both phased and unphased SNPs.
    """
    
    # If both SNPs are phased, use chromosome-specific calculations
    if ref_snp.is_phased() and snp_test.is_phased():
        # Compute the percentage for chromosome #1
        ref_snp_carriers_1 = ref_snp.get_alt_allele_carriers_chr_1()
        snp_test_carriers_1 = snp_test.get_alt_allele_carriers_chr_1()
        if ref_snp_carriers_1 is not None and snp_test_carriers_1 is not None:
            common_carriers_1 = set(ref_snp_carriers_1).intersection(snp_test_carriers_1)
            percentage_1 = (len(common_carriers_1) / len(ref_snp_carriers_1)) * 100 if ref_snp_carriers_1 else 0
        else:
            common_carriers_1 = []
            percentage_1 = 0

        # Compute the percentage for chromosome #2
        ref_snp_carriers_2 = ref_snp.get_alt_allele_carriers_chr_2()
        snp_test_carriers_2 = snp_test.get_alt_allele_carriers_chr_2()
        if ref_snp_carriers_2 is not None and snp_test_carriers_2 is not None:
            common_carriers_2 = set(ref_snp_carriers_2).intersection(snp_test_carriers_2)
            percentage_2 = (len(common_carriers_2) / len(ref_snp_carriers_2)) * 100 if ref_snp_carriers_2 else 0
        else:
            common_carriers_2 = []
            percentage_2 = 0

        # Compute the exact percentage including both chromosomes
        all_common_carriers = len(common_carriers_1) + len(common_carriers_2)
        all_ref_snp_carriers = len(ref_snp_carriers_1) + len(ref_snp_carriers_2)
        exact_percentage = (all_common_carriers / all_ref_snp_carriers) * 100 if all_ref_snp_carriers else 0

        # Compute the rough percentage for "phenotype"
        ref_snp_carriers = set(ref_snp_carriers_1).union(ref_snp_carriers_2)
        snp_test_carriers = set(snp_test_carriers_1).union(snp_test_carriers_2)
        common_carriers = ref_snp_carriers.intersection(snp_test_carriers)
        rough_percentage = (len(common_carriers) / len(ref_snp_carriers)) * 100 if ref_snp_carriers else 0

    # If either SNP is unphased, only use the combined carriers
    else:
        # Compute the rough percentage for "phenotype" using all individuals
        ref_snp_carriers = set(ref_snp.get_alt_allele_carriers())
        snp_test_carriers = set(snp_test.get_alt_allele_carriers())
        common_carriers = ref_snp_carriers.intersection(snp_test_carriers)
        rough_percentage = (len(common_carriers) / len(ref_snp_carriers)) * 100 if ref_snp_carriers else 0

        # Since there's no phasing, set exact, percentage_1, and percentage_2 to None
        exact_percentage = rough_percentage
        percentage_1, percentage_2 = None, None

    return rough_percentage, exact_percentage, percentage_1, percentage_2

def extract_snps_carried_by_reference(snp : SNPInformations, snp_region : list[SNPInformations], path : str, cp_cut : float, cp_mode : Literal['rough', 'exact']) -> tuple[list[SNPInformations], str] :
    """
    For given reference SNPInformations and a given list of SNPInformations to browse,
    returns a list of the SNPs that are carried by the reference up to a certain percentage, and the path to a textfile containing the info.
    """
    path_to_list_of_snps_carriers = f"{path}/Carriers.%carrier_{cp_cut}.txt"

    header = f"ref_snp_rsID\tref_snp_position\tsnp_tested_rsID\tsnp_tested_position\tcarrier_percentage_total_{cp_mode}\tcarrier_percentage_chromosome_1\tcarrier_percentage_chromosome_2"

    ref_snp_id = snp.get_snp().get_rs_id()
    ref_snp_pos = snp.get_snp().get_position()

    snps_carriers = []

    with open(path_to_list_of_snps_carriers, 'w', encoding='utf-8') as output_file :
        output_file.write(header+'\n')

        for snp_to_test in snp_region :
            rough_cp, exact_cp, cp_1, cp_2 = calculate_percentage_carrier(snp, snp_to_test)
            if cp_mode == 'rough' and rough_cp < cp_cut :
                continue
            if cp_mode == 'rough' :
                cp = rough_cp
            if cp_mode == 'exact' and exact_cp < cp_cut :
                continue
            if cp_mode == 'exact' :
                cp = exact_cp
            snps_carriers.append(snp_to_test)
            snp_test_id = snp_to_test.get_snp().get_rs_id()
            snp_test_pos = snp_to_test.get_snp().get_position()
            line = f"{ref_snp_id}\t{ref_snp_pos}\t{snp_test_id}\t{snp_test_pos}\t{cp}\t{cp_1}\t{cp_2}"
            output_file.write(line+'\n')

    return snps_carriers, path_to_list_of_snps_carriers

def split_genomic_regions(snp_data: list[SNPInformations], REGION_SIZE: int, MAX_BP_SIZE_FOR_HAPLOBLOCK: int) -> list[list[SNPInformations]]:
    """
    Split the genomic regions into parts of REGION_SIZE with an overlap of MAX_BP_SIZE_FOR_HAPLOBLOCK.
    """
    if snp_data and all(snp is not None for snp in snp_data[:5]):  # Check if snp_data is not empty and the first 5 elements are not None
        for i, snp_info in enumerate(snp_data[:5], start=1):
            logging.info(f"{i}th SNP in data : rs_ID: {snp_info.get_snp().get_rs_id()}, rs_POS: {snp_info.get_snp().get_position()}")
    else:
        logging.warning("snp_data is empty or one of the first 5 elements is None")

    if isinstance(MAX_BP_SIZE_FOR_HAPLOBLOCK, str):
        MAX_BP_SIZE_FOR_HAPLOBLOCK = int(MAX_BP_SIZE_FOR_HAPLOBLOCK)
        logging.info(f"MAX_BP_SIZE_FOR_HAPLOBLOCK type was incorrect, set to: {type(MAX_BP_SIZE_FOR_HAPLOBLOCK)}")
    if isinstance(REGION_SIZE, str):
        REGION_SIZE = int(REGION_SIZE)
        logging.info(f"REGION_SIZE type was incorrect, set to: {type(REGION_SIZE)}")


    # Get the minimum and maximum positions
    min_pos = min(snp_info.get_snp().get_position() for snp_info in snp_data)
    max_pos = max(snp_info.get_snp().get_position() for snp_info in snp_data)

    genomic_regions = []
    start = min_pos
    iteration = 0

    # Loop until the end of the chromosome is reached
    while start < max_pos:
        # Calculate the end position of the current region
        end = start + REGION_SIZE
        #logging.info(f"Iteration {iteration}, Start: {int(start)}, End: {int(end)}")
        # If the end position is beyond the end of the chromosome, set it to the end of the chromosome
        if end > max_pos:
            end = max_pos

        # Get the SNPs that fit in the current region
        region_snps = [snp_info for snp_info in snp_data if int(start) <= snp_info.get_snp().get_position() <= int(end)]

        # Add the current region to the list
        genomic_regions.append(region_snps)

        # Move the start position forward by REGION_SIZE - MAX_BP_SIZE_FOR_HAPLOBLOCK
        start += REGION_SIZE - MAX_BP_SIZE_FOR_HAPLOBLOCK
        iteration += 1

    return genomic_regions

def try_add_snp_to_haploblock_withLD(haploblock, snp, bloc_start, bloc_end, r_square_cut, d_prime_cut, maf_threshold, MAX_BP_SIZE_FOR_HAPLOBLOCK, cp_cut, cp_mode, haploblock_bounds):
    snp_pos = snp.get_snp().get_position()
    core_snp = haploblock.get_core_snp()

    if snp_pos < bloc_start and (bloc_start - snp_pos > MAX_BP_SIZE_FOR_HAPLOBLOCK):
        return False
    if snp_pos > bloc_end and (snp_pos - bloc_end > MAX_BP_SIZE_FOR_HAPLOBLOCK):
        return False

    maf_cut = maf_threshold * core_snp.get_maf() if maf_threshold else None
    is_in_high_ld, _, _, _, _ = caracterize_ld_with_reference(core_snp, snp, r_square_cut, d_prime_cut, maf_cut, cp_cut, cp_mode)

    if not is_in_high_ld:
        return False

    new_start = min(bloc_start, snp_pos)
    new_end = max(bloc_end, snp_pos)
    if (new_end - new_start) <= MAX_BP_SIZE_FOR_HAPLOBLOCK:
        haploblock.get_all_snps().append(snp)
        haploblock_bounds[haploblock] = (new_start, new_end)
        return True
    return False

def try_add_snp_to_haploblock_noLD(haploblock, snp, bloc_start, bloc_end, maf_threshold, cp_cut, cp_mode, MAX_BP_SIZE_FOR_HAPLOBLOCK, haploblock_bounds):
    snp_pos = snp.get_snp().get_position()
    core_snp = haploblock.get_core_snp()

    if (snp_pos < bloc_start and bloc_start - snp_pos > MAX_BP_SIZE_FOR_HAPLOBLOCK) or \
       (snp_pos > bloc_end and snp_pos - bloc_end > MAX_BP_SIZE_FOR_HAPLOBLOCK):
        return False

    maf_cut = maf_threshold * core_snp.get_maf() if maf_threshold else None
    if (maf_threshold is None or snp.get_maf() >= maf_cut) and \
       (cp_cut is None or calculate_percentage_carrier(core_snp, snp)[{'rough': 0, 'exact': 1}[cp_mode]] >= cp_cut):
        
        new_start = min(bloc_start, snp_pos)
        new_end = max(bloc_end, snp_pos)
        if new_end - new_start <= MAX_BP_SIZE_FOR_HAPLOBLOCK:
            haploblock.get_all_snps().append(snp)
            haploblock_bounds[haploblock] = (new_start, new_end)
            return True
    return False

def build_haploblocks_classic(snps: list[SNPInformations], r_square_cut: float, d_prime_cut: float, maf_threshold: float| None, cp_cut: float| None, 
                     cp_mode: Literal['rough', 'exact'] | None, MAX_BP_SIZE_FOR_HAPLOBLOCK: int | None, LD: int, region: int |None = -1) -> list[Haploblock]:
    """
    Build haploblocks from a list of SNPs based on LD, MAF threshold, and carrier percentage.
    
    Parameters:
    - snps: List of SNPInformations to process.
    - r_square_cut: float, parameter for LD calculations.
    - d_prime_cut: float, parameter for LD calculations.
    - maf_threshold: float, MAF threshold for SNP inclusion.
    - cp_cut: float, carrier percentage cut.
    - cp_mode: str, mode for calculating carrier percentage.
    - MAX_BP_SIZE_FOR_HAPLOBLOCK: int, maximum size for haploblock boundaries.
    - run_mode: int, 1 for LD-based processing, 2 for non-LD processing.
    
    Returns:
    - List of Haploblocks formed from the SNPs.
    """
    
    # Main function body
    snps_to_process = snps[:]
    processed_snps = []

    if not all(snp.is_phased() for snp in snps_to_process):
        logging.warning("Warning: not all data is phased. Using 'rough' carrier percentage.")
        cp_mode = 'rough'

    if LD == 1:
        if not r_square_cut and not d_prime_cut and not maf_threshold and not cp_cut:
            raise ValueError('At least one criterion of linkage has to be provided')

        if cp_cut and not cp_mode:
            logging.warning("Warning: cp cut but no cp mode: set to rough by default.")
            cp_mode = 'rough'

    if LD == 2:
        if not maf_threshold and not cp_cut:
            raise ValueError('At least one criterion of inclusion has to be provided')

    logging.info(f"Standard mode : {len(snps)} SNPs to process in the current region...")

    # Sort SNPs by ascending MAF
    snps_to_process.sort(key=lambda snp: snp.get_maf())
    haploblocks = []
    haploblock_bounds = {}

    # Main loop to process each SNP
    while snps_to_process:
        snp_processing = snps_to_process.pop(0)
        snp_pos = snp_processing.get_snp().get_position()
        snp_added = False

        for haploblock in haploblocks:
            bloc_start, bloc_end = haploblock_bounds[haploblock]
            if LD == 1:
                if try_add_snp_to_haploblock_withLD(
                    haploblock, snp_processing, bloc_start, bloc_end,
                    r_square_cut, d_prime_cut, maf_threshold,
                    MAX_BP_SIZE_FOR_HAPLOBLOCK, cp_cut, cp_mode, haploblock_bounds
                ):
                    snp_added = True

            elif LD == 2:
                if try_add_snp_to_haploblock_noLD(haploblock, snp_processing, bloc_start, bloc_end, 
                                                   maf_threshold, cp_cut, cp_mode, 
                                                   MAX_BP_SIZE_FOR_HAPLOBLOCK, haploblock_bounds):
                    snp_added = True

        # Create a new haploblock if SNP does not fit into any existing one
        if not snp_added:
            new_haploblock = Haploblock(snp_processing, [snp_processing], region)
            haploblocks.append(new_haploblock)
            haploblock_bounds[new_haploblock] = (snp_pos, snp_pos)

            # Re-scan processed SNPs for possible inclusion in the new haploblock
            for processed_snp in processed_snps:
                processed_snp_pos = processed_snp.get_snp().get_position()
                if abs(snp_pos - processed_snp_pos) <= MAX_BP_SIZE_FOR_HAPLOBLOCK:
                    if LD == 1:
                        try_add_snp_to_haploblock_withLD(
                            new_haploblock, processed_snp, snp_pos, snp_pos,
                            r_square_cut, d_prime_cut, maf_threshold,
                            MAX_BP_SIZE_FOR_HAPLOBLOCK, cp_cut, cp_mode, haploblock_bounds
                        )

                    elif LD == 2:
                        try_add_snp_to_haploblock_noLD(
                            new_haploblock, processed_snp, snp_pos, snp_pos, 
                            maf_threshold, cp_cut, cp_mode, 
                            MAX_BP_SIZE_FOR_HAPLOBLOCK, haploblock_bounds
                        )

        processed_snps.append(snp_processing)

    logging.info(f"nb haps : {len(haploblocks)}")
    return haploblocks

def build_haploblocks_listSNP(snps: list[SNPInformations], r_square_cut: float, d_prime_cut: float, maf_threshold: float| None, cp_cut: float| None,  
                     cp_mode: Literal['exact'] | None, MAX_BP_SIZE_FOR_HAPLOBLOCK: int | None, snp_list: list, all_hap: int, add_snp: bool, region_index: int |None = -1) -> list[Haploblock]:
    """
    Build haploblocks from a list of SNPs based on LD, MAF threshold, and carrier percentage.
    
    Parameters:
    - snps: List of SNPInformations to process.
    - r_square_cut: float, parameter for LD calculations.
    - d_prime_cut: float, parameter for LD calculations.
    - maf_threshold: float, MAF threshold for SNP inclusion.
    - cp_cut: float, carrier percentage cut.
    - cp_mode: str, mode for calculating carrier percentage.
    - MAX_BP_SIZE_FOR_HAPLOBLOCK: int, maximum size for haploblock boundaries.
    - LD: int, 1 for LD-based processing, 2 for non-LD processing.
    - snp_list: list, predefined SNPs of interest.
    - list_mode: int, mode for handling snp_list (see details below).

    list_mode:
    - 1: Classic mode + generate haploblocks for SNPs in snp_list (SNPs can be part of other haploblocks).
    - 2: Same as mode 1, but SNPs in snp_list will only generate their own haploblock and won't be added to other haploblocks.
    - 3: Only generates haploblocks for SNPs in snp_list (ignores other SNPs).

    Returns:
    - List of haploblocks formed from the SNPs.
    """

    # Main function body
    snps_to_process = snps[:]
    processed_snps = []
    all_hap= int(all_hap)
    snp_list_set = list(snp_list)
    cp_mode = cp_mode or 'exact'
    
    if cp_cut and not cp_mode:
        logging.warning("Warning: cp cut but no cp mode: set to exact by default.")
        cp_mode = 'exact'

    logging.info(f"List SNP mode: {len(snps)} SNPs to process in the current region...")

    # Sort SNPs by ascending MAF
    snps_to_process.sort(key=lambda snp: snp.get_maf())

    haploblocks = []
    haploblock_bounds = {}

    if all_hap == 2:
        # Filter SNPs that exactly match IDs in snp_list
        matching_snps = [snp for snp in snps_to_process if snp.get_snp().get_rs_id() in snp_list_set]
        logging.info(f"Matching SNPs count: {len(matching_snps)}")
        # Generate haploblocks for SNPs in snp_list
        for snp in snps_to_process:
            if snp.get_snp().get_rs_id() in snp_list_set:
                # Create a haploblock for the current SNP
                new_haploblock = Haploblock(snp, [snp], region_index)
                haploblock_bounds[new_haploblock] = (snp.get_snp().get_position(), snp.get_snp().get_position())
                haploblocks.append(new_haploblock)

        # Scan all SNPs and add them to existing haploblocks if they match criteria
        for haploblock in haploblocks:
            core_snp = haploblock.get_core_snp()
            bloc_start, bloc_end = haploblock_bounds[haploblock]
            
            for snp in snps_to_process:
                if snp == core_snp:  # Skip the core SNP of the current haploblock
                    continue

                if add_snp == False:
                    if snp.get_snp().get_rs_id() in snp_list_set :
                        continue
                
                snp_pos = snp.get_snp().get_position()

                try_add_snp_to_haploblock_withLD(
                    haploblock, snp, bloc_start, bloc_end,
                    r_square_cut, d_prime_cut, maf_threshold,
                    MAX_BP_SIZE_FOR_HAPLOBLOCK, cp_cut, cp_mode, haploblock_bounds
                )

        return haploblocks

    if all_hap == 1 :
        while snps_to_process:
            snp_processing = snps_to_process.pop(0)
            snp_id = snp_processing.get_snp().get_rs_id()
            snp_pos = snp_processing.get_snp().get_position()
            snp_added = False

            # Check if the SNP is in the predefined snp_list
            if snp_id in snp_list_set:
                # Create a new haploblock for the SNP as the core SNP
                new_haploblock = Haploblock(snp_processing, [snp_processing], region_index)
                new_haploblock_bounds = (snp_processing.get_snp().get_position(), snp_processing.get_snp().get_position())
                haploblocks.append(new_haploblock)  # Add the new haploblock to the list
                haploblock_bounds[new_haploblock] = new_haploblock_bounds  # Add bounds for the new haploblock
                
                # Skip further processing for list_mode == 2
                if add_snp == False:
                    continue

            for haploblock in haploblocks:
                bloc_start, bloc_end = haploblock_bounds[haploblock]
                if try_add_snp_to_haploblock_withLD(
                    haploblock, snp_processing, bloc_start, bloc_end,
                    r_square_cut, d_prime_cut, maf_threshold,
                    MAX_BP_SIZE_FOR_HAPLOBLOCK, cp_cut, cp_mode, haploblock_bounds
                ):
                    snp_added = True


            # Create a new haploblock if SNP does not fit into any existing one
            if not snp_added:
                new_haploblock = Haploblock(snp_processing, [snp_processing], region_index)
                haploblocks.append(new_haploblock)
                haploblock_bounds[new_haploblock] = (snp_pos, snp_pos)

                # Re-scan processed SNPs for possible inclusion in the new haploblock
                for processed_snp in processed_snps:
                    processed_snp_pos = processed_snp.get_snp().get_position()
                    if abs(snp_pos - processed_snp_pos) <= MAX_BP_SIZE_FOR_HAPLOBLOCK:
                        try_add_snp_to_haploblock_withLD(
                            new_haploblock, processed_snp, snp_pos, snp_pos,
                            r_square_cut, d_prime_cut, maf_threshold,
                            MAX_BP_SIZE_FOR_HAPLOBLOCK, cp_cut, cp_mode, haploblock_bounds
                        )

            processed_snps.append(snp_processing)
    logging.info(f"nb haps : {len(haploblocks)}")
    return haploblocks

def build_haploblocks_rescan(snps: list[SNPInformations], r_square_cut: float, d_prime_cut: float, maf_threshold: float | None, cp_cut: float | None,
                             cp_mode: Literal['rough', 'exact'] | None, MAX_BP_SIZE_FOR_HAPLOBLOCK: int | None, LD: int, region: int |None  = -1) -> list[Haploblock]:
    """
    Generates haploblocks for every SNP in the input list, scanning all SNPs for inclusion in each haploblock.
    """
    # Initial setup
    snps_to_process = snps[:]
   
    # Check if data is phased
    if not all(snp.is_phased() for snp in snps_to_process):
        logging.warning("Warning: not all data is phased. haploblocks will use rough carrier percentage.")
        r_square_cut, d_prime_cut = None, None
        cp_mode = 'rough'
   
    if not r_square_cut and not d_prime_cut and not maf_threshold and not cp_cut:
        raise ValueError('At least one criterion of linkage has to be provided')

    if cp_cut and not cp_mode:
        logging.warning("Warning: cp cut but no cp mode: set to rough by default.")
        cp_mode = 'rough'

    logging.info(f"Exhaustive mode: {len(snps)} SNPs to process in the current region...")

    haploblocks = []  # List of haploblocks
    haploblock_bounds = {}  # Tracks boundaries for each haploblock

    # Iterate over each SNP to create a new haploblock
    for snp_processing in snps_to_process:
        # Initialize a new haploblock with the current SNP
        snp_pos = snp_processing.get_snp().get_position()
        new_haploblock = Haploblock(snp_processing, [snp_processing], region)
        haploblock_bounds[new_haploblock] = (snp_pos, snp_pos)

        # Scan all SNPs for inclusion in the current haploblock
        for other_snp in snps_to_process:
            # Skip the SNP already included as the core SNP
            if other_snp == snp_processing:
                continue

            other_snp_pos = other_snp.get_snp().get_position()

            if abs(snp_pos - other_snp_pos) <= MAX_BP_SIZE_FOR_HAPLOBLOCK:
                if LD == 1:  # With LD
                    if try_add_snp_to_haploblock_withLD(
                        new_haploblock, other_snp, snp_pos, snp_pos,
                        r_square_cut, d_prime_cut, maf_threshold,
                        MAX_BP_SIZE_FOR_HAPLOBLOCK, cp_cut, cp_mode, haploblock_bounds
                    ):
                        pass

                elif LD == 2:  # Without LD
                    if try_add_snp_to_haploblock_noLD(
                        new_haploblock, other_snp, snp_pos, snp_pos,
                        maf_threshold, cp_cut, cp_mode,
                        MAX_BP_SIZE_FOR_HAPLOBLOCK, haploblock_bounds
                    ):
                        pass

        # Once the scan is complete, add the haploblock to the list
        haploblocks.append(new_haploblock)

    logging.info(f"Generated {len(haploblocks)} haploblocks.")
    return haploblocks

def filter_snp_chunk_local(haploblock_chunk, max_empty_snp_gap, dataset_snp_id):
    for haploblock in haploblock_chunk:
        haploblock.filter_snps_by_distance(max_empty_snp_gap, dataset_snp_id)
    return haploblock_chunk

def build_haploblocks_per_region(
    snps: list, REGION_SIZE: int, MAX_BP_SIZE_FOR_HAPLOBLOCK: int, 
    MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC: int, r_square_cut: float, d_prime_cut: float,
    maf_threshold: float | None, cp_cut: float | None, cp_mode: Literal['rough', 'exact'] | None, 
    ld_parameter: int, run_mode: int, list_snp: list | None, all_hap: int | None, add_snp: bool | None,
    extend_mode: bool | None = False, extend_threshold: float | None = 0.9, pruning_mode: bool | None = False, overlap_threshold: float | None = 0.95
) -> list:
    """
    Build haploblocks from SNPs divided into genomic regions using threads & multiprocessing, with CPU usage throttling and memory management.
    """
    run_mode = int(run_mode)
    ld_parameter = int(ld_parameter)

    logging.info("Slicing the data...")
    # Sort SNPs consistently based on position before splitting them into genomic regions
    snps_sorted = sorted(snps, key=lambda snp_info: snp_info.get_snp().get_position())
    
    # Split SNPs into genomic regions
    genomic_regions = split_genomic_regions(snps_sorted, REGION_SIZE, MAX_BP_SIZE_FOR_HAPLOBLOCK)

    haploblocks = []  # List to store all haploblocks
    total_regions = len(genomic_regions)  # Total number of genomic regions
    logging.info("Preparation of the haploblocks... (may take hours depending on the size of the dataset)")
    # Use a Queue to maintain determinism
    result_queue = Queue()
    TARGET_CPU_USAGE = 90  # Target CPU usage percentage to throttle usage

    def process_region(region: list, region_index: int):
        """
        Function to process a single genomic region and build haploblocks.
        """
        # Determine the appropriate processing function based on run_mode
        if run_mode == 1:  # Standard mode
            region_haploblocks = build_haploblocks_classic(
                region, r_square_cut, d_prime_cut, maf_threshold, cp_cut, cp_mode, MAX_BP_SIZE_FOR_HAPLOBLOCK, ld_parameter, region_index)
        elif run_mode == 2:  # Exhaustive mode
            region_haploblocks = build_haploblocks_rescan(
                region, r_square_cut, d_prime_cut, maf_threshold, cp_cut, cp_mode, MAX_BP_SIZE_FOR_HAPLOBLOCK, ld_parameter,region_index)
        elif run_mode == 3:  # List SNP mode
            region_haploblocks = build_haploblocks_listSNP(
                region, r_square_cut, d_prime_cut, maf_threshold, cp_cut, cp_mode, MAX_BP_SIZE_FOR_HAPLOBLOCK, list_snp, all_hap, add_snp, region_index)
        else:
            raise ValueError("Invalid run_mode.")

        # Safely add the region's haploblocks to the result queue
        result_queue.put((region_index, region_haploblocks))

        # Log progress
        thread_name = threading.current_thread().name
        logtime = utils.log_time()
        logging.info(f"{logtime} - [{thread_name}] Region n°{region_index + 1}/{total_regions} processed.")

        # Trigger garbage collection to free memory
        gc.collect()

    def throttle_cpu_usage(threads):
        """
        Function to monitor and throttle CPU usage if it exceeds the target limit.
        """
        while any(thread.is_alive() for thread in threads):
            current_cpu_usage = psutil.cpu_percent(interval=0.1)  # Check CPU usage
            if current_cpu_usage > TARGET_CPU_USAGE:
                time.sleep(0.1)  # Pause to reduce CPU load

    # Limit the number of threads to the available CPU cores or total regions
    max_threads = min(os.cpu_count()-1, total_regions)

    # Prepare threads to process regions concurrently
    region_threads = []
    start_time = time.time()
    for i, region in enumerate(genomic_regions):
        thread = threading.Thread(target=process_region, args=(region, i))
        region_threads.append(thread)

    # Start threads in batches limited by max_threads
    for i in range(0, len(region_threads), max_threads):
        batch = region_threads[i:i + max_threads]

        # Start the current batch of threads
        for thread in batch:
            thread.start()

        # Start the CPU throttling thread
        throttle_thread = threading.Thread(target=throttle_cpu_usage, args=(batch,))
        throttle_thread.start()

        # Wait for all threads in the current batch to finish
        for thread in batch:
            thread.join()

        # Ensure the throttling thread finishes
        throttle_thread.join()

    # Collect results from the queue
    haploblocks = []
    while not result_queue.empty():
        _, region_haploblocks = result_queue.get()
        haploblocks.extend(region_haploblocks)

    end_time = time.time()
    logging.info(f"Threading-based building of haploblocks complete in {end_time - start_time:.2f} seconds.")
    logging.info("Correcting the haploblocks using multiprocessing...")
    logging.info(f"Initial number of haploblocks: {len(haploblocks)}")

    # Define chunk size for correction
    num_workers = min(os.cpu_count()-1, total_regions)  # Limit workers to available CPUs
    chunk_size = max(1, len(haploblocks) // (num_workers * 10))
    haploblock_chunks = [haploblocks[i:i + chunk_size] for i in range(0, len(haploblocks), chunk_size)]
    filtered_results = []

    # Filtering haploblock chunks
    start_time = time.time()
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(filter_snp_chunk_local, chunk, MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC,
                            [snp.get_snp().get_rs_id() for snp in snps_sorted])
            for chunk in haploblock_chunks
        ]

        for future in as_completed(futures):
            try:
                result = future.result()
                filtered_results.extend(result)
            except Exception as e:
                logging.error(f"Error in chunk filtering: {e}")
            finally:
                gc.collect()

    logging.info(f"Multiprocessing correction completed in {time.time() - start_time:.2f}s\n")

    # --------------------------
    # Step 8: Optional extension of large haploblocks
    # --------------------------
    if extend_mode:
        print("Extending large haploblocks...")
        logging.info("Extending large haploblocks...")

        haploblocks_to_extend = [
            hb for hb in filtered_results
            if (
                max(snp.get_snp().get_position() for snp in hb.get_all_snps()) -
                min(snp.get_snp().get_position() for snp in hb.get_all_snps())
            ) >= extend_threshold * MAX_BP_SIZE_FOR_HAPLOBLOCK
        ]
        logging.info(f"Number of haploblocks eligible for extension: {len(haploblocks_to_extend)}")
        print(f"Number of haploblocks eligible for extension: {len(haploblocks_to_extend)}")
        extended_haploblocks, stats_list = extend_haploblocks(
            haploblocks_to_extend,
            genomic_regions,
            r_square_cut,
            d_prime_cut,
            maf_threshold,
            cp_cut,
            MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC
        )

        core_to_extended = {ehb.get_core_snp().get_snp().get_rs_id(): ehb for ehb in extended_haploblocks}
        final_haploblocks = [core_to_extended.get(hb.get_core_snp().get_snp().get_rs_id(), hb)
                            for hb in filtered_results]
        print(f"Extension complete. Number of haploblocks extended: {len(extended_haploblocks)}")
        logging.info(f"Extension complete. Number of haploblocks extended: {len(extended_haploblocks)}\n")
    else:
        final_haploblocks = filtered_results

    # --------------------------
    # Step 9: Optional pruning of overlapping haploblocks (will remove redundant blocks by keeping the largest)
    # --------------------------

    if pruning_mode:
        logging.info(f"Number of haploblocks before pruning: {len(final_haploblocks)}")
        print(f"Number of haploblocks before pruning: {len(final_haploblocks)}")
        final_haploblocks = prune_redundant_haploblocks(
            final_haploblocks,
            overlap_threshold,
            coordinate_margin_bp=0
        )
        print(f"Pruning complete. Number of haploblocks after pruning: {len(final_haploblocks)}")
        logging.info(f"Final number of haploblocks after pruning: {len(final_haploblocks)}\n")
    else:
        print(f"Final number of haploblocks: {len(final_haploblocks)}")
        logging.info(f"Final number of haploblocks: {len(final_haploblocks)}\n")
    
    gc.collect()  # Release memory
    return final_haploblocks

def caracterize_ld_with_reference(ref_snp : SNPInformations, snp_test : SNPInformations, r_square_cut : float | None, d_prime_cut : float | None, maf_cut : float | None, cp_cut : float | None, cp_mode : Literal['rough', 'exact'] | None) -> tuple[bool, float, float, float, float] :
    """
    For a reference SNP and a SNP to compare, compute their LD stats.
    If they are in sufficient LD considering the given parameters, returns True and the different statistics.
    Else, returns False followed by zeros.
    """
    # Don't test a snp with itself.
    if ref_snp.get_snp().get_rs_id() == snp_test.get_snp().get_rs_id() :
        return False, 0, 0, 0, 0

    # If a maf threshold is specified, don't test snps that have a maf below.
    if maf_cut and snp_test.get_maf() < maf_cut :
        return False, 0, 0, 0, 0

    r_square, d, d_prime, cp = 0, 0, 0, 0

    # If the data is not phased for both SNPs, only rough cp can be computed.
    if not ref_snp.is_phased() or not snp_test.is_phased() :
        r_square_cut, d_prime_cut = None, None
        cp_mode = 'rough'

    # If r^2 or D' cuts are specified, calculate r^2, D and D' and skip the SNP if r^2 or D' are too low.
    if r_square_cut or d_prime_cut or cp_mode == 'exact':
        r_square, d, d_prime = calculate_ld_between_two_snps(ref_snp, snp_test)
        if (r_square_cut and r_square < r_square_cut) or (d_prime_cut and abs(d_prime) < d_prime_cut) :
            return False, 0, 0, 0, 0

    # If a cp threshold is specified, calculate the requested carrier_percentage and skip snps that are below.
    if cp_cut and not cp_mode :
        logging.warning("Warning : cp cut but no cp mode : set to rough by default.")

    rough_cp, exact_cp, _, _ = calculate_percentage_carrier(ref_snp, snp_test)

    if cp_mode == 'rough' and cp_cut and rough_cp < cp_cut :
        return False, 0, 0, 0, 0
    if cp_mode == 'rough' :
        cp = rough_cp
    if cp_mode == 'exact' and cp_cut and exact_cp < cp_cut :
        return False, 0, 0, 0, 0
    if cp_mode == 'exact' :
        cp = exact_cp

    # If the snp is still here, it means it is in high ld with the ref snp according to the given criterion.
    return True, r_square, d, d_prime, cp

def _extend_single_haploblock_worker(
    index: int,
    haploblock,
    genomic_regions: List[List["SNPInformations"]],
    r_square_cut: float | None,
    d_prime_cut: float | None,
    maf_threshold: float | None,
    cp_cut: float | None
) -> Tuple[int, "Haploblock", str, int, int, int, int, int, int, int]:
    """
    Worker function for extending a single haploblock.

    Extension logic:
    ----------------
    - Start with regions [n-2 ... n+2].
    - If haploblock is close to borders (≤10% of n-2 region start or ≥90% of n+2 region end),
      extend further: [n-4 ... n+4].
    - Repeat recursively until:
        * no border condition is met, OR
        * the chromosome edges are reached (region 0 or last region).
    - Within candidate regions, SNPs are tested for LD with the core SNP,
      and added if they pass LD thresholds.
    """
    try:
        cp_mode = "exact"
        core_snp_info = haploblock.get_core_snp()
        core_rs_id = core_snp_info.get_snp().get_rs_id()
        region_index = haploblock.get_region()

        maf_cut = maf_threshold * core_snp_info.get_maf() if maf_threshold else None

        before_snps = haploblock.get_snp_size()
        before_bp = haploblock.get_bp_size()

        extended_snps = haploblock.get_all_snps()[:]
        existing_rsids = {s.get_snp().get_rs_id() for s in extended_snps}

        tested, added, rejected_ld = 0, 0, 0

        # --- Step 1: initial window [n-2 ... n+2] ---
        total_regions = len(genomic_regions)
        expansion = 2
        hb_start, hb_end = haploblock.get_start_bp(), haploblock.get_end_bp()

        while True:
            start_region = max(0, region_index - expansion)
            end_region = min(total_regions - 1, region_index + expansion)
            candidate_regions = list(range(start_region, end_region + 1))

            # Check border proximity
            expand_needed = False

            if genomic_regions[start_region]:
                left_start = genomic_regions[start_region][0].get_snp().get_position()
                left_end = genomic_regions[start_region][-1].get_snp().get_position()
                if hb_start <= left_start + 0.1 * (left_end - left_start):
                    expand_needed = True

            if genomic_regions[end_region]:
                right_start = genomic_regions[end_region][0].get_snp().get_position()
                right_end = genomic_regions[end_region][-1].get_snp().get_position()
                if hb_end >= right_end - 0.1 * (right_end - right_start):
                    expand_needed = True

            # Stop if no need to expand further or reached chromosome edges
            if not expand_needed or (start_region == 0 and end_region == total_regions - 1):
                break

            # Increase expansion window by 2 and retry
            expansion += 2

        # --- Step 2: scan candidate SNPs ---
        for ridx in candidate_regions:
            for candidate in genomic_regions[ridx]:
                cand_rsid = candidate.get_snp().get_rs_id()
                if cand_rsid in existing_rsids:
                    continue

                tested += 1
                is_in_ld, *_ = caracterize_ld_with_reference(
                    core_snp_info, candidate,
                    r_square_cut, d_prime_cut, maf_cut, cp_cut, cp_mode
                )

                if is_in_ld:
                    extended_snps.append(candidate)
                    existing_rsids.add(cand_rsid)
                    added += 1
                else:
                    rejected_ld += 1

        # Build new haploblock with extended SNPs
        new_haploblock = Haploblock(core_snp_info, extended_snps, region_index)

        after_snps = new_haploblock.get_snp_size()
        after_bp = new_haploblock.get_bp_size()

        return (index, new_haploblock, core_rs_id, before_snps, before_bp,
                after_snps, after_bp, tested, added, rejected_ld)

    except Exception as e:
        logging.error(f"Error in worker for haploblock index {index}: {e}", exc_info=True)
        return index, None, "", 0, 0, 0, 0, 0, 0, 0

def extend_haploblocks(haploblocks: List["Haploblock"], genomic_regions: List[List["SNPInformations"]], r_square_cut: float | None, d_prime_cut: float | None,
    maf_threshold: float | None, cp_cut: float | None, MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC: int | None) -> Tuple[List["Haploblock"], List[Tuple[str, int, int, int, int, int, int, int]]]:
    """
    Parallelized extension of haploblocks across multiple processes.

    Workflow:
    ---------
    1. Split haploblocks into chunks and distribute across workers.
    2. Each worker runs `_extend_single_haploblock_worker` to locally extend haploblocks
       using only SNPs in nearby regions (see worker for rules).
    3. Collect extended haploblocks and statistics.
    4. Post-process extended haploblocks:
       - Remove large internal gaps based on MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC.

    Parameters
    ----------
    haploblocks : List[Haploblock]
        Input haploblocks to extend.
    genomic_regions : List[List[SNPInformations]]
        Genomic regions (list of SNPs per region).
    r_square_cut, d_prime_cut, maf_threshold, cp_cut : thresholds for LD and filtering.
    MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC : int
        Max allowed gap size; larger gaps trigger haploblock cutting.

    Returns
    -------
    recut_results : List[Haploblock]
        Extended and post-processed haploblocks.
    stats_list : List[Tuple]
        Per-haploblock statistics (core SNP ID, before/after sizes, tested/added/rejected counts).
    """
    if not haploblocks:
        return [], []

    total = len(haploblocks)
    num_workers = min(os.cpu_count() - 1, total)
    chunk_size = max(1, len(haploblocks) // (num_workers * 10))

    logging.info(f"Extending {total} haploblocks using {num_workers} workers (chunk_size={chunk_size})")

    results_by_index: List[Optional["Haploblock"]] = [None] * total
    stats_list: List[Tuple[str, int, int, int, int, int, int, int]] = [None] * total

    # Split haploblocks into balanced chunks
    haploblock_chunks = [
        haploblocks[i:i + chunk_size] for i in range(0, len(haploblocks), chunk_size)
    ]

    # --- Run extension in parallel workers ---
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(
                _extend_single_haploblock_worker,
                batch_start + idx, hb, genomic_regions,
                r_square_cut, d_prime_cut, maf_threshold, cp_cut
            ): (batch_start + idx)
            for batch_start, chunk in enumerate(haploblock_chunks)
            for idx, hb in enumerate(chunk)
        }

        for fut in as_completed(futures):
            idx = futures[fut]
            try:
                (index, new_haploblock, core_rs_id, before_snps, before_bp,
                 after_snps, after_bp, tested, added, rejected_ld) = fut.result()

                logging.info(
                    f"[EXTEND] Core SNP {core_rs_id}: "
                    f"{before_snps} SNPs/{before_bp} bp → {after_snps} SNPs/{after_bp} bp "
                    f"(tested={tested}, added={added}, rejected={rejected_ld})"
                )

                results_by_index[index] = new_haploblock
                stats_list[index] = (core_rs_id, before_snps, before_bp,
                                     after_snps, after_bp, tested, added, rejected_ld)

            except Exception as e:
                logging.error(f"[ERROR] Exception while extending haploblock index={idx}: {e}", exc_info=True)

    extended_haploblocks = [hb for hb in results_by_index if hb is not None]

    # --- Post-extension filtering to remove large internal SNP gaps ---
    logging.info("Re-cutting extended haploblocks to remove large internal gaps...")

    haploblock_chunks = [
        extended_haploblocks[i:i + chunk_size]
        for i in range(0, len(extended_haploblocks), chunk_size)
    ]

    recut_results = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(
                filter_snp_chunk_local,
                chunk,
                MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC,
                [snp.get_snp().get_rs_id() for region in genomic_regions for snp in region]
            )
            for chunk in haploblock_chunks
        ]

        for future in as_completed(futures):
            try:
                recut_results.extend(future.result())
            except Exception as e:
                logging.error(f"Error in post-extension filtering: {e}")

    logging.info(f"Re-cut complete. Final number of haploblocks: {len(recut_results)}")
    return recut_results, stats_list

def prune_redundant_haploblocks(haploblocks: list, overlap_threshold: float | None, coordinate_margin_bp: int = 0, num_threads: int | None = None) -> list:
    """
    Remove redundant haploblocks where one haploblock (H1) is mostly contained within
    another (H2) based on shared SNP content.

    A haploblock H1 is considered redundant if:
        - The genomic intervals [start1, end1] and [start2, end2] overlap 
          (with an optional margin in base pairs),
        - At least `overlap_threshold * SNPs(H1)` SNPs of H1 are also present in H2,
        - H2 contains at least as many SNPs as H1.
    """

    if not haploblocks:
        return haploblocks

    if not (0.0 < overlap_threshold <= 1.0):
        raise ValueError("overlap_threshold must be between 0 and 1.")

    # --- Step 1: Precompute metadata for each haploblock ---
    haploblock_metadata = []
    for index, haploblock in enumerate(haploblocks):
        snps = haploblock.get_all_snps()
        rsids = {snp.get_snp().get_rs_id() for snp in snps}  # unique SNP IDs
        haploblock_metadata.append({
            "index": index,
            "haploblock": haploblock,
            "start": haploblock.get_start_bp(),
            "end": haploblock.get_end_bp(),
            "snp_set": rsids,
            "size": len(rsids)
        })

    # Sort haploblocks by genomic start for efficient interval scanning
    haploblock_metadata.sort(key=lambda m: m["start"])

    # Extract arrays for faster access in workers
    starts = [m["start"] for m in haploblock_metadata]
    ends = [m["end"] for m in haploblock_metadata]
    snp_sets = [m["snp_set"] for m in haploblock_metadata]
    sizes = [m["size"] for m in haploblock_metadata]
    haploblocks_sorted = [m["haploblock"] for m in haploblock_metadata]

    n = len(haploblock_metadata)
    is_redundant = [False] * n  # tracks whether each haploblock is redundant

    # --- Step 2: Define helper to check interval overlap ---
    def intervals_overlap(i: int, j: int) -> bool:
        """
        Return True if intervals [start_i, end_i] and [start_j, end_j]
        overlap within the coordinate margin.
        """
        return max(starts[i] - coordinate_margin_bp, starts[j] - coordinate_margin_bp) <= \
               min(ends[i] + coordinate_margin_bp, ends[j] + coordinate_margin_bp)

    # --- Step 3: Precompute candidate ranges using sweep-line algorithm ---
    # For each haploblock, only consider candidates whose start is before the end+margin.
    candidate_ranges: list[tuple[int, int]] = [(-1, -1)] * n
    right_ptr = 0
    for i in range(n):
        while right_ptr < n and starts[right_ptr] <= ends[i] + coordinate_margin_bp:
            right_ptr += 1
        candidate_ranges[i] = (0, right_ptr - 1)

    # --- Step 4: Worker to detect redundancy in a range of haploblocks ---
    def worker(index_range: range):
        local_redundant = []
        for i in index_range:
            # Skip if already marked redundant
            if is_redundant[i]:
                continue

            size_i = sizes[i]
            if size_i == 0:
                continue

            snps_i = snp_sets[i]
            required_overlap = math.ceil(overlap_threshold * size_i)  # minimum SNP overlap required

            j_start, j_end = candidate_ranges[i]
            if j_start < 0 or j_end < 0:
                continue

            found_redundant = False
            for j in range(j_start, j_end + 1):
                if j == i:
                    continue
                if not intervals_overlap(i, j):
                    continue

                size_j = sizes[j]
                # H2 must be at least as large as H1 and large enough to contain required overlap
                if size_j < required_overlap or size_j < size_i:
                    continue

                # --- Early-stopping intersection count ---
                shared_snp_count = 0
                remaining = size_i
                for snp in snps_i:
                    if snp in snp_sets[j]:
                        shared_snp_count += 1
                        if shared_snp_count >= required_overlap:
                            found_redundant = True
                            break
                    remaining -= 1
                    # Stop if even the maximum possible overlap can't reach threshold
                    if shared_snp_count + remaining < required_overlap:
                        break

                if found_redundant:
                    break

            if found_redundant:
                local_redundant.append(i)

        return local_redundant

    # --- Step 5: Parallel execution with threads ---
    if num_threads is None:
        num_threads = max(1, (os.cpu_count() or 2) - 1)

    # Divide work into small chunks (each chunk is a range of indices)
    chunk_size = max(1, n // (num_threads * 8) or 1)
    index_ranges = [range(i, min(i + chunk_size, n)) for i in range(0, n, chunk_size)]

    # Launch thread workers
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(worker, r) for r in index_ranges]
        for future in as_completed(futures):
            try:
                redundant_indices = future.result()
                # Mark indices found redundant by this worker
                for idx in redundant_indices:
                    is_redundant[idx] = True
            except Exception as e:
                logging.error(f"Error in pruning worker: {e}")

    # --- Step 6: Collect non-redundant haploblocks ---
    pruned_haploblocks = [haploblocks_sorted[i] for i in range(n) if not is_redundant[i]]

    logging.info(
        f"Pruned haploblocks: {n - len(pruned_haploblocks)} removed, "
        f"{len(pruned_haploblocks)} kept (threshold={int(overlap_threshold*100)}% of similarity)."
    )

    # Free memory explicitly
    gc.collect()
    return pruned_haploblocks

