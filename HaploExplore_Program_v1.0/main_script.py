# Copyright (c) 2025 HaploExplore

from typing import Literal
import os
from genetic_structures import *
import files
import genetic_computations
import utils
import plotting
import json
import logging
from math import ceil

########## FUNCTIONS THAT THE USER CAN CALL ##########

logging.getLogger('matplotlib').setLevel(logging.WARNING)

### Constants ###
script_dir = os.path.dirname(os.path.realpath(__file__))
config_path = os.path.join(script_dir, "config.json")

with open(config_path, 'r') as file:
    config = json.load(file)

REGION_SIZE = config["REGION_SIZE"][0]
MAX_BP_SIZE_FOR_haploblock = config["MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC"][0]
maf_threshold = config["maf_threshold"][0]
REF_VCF_FILEPATH = "reference_vcf_pathfile_to_use" 

def find_haploblocks_for_given_snps():
    """
    Asks the user for input paths for SNP file, haploblock composition file, and output directory, 
    and extracts only the haploblocks associated with each SNP in the list. Writes results immediately for each SNP
    without including SNPs from the haploblocks in the output.

        Workflow:
    1. Ask the user for an output directory where results will be saved.
       - If none is provided, defaults to the current directory.
       - If the directory does not exist, it will be created.

    2. Ask the user for the SNP file path.
       - Must exist and contain a list of SNPs to analyze.
       - Raise an error if not provided or not found.

    3. Ask the user for the haploblock composition file path.
       - Must exist and contain haploblock definitions.
       - Raise an error if not provided or not found.

    4. Generate an output file name based on the haploblock composition file name.
       - Replace "haploblocks_comp" in the filename with "SNPs_carried".

    5. Read SNPs from the SNP file into memory.

    6. For each SNP:
       - Write SNP label in the output file.
       - Search the haploblock composition file twice:
         (a) First pass: identify haploblocks where the SNP is the **tag SNP** (first SNP in the list).
         (b) Second pass: identify haploblocks where the SNP is present anywhere (not only as a tag).
       - Write results immediately for each SNP.
    """

    # ===========================
    # Step 1: Get and validate output directory
    # ===========================
    work_directory = input("Enter the path of the directory you want to output computations: ")
    logging.info("Checking directory…")
    if not work_directory.strip():  # Case: user pressed Enter → use current directory
        work_directory = os.curdir
        logging.info(f"No directory provided: set to current directory {work_directory}.")
    elif not os.path.exists(work_directory):  # Case: directory does not exist → create it
        os.makedirs(work_directory)
        logging.info(f"Created and set to {work_directory}.")
    else:  # Case: directory exists
        logging.info(f"Set to {work_directory}.")

    # ===========================
    # Step 2: Get SNP file path
    # ===========================
    snp_file_path = input("Enter the path for the SNPs file to analyze: ")
    logging.info("Checking SNP file…")
    if not snp_file_path.strip():  # Empty path
        raise ValueError("You did not provide a path for the SNP file.")
    elif not os.path.exists(snp_file_path):  # Path does not exist
        raise FileNotFoundError("The SNP file does not exist.")
    else:  # Path exists
        logging.info("SNP file found.")

    # ===========================
    # Step 3: Get haploblock composition file path
    # ===========================
    haploblock_file_path = input("Enter the path for the haploblock composition file: ")
    logging.info("Checking haploblock composition file…")
    if not haploblock_file_path.strip():  # Empty path
        raise ValueError("You did not provide a path for the haploblock composition file.")
    elif not os.path.exists(haploblock_file_path):  # Path does not exist
        raise FileNotFoundError("The haploblock composition file does not exist.")
    else:  # Path exists
        logging.info("haploblock composition file found.")

    # ===========================
    # Step 4: Generate output file name
    # ===========================
    haploblock_file_name = os.path.basename(haploblock_file_path)
    output_file_name = haploblock_file_name.replace("haploblocks_comp", "SNPs_carried")
    output_file_path = os.path.join(work_directory, output_file_name)
    logging.info(f"Output will be saved in {output_file_path}.")

    # ===========================
    # Step 5: Read SNPs from SNP file
    # ===========================
    with open(snp_file_path, 'r') as snp_file:
        snps_to_analyze = [line.strip() for line in snp_file.readlines()]  # Clean whitespace

    # ===========================
    # Step 6: Analyze each SNP
    # ===========================
    with open(output_file_path, 'w') as output_file:
        for snp in snps_to_analyze:
            # Write header for the SNP in the output
            output_file.write(f"SNP: {snp}\n")

            # ---- First pass: check where SNP is a tag SNP ----
            tag_snp_haploblocks = []
            with open(haploblock_file_path, 'r') as haplo_file:
                for line in haplo_file:
                    if line.startswith("haploblock ID"):  # Skip header
                        continue

                    # Parse line: first field = haploblock ID, second field = SNPs in haploblock
                    haploblock_id, snps_in_haploblock = line.strip().split('\t', 1)
                    snps_in_haploblock_list = snps_in_haploblock.split()

                    # A tag SNP = the first SNP in the haploblock list
                    if snp == snps_in_haploblock_list[0]:
                        tag_snp_haploblocks.append(haploblock_id)

            # Write results for tag SNP case
            if tag_snp_haploblocks:
                output_file.write("\thaploblocks where this SNP is a tag_SNP:\n")
                for haploblock_id in tag_snp_haploblocks:
                    output_file.write(f"\t\t{haploblock_id}\n")
            else:
                output_file.write("\tNo haploblock where this SNP is a tag_SNP.\n")

            # ---- Second pass: check where SNP is present in any haploblock ----
            haplotype_found = False
            with open(haploblock_file_path, 'r') as haplo_file:
                for line in haplo_file:
                    if line.startswith("haploblock ID"):  # Skip header
                        continue

                    haploblock_id, snps_in_haploblock = line.strip().split('\t', 1)
                    snps_in_haploblock_list = snps_in_haploblock.split()

                    if snp in snps_in_haploblock_list:  # SNP is present
                        haplotype_found = True
                        output_file.write(f"\t{haploblock_id}\n")

                # If no haploblock contained the SNP
                if not haplotype_found:
                    output_file.write("\tNo haploblock found for this SNP.\n")

    logging.info(f"haploblocks associated with SNPs have been saved to {output_file_path}.")

def list_snps_carrying_given_snps(carrier_percentage_cut: float, maf_threshold: float | None):
    """
    For each SNP in a given list, compute the percentage of individuals carrying the tested allele
    who also carry the reference allele in a genomic region.

        Workflow:
    1. Output directory setup → ask user for working directory and initialize logging
    2. Input VCF file → ask for source VCF file or fall back to reference
    3. Region selection (optional) → allow extraction of a subset of SNPs from the VCF
    4. SNP list to analyze (mandatory) → provide file with SNPs of interest
    5. VCF data processing → extract SNPs info for both region and SNP list
    6. Computation → for each SNP in the list, find carriers within the region and output results
    """

    # ===============================================================
    # 1. Output directory setup
    # ===============================================================
    work_directory = input("Enter the path of the directory you want to output computations : ")

    # Check that the directory is valid; if not, use current directory or create it
    logging.info("Checking directory…")
    if not work_directory.strip():
        work_directory = os.curdir
        logging.info(f"No directory provided : set to current directory {work_directory}.")
    elif not os.path.exists(work_directory):
        os.makedirs(work_directory)
        logging.info(f"Created and set to {work_directory}.")
    else:
        logging.info(f"Set to {work_directory}.")

    # Setup logging with dynamic log filename
    log_name = "List_snps_carrying_given_snps_log"
    log_filepath = utils.generate_log_filename(work_directory, log_name)
    utils.setup_logging(log_filepath)
    logging.info(f"Logging started. All output will be saved to {log_filepath}\n")

    # ===============================================================
    # 2. Input VCF file
    # ===============================================================
    source_vcf_filepath = input("Enter the path for the VCF file containing your data : ")

    logging.info("Checking source file…")
    if not source_vcf_filepath.strip():
        # No file provided → confirm fallback to reference data
        if not utils.ask_confirmation(
            "Warning : you did not provide a source vcf. Reference data will be used. Do you want to continue ? (yes/no) :"
        ):
            logging.info("You chose not to continue.")
            return
        source_vcf_filepath = REF_VCF_FILEPATH
    elif not files.is_vcf_file(source_vcf_filepath):
        # Invalid file → confirm fallback to reference data
        if not utils.ask_confirmation(
            "Warning : you did not provide a valid format for your vcf data. Reference data will be used. Do you want to continue ? (yes/no) :"
        ):
            logging.info("You chose not to continue.")
            return
        source_vcf_filepath = REF_VCF_FILEPATH
    else:
        logging.info("OK.")

    # ===============================================================
    # 3. Region selection (optional)
    # ===============================================================
    if utils.ask_confirmation("Do you want to explore a specific region of your source data ? (yes/no) :"):
        region_to_extract = input("Enter the path for the file containing the rsIDs of the region you want to explore : ")

        logging.info("Checking region to extract…")
        if not region_to_extract.strip():
            if not utils.ask_confirmation(
                "Warning : you did not provide a region to explore. All SNPs will be analyzed. Do you want to continue ? (yes/no) :"
            ):
                logging.info("You chose not to continue.")
                return
            path_to_vcf_with_all_work_snps = source_vcf_filepath
        elif not files.is_list_of_rsIds(region_to_extract):
            if not utils.ask_confirmation(
                "Warning : the region you want to explore is not a list of rsIds. All SNPs will be analyzed. Do you want to continue ? (yes/no) :"
            ):
                logging.info("You chose not to continue.")
                return
            path_to_vcf_with_all_work_snps = source_vcf_filepath
        else:
            logging.info("Extracting SNPs from the source vcf…")
            path_to_vcf_with_all_work_snps = files.extract_snps_from_vcf(
                source_vcf_filepath, region_to_extract, work_directory + '/extracted_region.vcf'
            )
            logging.info("OK.")
    else:
        logging.info("All SNPs will be analyzed.")
        path_to_vcf_with_all_work_snps = source_vcf_filepath

    # ===============================================================
    # 4. SNP list to analyze (mandatory)
    # ===============================================================
    file_with_significant_snps = input("Enter the path for the file containing the SNPs you want to analyze : ")

    logging.info("Checking SNPs to analyze…")
    if not file_with_significant_snps.strip():
        raise ValueError("You did not provide a path for the list of SNPs to analyze.")

    # SNP list can be a VCF or a list of rsIDs
    if files.is_vcf_file(file_with_significant_snps):
        path_to_vcf_with_snps_to_analyze = file_with_significant_snps
        logging.info("OK.")
    elif files.is_list_of_rsIds(file_with_significant_snps):
        path_to_vcf_with_snps_to_analyze = files.extract_snps_from_vcf(
            source_vcf_filepath, file_with_significant_snps, work_directory + '/snps_to_analyze.vcf'
        )
        logging.info("OK.")
    else:
        raise ValueError("The list of SNPs to analyze must either be a vcf or contain IDs only.")

    # ===============================================================
    # 5. VCF data processing
    # ===============================================================
    logging.info("Processing vcf and extracting SNPs data…")
    all_work_snps_data = files.retrieve_snp_data_from_vcf(path_to_vcf_with_all_work_snps, maf_threshold)
    snps_to_analyze_data = files.retrieve_snp_data_from_vcf(path_to_vcf_with_snps_to_analyze, maf_threshold)
    logging.info("OK.")

    # ===============================================================
    # 6. Computation
    # For each SNP of interest:
    # - Create output subdirectory
    # - Find SNPs in the region that carry it
    # - Save and log results
    # ===============================================================
    for snp_to_analyze in snps_to_analyze_data:
        # Create a SNP-specific subdirectory
        snp_directory = work_directory + '/' + snp_to_analyze.get_snp().get_rs_id()
        if not os.path.exists(snp_directory):
            os.makedirs(snp_directory)

        # Extract carriers and save list
        snps_carriers, path_to_list_of_snps_carriers = genetic_computations.extract_snps_carried_by_reference(
            snp_to_analyze, all_work_snps_data, snp_directory, carrier_percentage_cut
        )

        # Logging the result
        if not snps_carriers:
            logging.info(
                f"List of SNPs carrying SNP {snp_to_analyze.get_snp().get_rs_id()} has been created in file {path_to_list_of_snps_carriers} :\n"
                f"it is empty"
            )
        else:
            logging.info(
                f"List of SNPs carrying SNP {snp_to_analyze.get_snp().get_rs_id()} has been created in file {path_to_list_of_snps_carriers} :\n"
                f"it contains {len(snps_carriers)} SNPs, between "
                f"{min(snp_data.get_snp().get_position() for snp_data in snps_carriers)} bp and "
                f"{max(snp_data.get_snp().get_position() for snp_data in snps_carriers)} bp."
            )

def find_haploblocks_in_a_region(r_square_cut: float, d_prime_cut: float, maf_percentage_cut: float | None, carrier_percentage_cut: float | None, 
    REGION_SIZE: int | None, MAX_BP_SIZE_FOR_HAPLOBLOCK: int | None, MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC: int | None, maf_threshold: float | None):
    """
    Takes a genomic region and looks for haploblocks within it.

            Workflow :
        1. Setup working directory & logging
        2. Load input VCF file (or fallback to reference data)
        3. (Optional) Extract a specific region from the VCF
        4. (Optional) Provide a list of SNPs for haploblock generation
        5. Choose run mode (Standard, Exhaustive, or List-based)
        6. Configure optional haploblock extension & pruning
        7. Process VCF and extract SNP data
        8. Build haploblocks based on SNP data and parameters
        9. Generate output files & visualizations
    """

    try: 
        # ==============================
        # 1. Setup working directory & logging
        # ==============================
        work_directory = input("Enter the path of the directory you want to output computations: ")
        logging.info("Checking directory…\n")

        if not work_directory.strip():
            # Default: current directory if nothing entered
            work_directory = os.curdir
            logging.info(f"No directory provided: set to current directory {work_directory}.\n")
        elif not os.path.exists(work_directory):
            # Create directory if missing
            os.makedirs(work_directory)
            logging.info(f"Created and set to {work_directory}.\n")
        else:
            logging.info(f"Set to {work_directory}.\n")

        # Setup logging file for this run
        log_name = "Find_Haploblock_log"
        log_filepath = utils.generate_log_filename(work_directory, log_name)
        utils.setup_logging(log_filepath)
        logging.info(f"Logging started. All output will be saved to {log_filepath}\n")

        logging.info(f"Using r_square_cut: {r_square_cut}, d_prime_cut: {d_prime_cut}, maf_percentage_cut: {maf_percentage_cut}, carriere_percentage_cut: {carrier_percentage_cut}, Max_empty_gap: {MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC}, Max_BP_size: {MAX_BP_SIZE_FOR_HAPLOBLOCK}, REGION_SIZE: {REGION_SIZE}, maf_threshold: {maf_threshold}")

        # ==============================
        # 2. Load input VCF file (or reference data)
        # ==============================
        source_vcf_filepath = input("Enter the path for the VCF file containing your data: ")
        logging.info("Checking source file…")

        if not source_vcf_filepath.strip():
            # Fallback to reference if nothing provided
            if not utils.ask_confirmation("Warning: you did not provide a source vcf. Reference data will be used. Do you want to continue? (yes/no):"):
                logging.info("You chose not to continue.")
                return
            source_vcf_filepath = REF_VCF_FILEPATH

        elif not files.is_vcf_file(source_vcf_filepath):
            # Fallback to reference if invalid format
            if not utils.ask_confirmation("Warning: you did not provide a valid format for your vcf data. Reference data will be used. Do you want to continue? (yes/no):"):
                logging.info("You chose not to continue.\n")
                return
            source_vcf_filepath = REF_VCF_FILEPATH

        else:
            logging.info(f"Set to {source_vcf_filepath}.\n")
            logging.info("OK.\n")

        # ==============================
        # 3. (Optional) Extract a specific region
        # ==============================
        if utils.ask_confirmation("Do you want to explore a specific region of your source data? (yes/no):"):
            region_to_extract = input("Enter the path for the file containing the rsIDs of the region you want to explore: ")
            logging.info("Checking region to extract…\n")

            if not region_to_extract.strip():
                # No region file provided → use all SNPs
                if not utils.ask_confirmation("Warning: you did not provide a region to explore. All SNPs will be analyzed. Do you want to continue? (yes/no):"):
                    logging.info("You chose not to continue.")
                    return
                path_to_vcf_with_all_work_snps = source_vcf_filepath

            elif not files.is_list_of_rs_ids(region_to_extract):
                # Invalid file → fallback to all SNPs
                if not utils.ask_confirmation("Warning: the region you want to explore is not a list of rsIds. All SNPs will be analyzed. Do you want to continue? (yes/no):"):
                    logging.info("You chose not to continue.\n")
                    return
                path_to_vcf_with_all_work_snps = source_vcf_filepath

            else:
                # Extract SNPs into a smaller region file
                logging.info("Extracting SNPs from the source vcf…\n")
                path_to_vcf_with_all_work_snps = files.extract_snps_from_vcf(
                    source_vcf_filepath,
                    region_to_extract,
                    os.path.join(work_directory, 'extracted_region.vcf')
                )
                logging.info("OK.")
        else:
            logging.info("All SNPs will be analyzed.\n")
            path_to_vcf_with_all_work_snps = source_vcf_filepath

        # ==============================
        # 4. (Optional) Provide SNP list
        # ==============================
        all_hap = None
        add_snp = None
        specific_snps = None

        if utils.ask_confirmation("Do you want to generate haploblocks for specific SNPs? (yes/no):"):
            # Ask for SNP list file
            specific_snps_file = input("Enter the file path for the SNPs file to generate haploblocks for: ")

            # Ensure valid path
            while not os.path.exists(specific_snps_file):
                logging.info(f"Error: The file {specific_snps_file} does not exist. Please try again.\n")
                specific_snps_file = input("Enter the file path for the SNPs file to generate haploblocks for: ")

            # Load SNP list from file
            with open(specific_snps_file, 'r') as f:
                specific_snps = [line.strip() for line in f.readlines() if line.strip()]
            logging.info(f"\n Loaded SNP list (5 firsts SNPs): {specific_snps[:5]}")

            run_mode = "3"  # SNP list mode
            # Ask processing options
            while all_hap not in ["1", "2"]:
                all_hap = input("\n Enter the list processing mode (1 = Standard + list SNPs, 2 = Only list SNPs): \n")
                if all_hap not in ["1", "2"]:
                    logging.info("Invalid choice. Please choose 1 or 2.")

                while add_snp not in ["yes","y", "no", "n"]:
                    add_snp = input("Do you want to add the SNPs in the list to other haploblocks? (yes/no): ").strip().lower()
                    if add_snp not in ["yes","y", "no", "n"]:
                        logging.info("Invalid choice. Please choose Yes or No.")
            add_snp = add_snp == "yes" or add_snp == "y"

        else:
            # If no SNP list → choose between Standard or Exhaustive mode
            logging.info("No SNPs list will be used.")
            run_mode = None
            while run_mode not in ["1", "2"]:
                run_mode = input("\n Enter the run mode (1 = Standard, 2 = Exhaustive): \n")
                if run_mode not in ["1", "2"]:
                    logging.info("Invalid choice. Please choose 1 or 2.")

        # Log run mode
        if run_mode == "1":
            logging.info("Runmode : Classic.\n")
        elif run_mode == "2":
            logging.info("Runmode: Exhaustive.\n")
        elif run_mode == "3":
            logging.info("Runmode : List_SNP.\n")

        carrier_percentage_mode = 'exact'
        ld_parameter = 1
        # ==============================
        # 5. Configure haploblock extension
        # ==============================
        extend_mode = None
        while extend_mode not in ["yes", "y", "no", "n"]:
            extend_mode = input("\nDo you want to extend haploblocks reaching max size by adding SNPs in high LD? (yes/no): ").strip().lower()
            if extend_mode not in ["yes", "y", "no", "n"]:
                logging.info("Invalid choice. Please type 'yes' or 'no'.")

        extend_mode = extend_mode in ["yes", "y"]

        if extend_mode:
            logging.info("Extension mode enabled.")

            while True:
                extend_threshold = input("Enter the extension threshold (0.0 - 1.0, default = 0.9): ").strip()
                try:
                    extend_threshold = float(extend_threshold)
                    if 0 <= extend_threshold <= 1.0:
                        break
                    else:
                        logging.info("Invalid threshold. Must be between 0.0 and 1.0.")
                except ValueError:
                    logging.info("Invalid input. Please enter a numeric value.")
        else:
            extend_threshold = None

        # ==============================
        # 6. Configure haploblock pruning
        # ==============================
        pruning_mode = None
        while pruning_mode not in ["yes", "y", "no", "n"]:
            pruning_mode = input("\nDo you want to prune redundant haploblocks? (yes/no): ").strip().lower()
            if pruning_mode not in ["yes", "y", "no", "n"]:
                logging.info("Invalid choice. Please type 'yes' or 'no'.")

        pruning_mode = pruning_mode in ["yes", "y"]

        if pruning_mode:
            logging.info("Pruning mode enabled.")

            while True:
                overlap_threshold = input("Enter pruning overlap threshold (0.0 - 1.0, default = 0.95): ").strip()
                try:
                    overlap_threshold = float(overlap_threshold)
                    if 0 <= overlap_threshold <= 1.0:
                        break
                    else:
                        logging.info("Invalid threshold. Must be between 0.0 and 1.0.")
                except ValueError:
                    logging.info("Invalid input. Please enter a numeric value.")
        else:
            overlap_threshold = None

        # ==============================
        # 7. Process VCF & SNP data
        # ==============================
        logging.info("Processing vcf and extracting SNPs data…\n")
        all_work_snps_data = files.retrieve_snp_data_from_vcf(path_to_vcf_with_all_work_snps, maf_threshold, specific_snps)

        # Verify SNPs are phased
        for snp_info in all_work_snps_data:
            if not snp_info.is_phased():
                logging.error("❌ SNP data contains unphased SNPs. Exiting.")
                return False

        if not all_work_snps_data:
            logging.error("❌ No SNP data extracted. Exiting.")
            return False

        logging.info("✅ Successfully extracted SNP data.\n")

        # ==============================
        # 8. Build haploblocks
        # ==============================
        logging.info("Building haploblocks…\n")
        """haploblocks = genetic_computations.build_haploblocks_per_region(
            all_work_snps_data, REGION_SIZE, MAX_BP_SIZE_FOR_HAPLOBLOCK, MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC, r_square_cut,
            d_prime_cut, maf_percentage_cut, carrier_percentage_cut, run_mode, specific_snps, all_hap, add_snp,
            extend_mode, extend_threshold, pruning_mode, overlap_threshold
        )"""

        haploblocks = genetic_computations.build_haploblocks_per_region(
                all_work_snps_data, REGION_SIZE, MAX_BP_SIZE_FOR_HAPLOBLOCK, MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC,
                r_square_cut, d_prime_cut, maf_percentage_cut, carrier_percentage_cut, carrier_percentage_mode,
                ld_parameter, run_mode, specific_snps, all_hap, add_snp,
                extend_mode, extend_threshold, pruning_mode, overlap_threshold
            )

        if not haploblocks:
            logging.info("No haploblocks generated, end of the program.")
            return

        # ==============================
        # 9. Generate output files & plots
        # ==============================
        additional_info = ""
        chromosome = haploblocks[0].get_all_snps()[0].get_snp().get_chromosome()
        if chromosome: 
            additional_info += f"chr_{chromosome}"

        # Append run mode to filename
        if run_mode == "1":
            additional_info += ".runMode_standard"
        elif run_mode == "2":
            additional_info += ".runMode_exhaustive"
        elif run_mode == "3":
            additional_info += ".runMode_listSNP"
        if extend_mode:
            additional_info += ".extended"
        if extend_threshold:
            additional_info += f"_{extend_threshold}"
        if pruning_mode:
            additional_info += ".pruned"
        if overlap_threshold:
            additional_info += f"_{overlap_threshold}"

        # Append thresholds to filename
        if r_square_cut: 
            additional_info += f".r2_{r_square_cut}"
        if d_prime_cut: 
            additional_info += f".D_prime_{d_prime_cut}"
        if maf_percentage_cut:
            additional_info += f".maf_{maf_percentage_cut}"
        if carrier_percentage_cut:
            additional_info += f".%carrier_{carrier_percentage_cut}"
        if REGION_SIZE:
            additional_info += f".reg_size_{REGION_SIZE}"
        if MAX_BP_SIZE_FOR_HAPLOBLOCK:
            additional_info += f".reg_overlap_{MAX_BP_SIZE_FOR_HAPLOBLOCK}"
        if MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC:
            additional_info += f".max_empty_gap_{MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC}"

        # Output haploblock data
        logging.info('Generating the haploblock information files...\n')
        path_to_haploblock_data = f"{work_directory}/haploblocks_data.{additional_info}.txt"
        path_to_haploblock_comp = f"{work_directory}/haploblocks_comp.{additional_info}.txt"
        path_to_haploblock_seq = f"{work_directory}/HAPLOBLOCKS_seq.{additional_info}.txt"
        path_to_snp_maf = f"{work_directory}/SNP_MAF.{additional_info}.txt"
        files.generate_haploblock_files(
            all_work_snps_data, haploblocks, path_to_vcf_with_all_work_snps,
            path_to_haploblock_data, path_to_haploblock_comp,
            path_to_haploblock_seq, path_to_snp_maf
        )
        logging.info('Haploblock information files generated.\n')

        # Output histograms
        logging.info('Generating the histograms...\n')
        path_to_comp_hist = f"{work_directory}/Histogram_haploblocks-SNP.{additional_info}.pdf"
        path_to_size_hist = f"{work_directory}/Histogram_haploblocks-Size.{additional_info}.pdf"
        plotting.generate_haploblock_histograms(path_to_haploblock_data, path_to_comp_hist, path_to_size_hist)
        logging.info('Histograms generated.')

    except Exception as e:
        logging.error(f"An error occurred: {e}")

def print_haploblocks():
    """
    Process haploblock data, apply filtering criteria, and generate plots.

    Workflow :
    1. Ask the user for a working directory.
       - If none is provided, defaults to the current directory.
       - Create directory if it does not exist.
       - Set up logging in this directory with a dynamic filename.

    2. Ask the user for input files.

    3. Load input files into memory.

    4. Ask the user for filtering criteria:
       - Minimum and maximum SNP counts per haploblock
       - Minimum and maximum base pair size of haploblocks
       - Base pair range to include
       - Minimum MAF threshold

    5. Filter haploblocks based on criteria.

    6. Plot the filtered haploblocks:
       - If total haploblocks exceed 50, split into multiple plots.
    """

    try:
        # ===========================
        # Step 1: Setup working directory and logging
        # ===========================
        log_name = "Print_Haploblocks_log"

        # Ask user for working directory
        work_directory = input("Enter the working directory: ")
        if not work_directory:
            work_directory = os.curdir  # Use current directory if none provided

        # Generate dynamic log filename and initialize logging
        log_filepath = utils.generate_log_filename(work_directory, log_name)
        utils.setup_logging(log_filepath)
        logging.info(f"Output will be saved to {log_filepath}\n")

        logging.info("Checking the directory...")
        if not os.path.exists(work_directory):
            os.makedirs(work_directory)  # Create directory if it doesn't exist
            logging.info(f"- Directory created: {work_directory}.")
        else:
            logging.info(f"- Selected directory: {work_directory}.")

        # ===========================
        # Step 2: Get input files and validate the format
        # ===========================
        haploblock_file = input("Enter the path to the haploblock file: ")
        snp_file = input("Enter the path to the SNP file: ")

        logging.info("Checking input files...")
        # Validate that files exist and have the correct format
        if not haploblock_file or not files.is_haploblock_composition(haploblock_file):
            raise ValueError("- Invalid haploblock composition file.")
        if not snp_file or not files.is_SNP_MAF(snp_file):
            raise ValueError("- Invalid SNP-MAF file format.")

        # ===========================
        # Step 3: Load input files
        # ===========================
        logging.info("Loading input files...")
        haploblock_composition = files.read_haploblock_composition(haploblock_file)  # dict: haploblock_id -> (core_snp, snps)
        snp_maf_dict = files.read_snp_maf_file(snp_file)  # dict: snp_id -> (maf, position)
        logging.info("Input files loaded.")

        # ===========================
        # Step 4: Ask user for filtering criteria
        # ===========================
        def get_positive_input(prompt):
            """
            Prompt user until a valid non-negative integer is entered.
            """
            while True:
                try:
                    value = int(input(prompt))
                    if value < 0:
                        print("Value must be greater than or equal to 0.")
                    else:
                        return value
                except ValueError:
                    print("Please enter a valid integer.")

        # Get filtering criteria from the user
        min_size_snp = get_positive_input("Enter minimum SNP size (0 for no limit): ")
        max_size_snp = get_positive_input("Enter maximum SNP size (0 for no limit): ")
        min_size_bp = get_positive_input("Enter minimum base pair size (0 for no limit): ")
        max_size_bp = get_positive_input("Enter maximum base pair size (0 for no limit): ")
        bp_range_start = get_positive_input("Enter base pair range start (0 for no limit): ")
        bp_range_end = get_positive_input("Enter base pair range end (0 for no limit): ")

        maf_threshold = float(input("Enter the minimum MAF threshold (0.0 - 1.0) - will highlight in red SNPs that have a MAF equal to +/- MAF threshold * MAF core_snp : "))

        # ===========================
        # Step 5: Filter haploblocks
        # ===========================
        filtered_haploblock_composition = {}
        for haploblock_id, (core_snp, snps) in haploblock_composition.items():
            # Get the genomic positions of all SNPs in this haploblock
            snp_positions = [snp_maf_dict[snp][1] for snp in snps if snp in snp_maf_dict]
            if not snp_positions:
                logging.warning(f"No valid SNP positions found for haploblock {haploblock_id}")
                continue

            bloc_start_bp = min(snp_positions)  # Start position of the haploblock
            bloc_end_bp = max(snp_positions)    # End position of the haploblock
            bloc_size_bp = bloc_end_bp - bloc_start_bp  # Size in base pairs
            bloc_size_snp = len(snps)                   # Number of SNPs

            # Apply all filters
            if min_size_snp and bloc_size_snp < min_size_snp:
                continue
            if max_size_snp and bloc_size_snp > max_size_snp:
                continue
            if min_size_bp and bloc_size_bp < min_size_bp:
                continue
            if max_size_bp and bloc_size_bp > max_size_bp:
                continue
            if bp_range_start and bloc_end_bp < bp_range_start:
                continue
            if bp_range_end and bloc_start_bp > bp_range_end:
                continue

            # Save filtered haploblock (tuple for hashable storage)
            filtered_haploblock_composition[haploblock_id] = (core_snp, tuple(snps))

        if not filtered_haploblock_composition:
            raise ValueError("No haploblocks left after applying filters.")

        # ===========================
        # Step 6: Plot haploblocks in chunks if too many
        # ===========================
        total_haploblocks = len(filtered_haploblock_composition)
        logging.info(f"Total haploblocks after filtering: {total_haploblocks}")

        max_haploblocks_per_plot = 50  # Split into multiple plots if necessary
        num_parts = ceil(total_haploblocks / max_haploblocks_per_plot)
        output_file = None  # Initialize output file path

        for part in range(num_parts):
            start_idx = part * max_haploblocks_per_plot
            end_idx = min((part + 1) * max_haploblocks_per_plot, total_haploblocks)

            # Select subset of haploblocks for this plot
            chunk_haploblock_composition = {
                haploblock_id: filtered_haploblock_composition[haploblock_id]
                for haploblock_id in list(filtered_haploblock_composition.keys())[start_idx:end_idx]
            }

            output_file = os.path.join(work_directory, f"haploblocks_part_{part + 1}.png")

            try:
                # Generate plot using external plotting module
                plot = plotting.plot_haploblocks(output_file, snp_maf_dict, chunk_haploblock_composition, maf_threshold)
                if plot:
                    logging.info(f"Graph generated in: {output_file}")
                else:
                    logging.error("Error generating graph for haploblocks.")
                    break
            except Exception as e:
                logging.error(f"Error plotting haploblocks: {e}")
                break

        return os.path.exists(output_file) if output_file else False

    except Exception as e:
        logging.error(f"An error occurred: {e}")



