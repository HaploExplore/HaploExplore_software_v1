from genetic_structures import *
import files 
import genetic_computations
import utils
import plotting

from io import StringIO
from typing import Literal
from math import ceil
import sys
import traceback
import os
import logging
import streamlit as st

########## FUNCTIONS THAT THE USER CAN CALL ##########

logging.getLogger('matplotlib').setLevel(logging.WARNING)

### Constants ###
script_dir = os.path.dirname(os.path.realpath(__file__))

REF_VCF_FILEPATH = "data/reference_vcf_pathfile_to_use" 

def find_haploblocks_for_given_snps(work_directory: str, snp_file_path: str, haploblock_file_path: str):
    """
    Identifies haploblocks associated with each SNP from a given list and writes the results to an output file.
    """

    # Redirects standard output to capture print statements (for debugging purposes)
    output_capture = StringIO()
    sys.stdout = output_capture

    try:
        # Ensure the working directory is valid
        st.write("**Checking if work directory is valid...**")
        if not work_directory.strip():  # If empty, use the current directory
            work_directory = os.curdir
        elif not os.path.exists(work_directory):  # If it doesn't exist, create it
            os.makedirs(work_directory)
        st.write(f"- Using work directory: {work_directory}")

        # Verify the SNP file exists
        st.write("**Checking if SNP file exists...**")
        if not os.path.exists(snp_file_path):
            raise FileNotFoundError(f"The SNP file '{snp_file_path}' does not exist.")
        st.write("- SNP file found.")

        # Verify the Haploblock composition file exists
        st.write("**Checking if Haploblock composition file exists...**")
        if not os.path.exists(haploblock_file_path):
            raise FileNotFoundError(f"The haploblock composition file '{haploblock_file_path}' does not exist.")
        st.write("- Haploblock file found.")

        # Define the output file name by replacing a part of the input file name
        output_file_name = os.path.basename(haploblock_file_path).replace("HAPLOBLOCKS_comp", "SNPs_carried")
        output_file_path = os.path.join(work_directory, output_file_name)

        # Load the list of SNPs to analyze from the SNP file
        st.write("**Loading SNPs to analyze...**")
        with open(snp_file_path, 'r') as snp_file:
            # Read each line, strip whitespace, and store in a set (to avoid duplicates)
            snps_to_analyze = {line.strip() for line in snp_file if line.strip()}
        st.write(f"- Loaded {len(snps_to_analyze)} SNPs to analyze.")

        # Load the haploblock data from the haploblock composition file
        st.write("**Loading Haploblock data...**")
        haploblock_data = []
        with open(haploblock_file_path, 'r') as haplo_file:
            for line_number, line in enumerate(haplo_file, start=1):
                line = line.strip()

                # Skip header
                if line_number == 1:
                    continue

                # Ignore empty lines
                if not line:
                    continue

                # Corrected splitting logic
                parts = line.split('\t')
                if len(parts) < 3:
                    st.write(f"Warning: Line {line_number} is malformed and will be skipped.")
                    continue

                haploblock_id = parts[0]  # Haploblock ID
                core_snp = parts[1]  # Core SNP
                snps_in_haploblock = parts[2].split()  # List of SNPs

                # Append data correctly
                haploblock_data.append((haploblock_id, core_snp, snps_in_haploblock))

        st.write(f"- Loaded {len(haploblock_data)} haploblocks.")

        # Process each SNP to determine its associated haploblocks
        st.write("**Processing SNPs and identifying associated haploblocks...**")
        with open(output_file_path, 'w') as output_file:
            for snp in snps_to_analyze:
                output_file.write(f"SNP: {snp}\n")

                # Find haploblocks where this SNP is the core SNP
                core_snp_haploblocks = [haploblock_id for haploblock_id, core_snp, snps in haploblock_data if snp == core_snp]
                if core_snp_haploblocks:
                    output_file.write("\tHaploblocks where this SNP is a core SNP:\n")
                    for haploblock_id in core_snp_haploblocks:
                        output_file.write(f"\t\t{haploblock_id}\n")
                else:
                    output_file.write("\tNo haploblocks where this SNP is a core SNP.\n")

                # Find haploblocks where this SNP is present (anywhere in the haploblock)
                found_haploblocks = [haploblock_id for haploblock_id, core_snp, snps in haploblock_data if snp in snps]
                if found_haploblocks:
                    for haploblock_id in found_haploblocks:
                        output_file.write(f"\tHaploblock: {haploblock_id}\n")
                else:
                    output_file.write("\tNo haploblocks found for this SNP.\n")

        # Check if the output file was successfully created and contains data
        if os.path.exists(output_file_path) and os.path.getsize(output_file_path) > 0:
            st.write(f"✅ Results successfully saved in {output_file_path}")
            return True
        else:
            st.write("❌ Output file was not generated correctly.")
            return False

    except Exception as e:
        # Log and display any errors encountered during execution
        logging.error(f"An error occurred: {e}")
        st.write(f"❌ An error occurred: {e}")
        return False

    finally:
        # Restore standard output
        sys.stdout = sys.__stdout__

# TODO : annovar.
def list_snps_carrying_given_snps(carrier_percentage_cut: float, carrier_percentage_mode: Literal['rough', 'exact'], work_directory: str,
    source_vcf_filepath: str, file_with_snps_to_analyze: str, region_to_extract: str | None = None, maf_threshold: float | None = None):
    """
    Identifies SNPs that are present in a given genomic region and determines which of these are carryied by specific SNPs of interest (from a list).

    Parameters:
    - carrier_percentage_cut (float): The minimum percentage of individuals required to carry the SNP for it to be considered significant.
    - carrier_percentage_mode (Literal['rough', 'exact']): Specifies whether carrier percentage is calculated using rough or exact method.
    - work_directory (str): Path to the directory where results will be saved.
    - source_vcf_filepath (str): Path to the VCF file containing the full set of genomic variants.
    - file_with_snps_to_analyze (str): Path to a file containing a list of SNPs to analyze (either a VCF file or a list of rsIDs).
    - region_to_extract (str | None): Path to a file containing a list of rsIDs specifying a genomic region to extract.
                                      If None, all SNPs in the source VCF are considered.
    - maf_threshold (float | None): The minimum minor allele frequency required for a SNP to be added in the analysis.
                                      If None, all SNPs are considered.

    Outputs:
    - Extracted SNPs and analysis results are saved in the specified work directory.
    """

    # Redirect stdout to capture print statements
    output_capture = StringIO()
    sys.stdout = output_capture

    try:
        # Initialize logging for this execution
        log_filepath = utils.generate_log_filename()  # Generate a dynamic log filename
        utils.setup_logging(log_filepath)  # Initialize logging
        logging.info(f"Output will be saved to {log_filepath}\n")
        st.write(f"Output will be saved to {log_filepath}\n")

        # Ensure the work directory is valid
        logging.info("Checking work directory…")
        st.write("**Checking work directory…**")
        if not work_directory.strip():
            work_directory = os.curdir
            logging.info(f"No directory provided: defaulting to current directory {work_directory}.")
            st.write(f"- No directory provided: defaulting to current directory {work_directory}.")
        elif not os.path.exists(work_directory):
            os.makedirs(work_directory)
            logging.info(f"Work directory created: {work_directory}.")
            st.write(f"- Work directory created: {work_directory}.")
        else:
            logging.info(f"Using existing work directory: {work_directory}.")
            st.write(f"- Using existing work directory: {work_directory}.")

        # Validate the source VCF file
        logging.info("Checking source VCF file…")
        st.write("**Checking source VCF file…**")
        if not source_vcf_filepath.strip():
            st.write("- No source VCF file provided.")
            raise ValueError("- No source VCF file provided.")
        elif not files.is_vcf_file(source_vcf_filepath):
            st.write("- Invalid VCF file format provided.")
            raise ValueError("- Invalid VCF file format provided.")
        else:
            logging.info("- Source VCF file verified.")
            st.write("✅ Source VCF file verified.")

        # Extract SNPs from a specific region if provided, otherwise analyze all SNPs
        logging.info("Checking region to extract…")
        st.write("**Checking region to extract…**")
        if not region_to_extract:
            logging.info("No specific region provided. Analyzing all SNPs in the source VCF.")
            st.write("- No specific region provided. Analyzing all SNPs in the source VCF.")
            path_to_vcf_with_all_work_snps = source_vcf_filepath
        elif not files.is_list_of_rsIds(region_to_extract):
            logging.warning("Invalid region format. Defaulting to analyzing all SNPs.")
            st.write("- Invalid region format. Defaulting to analyzing all SNPs.")
            path_to_vcf_with_all_work_snps = source_vcf_filepath
        else:
            logging.info("Extracting SNPs from the specified region in the source VCF…")
            st.write("- Extracting SNPs from the specified region in the source VCF…")
            path_to_vcf_with_all_work_snps = files.extract_snps_from_vcf(
                source_vcf_filepath, region_to_extract, f"{work_directory}/extracted_region.vcf"
            )
            logging.info("Region extraction complete.")
            st.write("**Region extraction complete.**")

        # Validate the file containing significant SNPs
        logging.info("Checking SNP list to analyze…")
        st.write("**Checking SNP list to analyze…**")
        if not file_with_snps_to_analyze.strip():
            raise ValueError("- No file provided containing the list of SNPs to analyze.")

        if files.is_vcf_file(file_with_snps_to_analyze):
            path_to_vcf_with_snps_to_analyze = file_with_snps_to_analyze
            logging.info("- SNP list is a valid VCF file.")
            st.write("- SNP list is a valid VCF file.")
        elif files.is_list_of_rsIds(file_with_snps_to_analyze):
            path_to_vcf_with_snps_to_analyze = files.extract_snps_from_vcf(
                source_vcf_filepath, file_with_snps_to_analyze, f"{work_directory}/snps_to_analyze.vcf"
            )
            logging.info("- SNP list extracted successfully.")
            st.write("✅ SNP list extracted successfully.")
        else:
            raise ValueError("❌ The list of SNPs to analyze must be either a VCF file or a text file containing rsIDs.")

        # Process the VCF files to retrieve SNP data
        logging.info("Extracting SNP data from VCF files…")
        st.write("**Extracting SNP data from VCF files…**")
        all_work_snps_data = files.retrieve_snp_data_from_vcf(path_to_vcf_with_all_work_snps, maf_threshold)
        snps_to_analyze_data = files.retrieve_snp_data_from_vcf(path_to_vcf_with_snps_to_analyze, maf_threshold)
        logging.info("SNP data extraction complete.")
        st.write("- SNP data extraction complete.")

        # Identify SNPs that are carried by the given SNPs
        for snp_to_analyze in snps_to_analyze_data:
            snp_directory = os.path.join(work_directory, snp_to_analyze.get_snp().get_rs_id())
            os.makedirs(snp_directory, exist_ok=True)

            snps_carriers, path_to_list_of_snps_carriers = genetic_computations.extract_snps_carried_by_reference(
                snp_to_analyze, all_work_snps_data, snp_directory, carrier_percentage_cut, carrier_percentage_mode
            )

            # Log results
            if not snps_carriers:
                logging.info(f"No SNPs found carrying SNP {snp_to_analyze.get_snp().get_rs_id()}.")
                st.write(f"❌ No SNPs found carrying SNP {snp_to_analyze.get_snp().get_rs_id()}.")
            else:
                logging.info(f"SNPs carrying SNP {snp_to_analyze.get_snp().get_rs_id()} saved at {path_to_list_of_snps_carriers}.")
                st.write(f"SNPs carrying SNP {snp_to_analyze.get_snp().get_rs_id()} saved at {path_to_list_of_snps_carriers}.")
        return True
    
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        st.write(f"An error occurred: {e}")
        return False
    
    finally:
        # Restore standard output
        sys.stdout = sys.__stdout__

def find_haploblocks_in_a_region(r_square_cut: float, d_prime_cut: float, maf_percentage_cut: float | None, carrier_percentage_cut: float | None, carrier_percentage_mode: Literal['rough', 'exact'] | None,
    REGION_SIZE: int | None, MAX_BP_SIZE_FOR_HAPLOBLOCK: int | None, MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC: int | None, work_directory: str,
    source_vcf_filepath: str, specific_snps_file: str | None, run_mode: str, ld_parameter: str, all_hap: str | None, add_snp: bool | None, maf_threshold: float | None,
    extend_mode: bool| None, extend_threshold: float | None, pruning_mode: bool |None, overlap_threshold: float | None):
    """
    Identifies haploblocks within a genomic region.

    Parameters:
    - r_square_cut (float): Minimum r² threshold to add the SNP to the haploblocks.
    - d_prime_cut (float): Minimum D' threshold to add the SNP to the haploblocks.
    - maf_percentage_cut (float | None): MAF percentage cutoff to add the SNP to the haploblocks.
    - carrier_percentage_cut (float | None): Minimum percentage of individuals carrying the haploblock.
    - carrier_percentage_mode (Literal['rough', 'exact'] | None): Defines how carrier percentage is calculated.
    - REGION_SIZE (int | None): Maximum base-pair (bp) size of the genomic region.
    - MAX_BP_SIZE_FOR_HAPLOBLOCK (int | None): Maximum allowed bp size for a haploblock (also the overlap between regions).
    - MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC (int | None): Maximum allowed empty SNP gap within a haploblock.
    - maf_threshold (float | None): Minimum allele frequency threshold for SNP selection. If None, all SNPs are considered.
    - work_directory (str): Directory for storing output files.
    - source_vcf_filepath (str): Path to the source VCF file containing SNP data.
    - specific_snps_file (str | None): Path to a file containing specific SNPs to analyze.
    - run_mode (str): Defines the mode of execution.
    - ld_parameter (str): Specifies the LD measure used (e.g., r² or D').
    - all_hap (str | None): Additional parameter for haploblock computation.
    - add_snp (bool | None): Whether to add SNPs dynamically.

    Outputs:
    - Generates multiple output files and histograms and returns their paths or False if failed.
    """

    # Redirect stdout to capture print statements
    output_capture = StringIO()
    sys.stdout = output_capture

    try:
        # Logging setup
        log_filepath = utils.generate_log_filename()

        # Avoid duplicate logging setup
        if not logging.getLogger().hasHandlers():
            utils.setup_logging(log_filepath)

        logging.info(f"Logging started. All output will be saved to {log_filepath}\n")
        st.write(f"All output will be saved to {log_filepath}\n")

        # Log the parameters used in the analysis
        logging.info(f"Using parameters: r_square_cut={r_square_cut}, d_prime_cut={d_prime_cut}, maf_percentage_cut={maf_percentage_cut}, "
                     f"carrier_percentage_cut={carrier_percentage_cut}, carrier_percentage_mode={carrier_percentage_mode}, "
                     f"MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC={MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC}, "
                     f"MAX_BP_SIZE_FOR_HAPLOBLOC={MAX_BP_SIZE_FOR_HAPLOBLOCK}, REGION_SIZE={REGION_SIZE}, maf_threshold={maf_threshold}")

        # Ensure that the working directory exists, create if necessary
        if not os.path.exists(work_directory):
            os.makedirs(work_directory)
            logging.info(f"Created work directory: {work_directory}\n")
            st.write(f"Created work directory: {work_directory}\n")
        else:
            logging.info(f"Using existing work directory: {work_directory}\n")
            st.write(f"Using existing work directory: {work_directory}\n")

        # Validate existence of the source VCF file
        if not os.path.exists(source_vcf_filepath):
            logging.error(f"❌ CRITICAL ERROR: The VCF file {source_vcf_filepath} does not exist. Exiting.")
            st.write(f"❌ CRITICAL ERROR: The VCF file {source_vcf_filepath} does not exist. Exiting.")
            return False

        # Check if the file has a valid VCF format
        if not files.is_vcf_file(source_vcf_filepath):
            logging.error(f"❌ Invalid VCF file format: {source_vcf_filepath}. Exiting.")
            st.write(f"❌ Invalid VCF file format: {source_vcf_filepath}. Exiting.")
            return False

        logging.info("VCF file validated successfully.")
        st.write("VCF file validated successfully.")

        # Set path to the VCF with all SNP data
        path_to_vcf_with_all_work_snps = source_vcf_filepath
        specific_snps = None

        # Load specific SNPs from file (if provided)
        if specific_snps_file:
            if not os.path.exists(specific_snps_file):
                logging.error(f"❌ Error: The specific SNPs file {specific_snps_file} does not exist. Exiting.")
                st.write(f"❌ Error: The specific SNPs file {specific_snps_file} does not exist. Exiting.")
                return False
            with open(specific_snps_file, 'r') as f:
                specific_snps = [line.strip() for line in f.readlines() if line.strip()]
            logging.info(f"Loaded SNP list (first 5 SNPs): {specific_snps[:5]}\n")
            st.write(f"Loaded SNP list (first 5 SNPs): {specific_snps[:5]}\n")

        # Extract SNP data from VCF
        logging.info("Processing VCF and extracting SNP data...\n")
        st.write("Processing VCF and extracting SNP data...\n")

        try:
            all_work_snps_data = files.retrieve_snp_data_from_vcf(path_to_vcf_with_all_work_snps, maf_threshold, specific_snps)

            if not all_work_snps_data:
                logging.error("❌ No SNP data extracted. Exiting.")
                st.write("❌ No SNP data extracted. Exiting.")
                return False

            logging.info("✅ Successfully extracted SNP data.")
            st.write("✅ Successfully extracted SNP data.")

        except Exception as e:
            logging.error(f"❌ Error processing VCF file {path_to_vcf_with_all_work_snps}: {e}")
            st.write(f"❌ Error processing VCF file {path_to_vcf_with_all_work_snps}: {e}")
            return False

        # Construct haploblocks based on LD and filtering criteria
        logging.info("Building haploblocks...")
        st.write("Building haploblocks...")
        try:
            haploblocks = genetic_computations.build_haploblocks_per_region(
                all_work_snps_data, REGION_SIZE, MAX_BP_SIZE_FOR_HAPLOBLOCK, MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC,
                r_square_cut, d_prime_cut, maf_percentage_cut, carrier_percentage_cut, carrier_percentage_mode,
                ld_parameter, run_mode, specific_snps, all_hap, add_snp,
                extend_mode, extend_threshold, pruning_mode, overlap_threshold
            )
            if not haploblocks:
                logging.warning("⚠ No haploblocks generated. End of the program.")
                st.write("⚠ No haploblocks generated. End of the program.")
                return False
            logging.info(f"✅ Successfully built {len(haploblocks)} haploblocks.")
            st.write(f"✅ Successfully built {len(haploblocks)} haploblocks.")
        except Exception as e:
            logging.error(f"❌ Error generating haploblocks: {e}")
            st.write(f"❌ Error generating haploblocks: {e}")
            return False

        # Create output file containing haploblock information
        additional_info = ""
        chromosome = haploblocks[0].get_all_snps()[0].get_snp().get_chromosome()
        if chromosome:
            additional_info += f"chr_{chromosome}"

        # Adjust additional_info based on parameters
        if ld_parameter == "1":
            additional_info += ".withLD"
        elif ld_parameter == "2":
            additional_info += ".noLD"

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

        logging.info('Generating the Haploblock information files...')
        st.write('Generating the Haploblock information files...')
        path_to_haploblock_data = f"{work_directory}/HAPLOBLOCKS_data.{additional_info}.txt"
        path_to_haploblock_comp = f"{work_directory}/HAPLOBLOCKS_comp.{additional_info}.txt"
        path_to_haploblock_seq = f"{work_directory}/HAPLOBLOCKS_seq.{additional_info}.txt"
        path_to_snp_maf = f"{work_directory}/SNP_MAF.{additional_info}.txt"

        files.generate_haploblock_files(all_work_snps_data, haploblocks, path_to_vcf_with_all_work_snps, path_to_haploblock_data, path_to_haploblock_comp, path_to_haploblock_seq, path_to_snp_maf)
        logging.info(f"✅ Haploblock information files generated: {path_to_haploblock_data}, {path_to_haploblock_comp}, {path_to_haploblock_seq},{path_to_snp_maf}")
        st.write(f"✅ Haploblock information files generated")
        st.write(f"- Haploblocks data : {path_to_haploblock_data}")
        st.write(f"- Haploblocks composition : {path_to_haploblock_comp}")
        st.write(f"- Haploblocks sequence : {path_to_haploblock_seq}")
        st.write(f"- SNPs informations : {path_to_snp_maf}")

        logging.info('Generating the histograms...')
        st.write('Generating the histograms...')
        path_to_comp_hist = f"{work_directory}/Histogram_Haploblocks-SNP.{additional_info}.pdf"
        path_to_size_hist = f"{work_directory}/Histogram_Haploblocks-Size.{additional_info}.pdf"
        snp_success, size_success = plotting.generate_haploblock_histograms(path_to_haploblock_data, path_to_comp_hist, path_to_size_hist)
        logging.info(f"✅ Histograms generated: {path_to_comp_hist}, {path_to_size_hist}")
        st.write(f"✅ Histograms generated")
        st.write(f"- Haploblocks SNPs size Distribution : {path_to_comp_hist}")
        st.write(f"- Haploblocks Bp size Distribution : {path_to_size_hist}")

        if snp_success and size_success:
            return path_to_comp_hist, path_to_size_hist
        else:
            return False

    except Exception as e:
        logging.critical(f"Unexpected error: {e}")
        logging.debug(traceback.format_exc())
        st.write(f"Unexpected error: {e}")
        return False

    finally:
        # Retrieve the captured output
        output = output_capture.getvalue()
        # Display the output in Streamlit
        st.text(output)
        # Restore stdout
        sys.stdout = sys.__stdout__
        for handler in logging.getLogger().handlers:
            handler.flush()

def print_haploblocks(work_directory: str, haploblock_file: str, snp_file: str, min_size_snp: int | None = None, max_size_snp: int | None = None,
                      min_size_bp: int | None = None, max_size_bp: int | None = None, bp_range_start: int | None = None,
                      bp_range_end: int | None = None, maf_threshold: float | None = None):
    """
    Process haploblock data, apply filtering criteria, and generate plots.
    """

    # Redirect stdout to capture print statements
    output_capture = StringIO()
    sys.stdout = output_capture

    try:
        logging.info("Checking the directory...")
        st.write("**Checking the directory...**")

        if not work_directory:
            work_directory = os.curdir
            logging.info(f"- Default directory set to: {work_directory}.")
        elif not os.path.exists(work_directory):
            os.makedirs(work_directory)
            logging.info(f"- Directory created: {work_directory}.")
        else:
            logging.info(f"- Selected directory: {work_directory}.")

        st.write(f"**Checking input files...**")
        if not haploblock_file or not files.is_haploblock_composition(haploblock_file):
            raise ValueError("- Invalid haploblock composition file.")
        if not snp_file or not files.is_SNP_MAF(snp_file):
            raise ValueError("- Invalid SNP-MAF file format.")

        st.write("**Loading haploblock composition...**")
        haploblock_composition = files.read_haploblock_composition(haploblock_file)

        st.write("**Loading SNP-MAF file...**")
        snp_maf_dict = files.read_snp_maf_file(snp_file)

        st.write("✅ Input files loaded.")

        # Filtering haploblocks
        filtered_haploblock_composition = {}
        for haploblock_id, (core_snp, snps) in haploblock_composition.items():
            if not isinstance(haploblock_id, str):
                raise TypeError(f"Haploblock ID must be a string, got {type(haploblock_id)}")

            snp_positions = [snp_maf_dict[snp][1] for snp in snps if snp in snp_maf_dict]
            if not snp_positions:
                st.write(f"Warning: No valid SNP positions found for haploblock {haploblock_id}")
                continue

            bloc_start_bp = min(snp_positions)
            bloc_end_bp = max(snp_positions)
            bloc_size_bp = bloc_end_bp - bloc_start_bp
            bloc_size_snp = len(snps)

            if min_size_snp and bloc_size_snp < int(min_size_snp):
                continue
            if max_size_snp and bloc_size_snp > int(max_size_snp):
                continue
            if min_size_bp and bloc_size_bp < int(min_size_bp):
                continue
            if max_size_bp and bloc_size_bp > int(max_size_bp):
                continue
            if bp_range_start and bloc_end_bp < int(bp_range_start):
                continue
            if bp_range_end and bloc_start_bp > int(bp_range_end):
                continue

            # Convert snps list to a tuple (to ensure it's hashable)
            filtered_haploblock_composition[haploblock_id] = (core_snp, tuple(snps))  

        # Check if there are any haploblocks left after filtering
        if not filtered_haploblock_composition:
            raise ValueError("No haploblocks left after applying filters.")

        # Plotting the haploblocks
        total_haploblocks = len(filtered_haploblock_composition)
        st.write(f"Total haploblocks after filtering: {total_haploblocks}")
        max_haploblocks_per_plot = 50
        num_parts = ceil(total_haploblocks / max_haploblocks_per_plot)

        output_file = None  # Initialize output_file to avoid reference before assignment

        for part in range(num_parts):
            start_idx = part * max_haploblocks_per_plot
            end_idx = min((part + 1) * max_haploblocks_per_plot, total_haploblocks)
            chunk_haploblock_composition = {haploblock_id: filtered_haploblock_composition[haploblock_id]
                                            for haploblock_id in list(filtered_haploblock_composition.keys())[start_idx:end_idx]}

            output_file = os.path.join(work_directory, f"haploblocks_part_{part + 1}.png")

            try:
                plot = plotting.plot_haploblocks(output_file, snp_maf_dict, chunk_haploblock_composition, maf_threshold)
                if plot:
                    st.write(f"**Graph generated in: {output_file}**")
                else:
                    st.write(f"**Error generating graph for haploblocks.**")
                    break
            except Exception as e:
                st.write(f"**Error plotting haploblocks: {e}**")
                break

        return os.path.exists(output_file) if output_file else False

    except Exception as e:
        st.write(f"**An error occurred: {e}**")

    finally:
        output = output_capture.getvalue()
        st.write(output)
        sys.stdout = sys.__stdout__
        for handler in logging.getLogger().handlers:
            handler.flush()



