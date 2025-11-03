import subprocess
import os
import gzip
import re
from genetic_structures import SNP, SNPInformations
import utils
import logging
import subprocess
import multiprocessing
from typing import List, Optional

# Path to local version of bcftools.
BCFTOOLS = "bcftools"

# Patterns for file containing SNP or haploblock data.
RS_PATTERN = r'rs\d+|\.'
CHR_PATTERN = r'chr\d+:\d+:[A-Z]:[A-Z]'
HLA_PATTERN = r'HLA_[A-Z0-9_]+' 
CHROM_POS_REF_ALT_PATTERN = r'(?:[1-9]|1[0-9]|2[0-2]|X|Y):\d+:[A-Z]+:[A-Z]+'

SNP_ID_PATTERN = re.compile(f'^(?:{RS_PATTERN}|{CHR_PATTERN}|{HLA_PATTERN}|{CHROM_POS_REF_ALT_PATTERN})$')
HAPLOBLOCK_PATTERN = [int, SNP_ID_PATTERN, int, SNP_ID_PATTERN, int, SNP_ID_PATTERN, int, int]
HAPLOBLOCK_HEADER_PATTERN = ['haploblock_index', 'core_snp_id', 'core_snp_bp', 'start_id', 'start_bp', 'end_id', 'end_bp', 'size_in_snps']

def is_vcf_file(filename: str) -> bool:
    """
    Checks if a file is a VCF.
    """
    # Check if the file exists
    if not os.path.isfile(filename):
        logging.error(f"File not found: {filename}")
        return False

    # Check if the file has a .vcf or .vcf.gz extension
    if not (filename.lower().endswith('.vcf') or filename.lower().endswith('.vcf.gz')):
        return False

    # Check the content of the file
    try:
        open_func = gzip.open if filename.lower().endswith('.vcf.gz') else open
        mode = 'rt' if filename.lower().endswith('.vcf.gz') else 'r'

        with open_func(filename, mode, encoding='utf-8', errors='ignore') as file:
            # Read the first few lines to check for the VCF header
            for line in file:
                if line.startswith('##fileformat=VCF'):
                    return True
                if not line.startswith('##'):
                    # Stop checking if we've passed the header lines
                    break
    except Exception as e:
        logging.error(f"Error reading file: {e}")
        return False

    return False

def is_list_of_rsIds(filename : str) -> bool :
    """
    Checks if a given text file contains a list of rsIds that can be used as a vcf/bcf tools parameter.
    """
    with open(filename, 'r', encoding='utf-8', errors='ignore') as file:
        # Read the first line.
        first_line = file.readline().strip()
        # Check if the first line matches the pattern of a SNP rsId.
        if not SNP_ID_PATTERN.match(first_line):
            # If the first line doesn't match, treat it as a header and skip it.
            # check the rest of the file
            lines = file
        else:
            # If the first line matches, process it as a data line.
            lines = [first_line] + file.readlines()

        for line in lines:
            line = line.strip()
            if not SNP_ID_PATTERN.match(line):
                return False
    return True

def is_SNP_MAF(filename):
    """
    Check if the SNP-MAF file is in the correct format: ID MAF position.
    The ID format should follow alphanumeric with colons, and MAF should be a decimal in [0, 1].
    """
    snp_maf_pattern = re.compile(r'^\d+:\d+:[ACGT]+:(?:[ACGT]+|<[^<>]+>) \d\.\d{2} \d+$')  # Line pattern

    try:
        with open(filename, 'r') as f:
            for line_number, line in enumerate(f, start=1):
                line = line.strip()
                if not snp_maf_pattern.match(line):
                    print(f"Line {line_number}: Invalid format '{line}'")
                    return False  # Return False if any line doesn't match the format
                
                # Further check that MAF is between 0 and 1
                snp_id, maf, pos = line.split(' ')
                if not (0 <= float(maf) <= 1):
                    print(f"Line {line_number}: Invalid MAF value '{maf}', should be between 0 and 1.")
                    return False
        return True  # All lines match the format
    except Exception as e:
        print(f"Error when validating the file: {e}")
        return False

def read_snp_maf_file(maf_file_path):
    """
    Reads the SNP MAF file and returns a dictionary with SNP ID as keys and a tuple of (MAF, position) as values.
    The file format: SNP_ID MAF position (separated by one space).
    """
    snp_maf_dict = {}
    try:
        with open(maf_file_path, 'r') as f:
            for line_number, line in enumerate(f, start=1):
                line = line.strip()
                if " " in line:  # Ensure the format is ID MAF position with one space
                    try:
                        snp_id, maf, pos = line.split(' ')  # Split by a single space
                        snp_maf_dict[snp_id] = (float(maf), int(pos))  # Store MAF as float and position as int
                    except ValueError:
                        print(f"Line {line_number}: Format error in '{line}'")
    except Exception as e:
        print(f"Error when reading the file: {e}")
    return snp_maf_dict

def is_haploblock_composition(file_path: str) -> bool:
    """
    Validates if the haploblock composition file is in the expected format:

    Expected format:
    Haploblock ID\tCore SNP\tHaploblock composition
    Where:
    - Haploblock ID matches "Haploblock X" (X is an integer).
    - Core SNP and all SNPs in the composition follow CHR:POS:REF:ALT format.

    Args:
        file_path (str): Path to the file to validate.

    Returns:
        bool: True if the file format is valid, False otherwise.
    """
    haploblock_pattern = re.compile(r'^Haploblock \d+$')  # Pattern for haploblock identifier
    snp_pattern = re.compile(r'^\d+:\d+:[ACGT]+:(?:[ACGT]+|<[^<>]+>)$')  # SNP format (CHR:POS:REF:ALT)

    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            for line_number, line in enumerate(file, start=1):
                line = line.strip()

                # Skip the header
                if line_number == 1:
                    continue

                # Ignore empty lines
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) != 3:
                    print(f"Line {line_number}: Invalid format (missing Core SNP or Haploblock composition).")
                    return False

                haploblock_id, core_snp, haploblock_composition = parts

                # Validate Haploblock ID
                if not haploblock_pattern.match(haploblock_id):
                    print(f"Line {line_number}: Invalid Haploblock ID '{haploblock_id}'.")
                    return False

                # Validate Core SNP
                if not snp_pattern.match(core_snp):
                    print(f"Line {line_number}: Invalid Core SNP format '{core_snp}'.")
                    return False

                # Validate each SNP in the Haploblock composition
                haploblock_snps = haploblock_composition.split()
                for snp in haploblock_snps:
                    if not snp_pattern.match(snp):
                        print(f"Line {line_number}: Invalid SNP format in Haploblock composition '{snp}'.")
                        return False

        return True

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' does not exist.")
        return False
    except IOError as e:
        print(f"Error reading file '{file_path}': {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False

def read_haploblock_composition(file_path):
    """
    Reads the haploblock composition file and returns a dictionary:
    {Haploblock_ID: (Core SNP, [Composition SNPs])}

    Args:
        file_path (str): Path to the haploblock composition file.

    Returns:
        dict: Dictionary mapping each Haploblock_ID to its (Core SNP, List of SNPs).
    """
    haploblock_composition = {}

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_number, line in enumerate(f, start=1):
                line = line.strip()

                # Skip header line
                if line_number == 1:
                    continue

                # Skip empty lines
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) != 3:
                    print(f"Warning: Line {line_number} is malformed and will be skipped.")
                    continue

                haploblock_id = parts[0]
                core_snp = parts[1]
                snps = parts[2].split()

                haploblock_composition[haploblock_id] = (core_snp, snps)

    except Exception as e:
        print(f"Error when reading the composition file: {e}")

    return haploblock_composition

def extract_snps_from_vcf(vcf_input : str, snp_list : str, output_filepath : str) -> str | None :
    """
    Extract a list of SNPs from a vcf file into a new vcf file, using bcftools.
    """
    bcftools_command = [BCFTOOLS, 'view', '-i', f'ID=@{snp_list}', '-Oz', '-o', output_filepath, vcf_input]
    try:
        subprocess.run(bcftools_command, check=True)
        return output_filepath
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred while running bcftools : {e}")
        return None
    
def replace_dot_with_pattern(rs_id: str, snp: SNP) -> SNP:
    if rs_id == '.' or rs_id == None:
        new_rs_id = f"{snp.get_chromosome()}:{snp.get_position()}:{snp.get_ref_allele()}:{snp.get_alt_allele()}"
        return SNP(new_rs_id, snp.get_chromosome(), snp.get_position(), snp.get_ref_allele(), snp.get_alt_allele())
    else:
        return SNP(rs_id, snp.get_chromosome(), snp.get_position(), snp.get_ref_allele(), snp.get_alt_allele())

def process_snp_line(line, snp_data, maf_threshold: Optional[float] = None, snp_list : Optional[list] = None):
    """
    Process a single SNP line and add the result to the shared data structure.
    Handles parsing and calculations for SNP information.
    """
    if line.startswith('#') or not line.strip():
        return  # Skip header lines or empty lines

    fields = line.strip().split('\t')
    if len(fields) < 10:  # Ensure sufficient fields, genotypes start at column 10
        logging.warning(f"Skipping line due to insufficient columns: {line.strip()}")
        return

    rs_id, chrom, pos, ref, alt = fields[:5]
    genotypes = fields[5:]  # Genotypes are in columns starting from the 6th column onward

    try:
        pos = int(pos)  # Ensure position is an integer
    except ValueError:
        logging.warning(f"Skipping line with invalid position: {line.strip()}")
        return

    alt_count = 0
    individuals_genotyped = []
    phasing = []
    alt_carriers = []
    alt_carriers_1 = []
    alt_carriers_2 = []
    ref_carriers = []
    ref_carriers_1 = []
    ref_carriers_2 = []

    for i, gt in enumerate(genotypes):
        if gt == '.' or not gt:
            continue
        if '|' in gt:
            alleles = gt.split('|')
            phasing.append(True)
        else:
            alleles = gt.split('/')
            phasing.append(False)

        if len(alleles) == 2 and all(allele.isdigit() for allele in alleles):
            if int(alleles[0]) == 1:
                alt_carriers.append(i)
                if phasing[-1]:
                    alt_carriers_1.append(i)
            else:
                ref_carriers.append(i)
                if phasing[-1]:
                    ref_carriers_1.append(i)
            if int(alleles[1]) == 1:
                alt_carriers.append(i)
                if phasing[-1]:
                    alt_carriers_2.append(i)
            else:
                ref_carriers.append(i)
                if phasing[-1]:
                    ref_carriers_2.append(i)
            count = sum(map(int, alleles))
            alt_count += count
            individuals_genotyped.append(i)

    snp_phased = all(phasing) if phasing else False
    if not snp_phased:
        alt_carriers_1, alt_carriers_2 = [], []
        ref_carriers_1, ref_carriers_2 = [], []

    maf = alt_count / (2 * len(individuals_genotyped)) if individuals_genotyped else 0
    if maf > 0.5:
        alt_carriers, ref_carriers = ref_carriers, alt_carriers
        alt_carriers_1, ref_carriers_1 = ref_carriers_1, alt_carriers_1
        alt_carriers_2, ref_carriers_2 = ref_carriers_2, alt_carriers_2
        maf = 1 - maf
        ref, alt = alt, ref

    # Create the SNP object
    snp = SNP(rs_id, str(chrom), pos, ref, alt)
    snp = replace_dot_with_pattern(rs_id, snp)

    if snp_list and snp.get_rs_id().strip() in [s.strip() for s in snp_list]:
        snp_data.append(SNPInformations(snp, maf, snp_phased, individuals_genotyped,
                                        alt_carriers, alt_carriers_1, alt_carriers_2,
                                        ref_carriers, ref_carriers_1, ref_carriers_2))
        return

    elif maf_threshold is None or maf >= maf_threshold:
        snp_data.append(SNPInformations(snp, maf, snp_phased, individuals_genotyped,
                                        alt_carriers, alt_carriers_1, alt_carriers_2,
                                        ref_carriers, ref_carriers_1, ref_carriers_2))

def process_snp_chunks(chunk: List[str], maf_threshold: Optional[float], snp_list: Optional[list]) -> List:
    """
    Process a chunk of SNP lines in a way that can be pickled
    
    Args:
        chunk (List[str]): Lines to process
        maf_threshold (Optional[float]): Minor Allele Frequency threshold
        snp_list (Optional[list]): List of specific SNPs to filter
    
    Returns:
        List: Processed SNP data
    """
    snp_data = []
    for line in chunk:
        process_snp_line(line, snp_data, maf_threshold, snp_list)
    
    return snp_data

def retrieve_snp_data_from_vcf(vcf_filepath: str, maf_threshold: Optional[float] = None , snp_list: Optional[list] = None) -> List:
    """
    Retrieve SNP data from VCF file using bcftools and multiprocessing
    
    Args:
        vcf_filepath (str): Path to VCF file
        maf_threshold (float): Minor Allele Frequency threshold
        snp_list (Optional[list]): List of specific SNPs to filter
    
    Returns:
        List: Processed SNP data
    """
    # Check if bcftools is available
    try:
        subprocess.run(['bcftools', '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        logging.error("❌ BCFtools is not installed. Please install it to proceed.")
        print("❌ BCFtools is not installed. Please install it to proceed.")
        return []

    # Construct bcftools command
    bcftools_command = ['bcftools', 'query', '-f', '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n', vcf_filepath]
    
    try:
        # Run bcftools and capture output
        result = subprocess.run(bcftools_command, capture_output=True, text=True, check=True)
        
        # Split output into lines
        lines = result.stdout.strip().split('\n')
        
        # Determine number of processes
        num_processes = max(1, multiprocessing.cpu_count() - 1)
        
        # Split lines into chunks of lines for processing
        chunk_size = max(1, len(lines) // num_processes)
        chunks = [lines[i:i+chunk_size] for i in range(0, len(lines), chunk_size)]
        
        # Use multiprocessing to process chunks
        with multiprocessing.Pool(processes=num_processes) as pool:
            # Map the process_snp_chunks function to the chunks
            results = pool.starmap(
                process_snp_chunks, 
                [(chunk, maf_threshold, snp_list) for chunk in chunks]
            )
        
        # Flatten the results
        snp_data = [item for sublist in results for item in sublist]
        
        logging.info(f"Extracted {len(snp_data)} SNP entries")        
        return snp_data

    except subprocess.CalledProcessError as e:
        logging.error(f"BCFtools error: {e.stderr}")
        print(f"BCFtools error: {e.stderr}")
        return []
    except Exception as e:
        logging.error(f"Error while processing VCF: {e}")
        print(f"Error while processing VCF: {e}")
        return []
 
def generate_haploblock_files(all_work_snps_data: list[SNPInformations], haploblocks, path_to_vcf_with_all_work_snps, path_to_haploblock_data,
    path_to_haploblock_comp, path_to_haploblock_seq ,path_to_snp_maf_file) -> None:

    def ensure_txt_extension(file_path: str) -> str:
        """Ensure the file path ends with `.txt`."""
        if file_path is None:
            raise ValueError("File path is None")
        return file_path if file_path.endswith(".txt") else file_path + ".txt"

    def format_snp_info(snp: SNP) -> str:
        """Format SNP information as a string."""
        return f"{snp.get_chromosome()}:{snp.get_position()}:{snp.get_ref_allele()}:{snp.get_alt_allele()}"

    # Ensure file paths have the correct extensions
    try:
        logging.info(f"Ensuring file extensions: path_to_haploblock_comp={path_to_haploblock_comp}, path_to_haploblock_data={path_to_haploblock_data}, path_to_snp_maf_file={path_to_snp_maf_file}")
        path_to_haploblock_comp = ensure_txt_extension(path_to_haploblock_comp)
        path_to_haploblock_data = ensure_txt_extension(path_to_haploblock_data)
        path_to_snp_maf_file = ensure_txt_extension(path_to_snp_maf_file)
    except ValueError as e:
        logging.error(f"Error ensuring file extensions: {e}")
        return

    try:
        # Generate haploblock composition
        with open(path_to_haploblock_comp, 'w', encoding='utf-8') as comp_file:
            # Write header
            comp_file.write("Haploblock ID\tCore SNP\tHaploblock composition\n")

            for i, haploblock in enumerate(haploblocks):
                # Get core SNP and format its information
                core_snp = haploblock.get_core_snp()
                core_snp_info = format_snp_info(core_snp.get_snp())

                # Get all SNPs (including the core SNP) sorted by position
                sorted_snps = sorted(haploblock.get_all_snps(), key=lambda snp_info: snp_info.get_snp().get_position())

                # Format all SNPs for output
                haploblock_composition = ' '.join(format_snp_info(snp.get_snp()) for snp in sorted_snps)

                # Write haploblock data to file
                comp_file.write(f"Haploblock {i+1}\t{core_snp_info}\t{haploblock_composition}\n")

        logging.info(f"Haploblock composition for region {path_to_vcf_with_all_work_snps} created in {path_to_haploblock_comp}, containing {len(haploblocks)} blocks.")

        # Generate haploblock data
        header = "haploblock_index\tcore_snp_id\tcore_snp_bp\tstart_id\tstart_bp\tend_id\tend_bp\tsize_in_snps"
        with open(path_to_haploblock_data, 'w') as data_file:
            data_file.write(header + '\n')
            for i, haploblock in enumerate(haploblocks, start=1):
                line = (
                    f"{i}\t{haploblock.get_core_snp().get_snp().get_rs_id()}\t"
                    f"{haploblock.get_core_snp().get_snp().get_position()}\t"
                    f"{haploblock.get_start_snp().get_snp().get_rs_id()}\t"
                    f"{haploblock.get_start_bp()}\t"
                    f"{haploblock.get_end_snp().get_snp().get_rs_id()}\t"
                    f"{haploblock.get_end_bp()}\t"
                    f"{haploblock.get_snp_size()}"
                )
                data_file.write(line + '\n')

        logging.info(f"Haploblock data for region {path_to_vcf_with_all_work_snps} created in {path_to_haploblock_data} containing {len(haploblocks)} blocks.")

        # Generate MAF file for all SNPs
        with open(path_to_snp_maf_file, 'w', encoding='utf-8') as maf_file:
            for snp_info in all_work_snps_data:
                snp_id = format_snp_info(snp_info.get_snp())
                maf = f"{snp_info.get_maf():.2f}"  # Format MAF to 2 decimal places
                pos = f"{snp_info.get_snp().get_position()}"
                maf_file.write(f"{snp_id} {maf} {pos}\n")
        logging.info(f"SNP data for region {path_to_vcf_with_all_work_snps} created in {path_to_snp_maf_file}")

        # Generate haploblock sequence file
        with open(path_to_haploblock_seq, 'w', encoding='utf-8') as hap_seq_file:
            for i, haploblock in enumerate(haploblocks):
                # Sort SNPs in the haploblock by position
                sorted_snps = sorted(haploblock.get_all_snps(), key=lambda snp_info: snp_info.get_snp().get_position())
                alt_alleles = [snp_info.get_snp().get_alt_allele() for snp_info in sorted_snps]
                hap_seq_file.write(f"Haploblock {i+1}\n")
                hap_seq_file.write(" ".join(alt_alleles) + "\n\n")

    except IOError as e:
        logging.error(f"An I/O error occurred: {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")