import streamlit as st
import logging
import traceback
import sys
from io import StringIO

from main_script import *
from contextlib import redirect_stdout

# Set page configuration and title
st.set_page_config(
    page_title="HaploExplore",
    page_icon="🧬",
    layout="wide"
)

class StreamlitLogHandler(logging.Handler):
    def __init__(self, log_queue):
        super().__init__()
        self.log_queue = log_queue
        self.formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    def emit(self, record):
        try:
            msg = self.format(record)
            self.log_queue.put(msg)
        except Exception:
            self.handleError(record)


# Center the title at the top with dark blue color and larger size
st.markdown(
    """
    <style>
    .title {
        font-size: 50px; /* Slightly larger font size */
        font-weight: bold;
        text-align: center;
        margin-top: -50px; /* Adjust vertical alignment */
        color: #003366; /* Dark blue color */
    }
    </style>
    <div class="title">HaploExplore</div>
    """,
    unsafe_allow_html=True,
)

# Sidebar navigation for all 4 functions
st.sidebar.title("Menu")

option = st.sidebar.selectbox("Choose your analysis", [
    "Home Page",
    "Find Haploblocks in a Region",
    "Find Haploblocks for SNPs", 
    "List SNPs Carrying Given SNPs",
    "Plot Haploblocks"
])

footer = """
    <style>
    .sidebar-footer {
        position: fixed;
        bottom: 0;
        left: -41.5%;
        width: 100%;
        background-color: transparent;
        text-align: center;
        padding: 10px 0;
        font-size: 12px;
        color: gray;
    }
    .sidebar-footer a {
        text-decoration: none;
        color: #4F8BF9; /* Couleur du lien */
    }
    .sidebar-footer a:hover {
        text-decoration: underline;
    }
    </style>
    <div class="sidebar-footer">
        <b>© 2025 - GBCM Lab., CNAM</b><br>
        Developped by <a href="https://github.com/HaploExplore/HaploExplore_v1/" target="_blank">Samuel HIET</a>.
    </div>
"""
st.sidebar.markdown(footer, unsafe_allow_html=True)

def run_homepage():
    # Centered title and subtitle using custom HTML
    st.markdown(
        """
        <style>
        .subtitle {
            font-size: 30px; /* Smaller font for the subtitle */
            text-align: center;
            color: #003366; /* Same dark blue color */
            margin-top: 0; /* Reduce gap between title and subtitle */
        }
        .content {
            font-size: 16px;
            text-align: justify;
            margin: 20px auto;
            max-width: 800px; /* Restrict content width for better readability */
        }
        </style>
        <div class="subtitle">Home Page</div>
        """,
        unsafe_allow_html=True,
    )

    st.text("This tool is designed to analyze and visualize genetic data, specifically focusing on haploblocks and their associated SNPs (Single Nucleotide Polymorphisms). The tool offers functionalities to generates Minor alleles haploblocks from genotypic data, plot the distribution of haploblocks on a chromosome, visualize SNPs within a specific haploblock, and generate histograms illustrating SNP count distribution.")

    st.markdown("### **Usage**")
    st.text("To use one of the function use the sidebar and select the needed one.")

    st.markdown("## **Functions**")
    st.markdown("For more detail, go on the functions pages.")
    st.markdown("### **Data Analysis**")
    st.markdown("""
    - **`find_Haploblocks_in_a_region`**: Takes a genomic region (from a VCF file) and identifies haploblocks within it.  
      Generates output files and histograms for haploblock distribution.
    - **`find_Haploblocks_for_given_snps`**: Extracts Haploblocks associated with each SNP in the given list and writes the results to a file.
    - **`list_snps_carrying_given_snps`**: Given a genomic region and a list of SNPs, extracts all SNPs in the region that carry the same minor allele as each SNP in the list.
    """)

    st.markdown("### **Haploblock Visualization**")
    st.markdown("""
    - **`print_haploblock`**: Visualizes haploblocks on a chromosome. Each haploblock is represented as a rectangle, with its core SNP marked by a star.  
    """)

    st.markdown("### **Computational Modes**")
    st.markdown("""
        For the **`find_Haploblocks_in_a_region`**, HaploExplore offers three computational modes for different analysis needs:

        #### 1. Standard Mode
        Constructs haploblocks iteratively, with or without LD parameters. For each SNP, the algorithm will try to add this SNP to any existing haploblock (can be several), if the SNP can't be add to any existing haploblock it will generates its own haploblock and it will be the coreSNP of it.
        - **With LD :** SNPs are grouped based on MAF, LD thresholds (r² ≥ 0.1, D' ≥ 0.7), and a Carrier Percentage (CP) ≥ 75%.  
        - **Without LD :** SNPs are grouped based only on MAF and CP thresholds. This mode is faster but less precise.

        #### 2. ListSNP Mode
        Constructs haploblocks based on a predefined list of SNPs.  
        **Options:**  
        - **All Haploblocks:**          
            - **1** : Generates haploblocks for the SNPs in the list and for the SNPs that can't fit in any existing haploblock.
                
            - **2** : Only generates haploblocks for the SNPs in the list.

        - **Add SNP:** 
            - **1** : Adds the SNPs of the list in the haploblocks (not only as core SNPs).  
                
            - **2** : Does not add the SNPs of the list in the haploblocks (theses SNPs will be use only to generate the Haploblocks as tcoreag SNP).

        Supports both LD-based and non-LD-based grouping.

        #### 3. Exhaustive Mode
        Iteratively scans all SNPs and assigns each to its own haploblock initially. Refines haploblocks based on LD and carrier percentage criteria.  
        **Use Case :** Ideal for regions with complex LD patterns or high genomic diversity.  

        Supports both LD-based and non-LD-based grouping.
    """)
    return

# Function to run Find Haploblocks for Given SNPs
def run_find_Haploblocks_for_given_snps():
    st.header("Find Haploblocks for Given SNPs")

    st.markdown(
        """
        <style>
        .description-box {
            background-color: #eaf4fc; /* Light grey background */
            padding: 15px;
            border-radius: 10px;
            border: 1px solid #ddd; /* Optional: Add a border */
            font-size: 16px;
            line-height: 1.6;
        }
        </style>
        <div class="description-box">
            <p>Extracts Haploblocks associated with each SNP in the given list. This function won't generate the haploblocks.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("""### **Input file**""")
    st.markdown("""
    - A list of SNP IDs in chromosome-position-based SNP identifier format (Chromosome:Position:Position:RefAllele:AltAllele) in a text file, such as:
    """)
    
    snp_list = """SNP_ID1\nSNP_ID2\nSNP_ID3\n..."""
    st.code(snp_list, language="text")

    st.markdown("""
    - A composition file generated by `build_Haploblocks_per_region`.
    """)

    st.markdown("""### **Output file**""")
    st.markdown("""
    - Text file containing for each SNP the different haploblocks it can be found in.
    """)

    st.markdown("""### **Functional part**""")

    # Work directory input for output
    work_directory = st.text_input("Enter the output directory for results", value="")
    if not work_directory.strip():
        st.warning("Output directory is mandatory. Please specify a valid output directory.")
        return

    # User enters full paths for the SNP file and Haploblock composition file
    snp_file_path = st.text_input("Enter the full file path for the SNP list file", value="")
    haploblock_file_path = st.text_input("Enter the full file path for the Haploblock composition file", value="")

    # Check if both file paths are provided
    if not snp_file_path.strip() or not haploblock_file_path.strip():
        st.warning("Both SNP file path and Haploblock composition file path are mandatory.")
        return

    # Button to run the analysis
    if st.button("Run Analysis"):
        st.write("### **Run logs :**")
        # Create a placeholder for live output
        output_placeholder = st.empty()

        # Redirect stdout to capture print statements
        output_capture = StringIO()
        sys.stdout = output_capture

        try:
            with st.spinner('Searching Haploblocks for the given SNPs...'):
                # Call the function with user-provided paths
                success = find_haploblocks_for_given_snps(work_directory, snp_file_path, haploblock_file_path)
                if success:
                    st.success("✅ Analysis Complete!")
                else:
                    st.error("❌ Error: Analysis failed or returned incomplete results.")
        except Exception as e:
            st.error(f"An error occurred: {e}")

        # Retrieve the captured output
        output = output_capture.getvalue()
        # Display the output in Streamlit
        output_placeholder.text(output)
        # Restore stdout
        sys.stdout = sys.__stdout__

# Function to run List SNPs Carrying Given SNPs
def run_list_snps_carrying_given_snps():
    st.header("List SNPs Carrying Given SNPs")
    st.markdown(
        """
        <div class="description-box">
            <p>This function calculates, for a genomic region, the percentage of individuals carrying
            a reference allele who also carry the tested allele.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("### **Input Files**")
    st.markdown("- **VCF file** (phased/unphased & compressed/not compressed)")
    st.markdown("- **List of SNP IDs** (same format as the IDs in the VCF file)")

    st.markdown("### **Output Files**")
    st.markdown("Tables containing SNPs and their Carrier Percentage.")

    st.markdown("### **Functional part**")

    carrier_percentage_cut = st.number_input("Carrier Percentage Cut (default: 80%)", min_value=0, max_value=100, value=80)
    carrier_percentage_mode = st.selectbox("Carrier Percentage Mode", ["rough", "exact"], index=1)
    maf_threshold = st.number_input("MAF Threshold (default: 0.01)", value=0.01)

    work_directory = st.text_input("Enter the output directory for results", value="")
    if not work_directory.strip():
        st.warning("Output directory is mandatory.")
        return

    source_vcf_path = st.text_input("Enter the full file path for the VCF file", value="")
    snp_file_path = st.text_input("Enter the file path for the SNP list file", value="")

    # Ask for region to extract
    use_specific_region = st.checkbox("Extract SNPs from the vcf file for a specific analyze ? (will analyze the SNPs provided in the list)")
    region_to_extract = None
    if use_specific_region:
        region_to_extract = st.text_input("Enter the file path for the SNP list", value="")

    if st.button("Run Analysis"):
        # Create a placeholder for live output
        st.write("### **Run logs :**")
        output_placeholder = st.empty()

        # Redirect stdout to capture print statements
        output_capture = StringIO()
        sys.stdout = output_capture

        try:
            with st.spinner('Running SNPs Analysis...'):
                # Call the function with user-provided paths
                success = list_snps_carrying_given_snps(
                    carrier_percentage_cut, carrier_percentage_mode, work_directory,
                    source_vcf_path, snp_file_path, region_to_extract, maf_threshold
                    )
                if success:
                    st.success("✅ Analysis Complete! Results saved in the specified output directory.")
                else:
                    st.error("❌ Error: Analysis failed or returned incomplete results.")
        except Exception as e:
            st.error(f"An error occurred: {e}")

        # Retrieve the captured output
        output = output_capture.getvalue()
        # Display the output in Streamlit
        output_placeholder.text(output)
        # Restore stdout
        sys.stdout = sys.__stdout__

# Function to run Find Haploblocks in a Region
def run_find_Haploblocks_in_a_region():
    """Streamlit UI to find haploblocks in a genomic region."""

    # Initialize session state for toggling additional notes visibility
    if "show_note_cp" not in st.session_state:
        st.session_state["show_note_cp"] = False
    if "show_max_empty_snp_gap_note" not in st.session_state:
        st.session_state["show_max_empty_snp_gap_note"] = False
    if "show_maf_thr_note" not in st.session_state:
        st.session_state["show_maf_thr_note"] = False

    # Functions to toggle the visibility of additional explanatory notes
    def toggle_note_cp():
        st.session_state["show_note_cp"] = not st.session_state["show_note_cp"]
    def toggle_max_empty_snp_gap_note():
        st.session_state["show_max_empty_snp_gap_note"] = not st.session_state["show_max_empty_snp_gap_note"]
    def toggle_maf_thr_note():
        st.session_state["show_maf_thr_note"] = not st.session_state["show_maf_thr_note"]

    # Display the application title and input parameter section
    st.header("Find Haploblocks in a Region")

    st.markdown(
        """
        <div class="description-box">
            <p>This function builds the haploblocks within a genomic region.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("### **Input Files**")
    st.markdown("- **VCF file** (phased/unphased & compressed/not compressed)")
    st.markdown("- **List of SNP IDs** (Chromosome:Position:Position:RefAllele:AltAllele format) - **Optional**")

    st.markdown("### **Output Files**")

    st.markdown("*Tables :*")
    st.markdown("   - Basics haploblocks informations")
    st.markdown("   - Haploblocks SNP composition")
    st.markdown("   - Haploblocks sequence")
    st.markdown("   - SNPs informations")

    st.markdown("*Histograms :*")

    st.markdown("   - Haploblocks distribution based on the size in bp")
    st.markdown("   - Haploblocks distribution based on the size in SNPs")

    st.markdown("### **Functional part**")

    st.markdown("### **Input Parameters**")

    # Input fields for various haploblock search parameters
    r_square_cut = st.number_input("R-Square Cutoff (default: 0.1)", min_value=0.0, max_value=1.0, value=0.1)
    d_prime_cut = st.number_input("D Prime Cutoff (default: 0.7)", min_value=0.0, max_value=1.0, value=0.7)
    maf_percentage_cut = st.number_input("MAF Percentage Cutoff (default: 0.8)", min_value=0.0, max_value=1.0, value=0.8)
    carrier_percentage_cut = st.number_input("Carrier Percentage Cutoff (default: 80%)", min_value=0, max_value=100, value=80)
    carrier_percentage_mode = "exact"

    # Additional configuration parameters
    region_size = st.number_input("Region Size (default: 10Mb)", min_value=1, value=10000000)
    max_bp_size_for_haploblock= st.number_input("Max BP Size for Haploblock (default: 5Mb)", min_value=10, value=5000000)
    max_empty_snp_gap = st.number_input("Max Empty SNP Gap (default: 200)", min_value=10, value=200)
    maf_threshold = st.number_input("MAF Threshold (default: 0.01)", min_value=0.0, max_value=1.0, value=0.01)

    # Option to generate haploblocks for specific SNPs
    generate_specific_snps = st.selectbox("Do you want to generate Haploblocks for specific SNPs?", ["No", "Yes"], index=0)
    specific_snps_file = None
    run_mode = "1"  # Default mode set to "Classic"
    all_hap = None
    add_snp = None

    extend_mode = st.checkbox("Extend Mode (extend large haploblocks)", value=False)
    if extend_mode:
        extend_threshold = st.number_input("Extend Threshold (default: 0.9)", min_value=0.0, max_value=1.0, value=0.9)
    else:
        extend_threshold = None
    pruning_mode = st.checkbox("Pruning Mode (remove overlapping haploblocks)", value=False)
    if pruning_mode:
        overlap_threshold = st.number_input("Overlap Threshold (default: 0.95)", min_value=0.0, max_value=1.0, value=0.95)
    else:
        overlap_threshold = None

    if generate_specific_snps == "Yes":
        snp_file_path = st.text_input("Enter the full file path for the SNP list file (txt)", value="")
        if snp_file_path:
            specific_snps_file = snp_file_path
            run_mode = "3"  # List SNP mode activated
            st.info("Run mode set to 'List SNP' as you've provided a SNP list.")
            all_hap = st.selectbox("Choose List Mode", ["1 - Generate for list SNPs & ungrouped SNPs", "2 - Only generate for list SNPs"], index=0)
            all_hap = 1 if all_hap.startswith("1") else 2
            add_snp = st.selectbox("Include SNPs in Haploblocks (not as core SNPs)?", ["Yes", "No"], index=0) == "Yes"
    else:
        st.info("You chose not to generate Haploblocks for specific SNPs.")
        run_mode = st.selectbox("Choose Run Mode", ["Standard", "Exhaustive"], index=0)
        run_mode = "1" if run_mode == "Standard" else "2"

    ld_parameter = "1" 

    # User inputs for output directory and VCF file path
    work_directory = st.text_input("Enter the output directory for results", value="")
    source_vcf_filepath = st.text_input("Enter the full file path for the VCF file", value="")

    # Ensure mandatory inputs are provided
    if not work_directory.strip() or not source_vcf_filepath.strip():
        st.warning("⚠️ Both output directory and VCF file path are mandatory.")
        return

    # Execute haploblock search when the user clicks the button
    if st.button("Run Haploblock Search"):
        # Create a placeholder for live output
        st.write("### **Run logs :**")
        output_placeholder = st.empty()

        # Redirect stdout to capture print statements
        output_capture = StringIO()
        sys.stdout = output_capture

        try:
            result = find_haploblocks_in_a_region(
                        r_square_cut, d_prime_cut, maf_percentage_cut, carrier_percentage_cut,
                        carrier_percentage_mode, region_size, max_bp_size_for_haploblock,
                        max_empty_snp_gap, work_directory,
                        source_vcf_filepath, specific_snps_file, run_mode, ld_parameter, all_hap, add_snp, maf_threshold,
                        extend_mode=extend_mode, extend_threshold=extend_threshold,
                        pruning_mode=pruning_mode, overlap_threshold=overlap_threshold
                    )

            st.session_state['analysis_result'] = result
        except Exception as e:
            st.session_state['analysis_error'] = e

        # Retrieve the captured output
        output = output_capture.getvalue()
        # Display the output in Streamlit
        output_placeholder.text(output)
        # Restore stdout
        sys.stdout = sys.__stdout__

        # Display results or error messages
        if 'analysis_error' in st.session_state:
            st.error(f"❌ An error occurred: {st.session_state['analysis_error']}")
            st.text_area("Debug Traceback", traceback.format_exc(), height=300)
        elif 'analysis_result' in st.session_state:
            result = st.session_state['analysis_result']
            if result:
                st.success("✅ Haploblock search completed successfully!")
            else:
                st.error("❌ Error: Haploblock search failed or returned incomplete results.")

        # Cleanup session state to avoid storing unnecessary data
        if 'analysis_result' in st.session_state:
            del st.session_state['analysis_result']
        if 'analysis_error' in st.session_state:
            del st.session_state['analysis_error']

# Function to plot Haploblocks
def run_print_Haploblocks():
    st.header("Plot Haploblocks")
    
    st.markdown(
        """
        <div class="description-box">
            <p>Generates images to visualize haploblocks on a chromosome. Each haploblock is represented as a rectangle, with its core SNP marked by a star.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
    
    st.markdown("""### **Input file**""")
    st.markdown("""Haploblock data file generated by the `find_haploblocks_in_a_region` function.""")

    # Input for the output directory
    work_directory = st.text_input("Enter the output directory for results", value="")

    # Input for the Haploblock file path
    haploblock_file_path = st.text_input("Enter the full file path for the Haploblock composition file", value="")
    snp_file_path = st.text_input("Enter the full file path for the SNP-MAF file", value="")

    # Inputs for filtering criteria with validation
    use_min_size_snp = st.checkbox("Set a minimum SNP size for the haploblocks to plot.")
    if use_min_size_snp:
        min_size_snp = st.number_input("Minimum size in SNPs :", min_value=0, step=1)
    else:
        min_size_snp = None

    use_max_size_snp = st.checkbox("Set a maximum SNP size for the haploblocks to plot.")
    if use_max_size_snp:
        max_size_snp = st.number_input("Maximum size in SNPs :", min_value=0, step=1)
    else:
        max_size_snp = None

    use_min_size_bp = st.checkbox("Set a minimum base pairs size for the haploblocks to plot.")
    if use_min_size_bp:
        min_size_bp = st.number_input("Minimum size in base pairs :", min_value=0, step=1)
    else:
        min_size_bp = None

    use_max_size_bp = st.checkbox("Set a maximum base pairs size for the haploblocks to plot.")
    if use_max_size_bp:
        max_size_bp = st.number_input("Maximum size in base pairs :", min_value=0, step=1)
    else:
        max_size_bp = None

    use_start_bp = st.checkbox("Set a starting position of the range in bp to plot the haploblocks.")
    if use_start_bp:
        bp_range_start = st.number_input("Starting position in base pairs :", min_value=0, step=1)
    else:
        bp_range_start = None

    use_end_bp = st.checkbox("Set an ending position of the range in bp to plot the haploblocks.")
    if use_end_bp:
        bp_range_end = st.number_input("Ending position in base pairs :", min_value=0, step=1)
    else:
        bp_range_end = None

    use_maf = st.checkbox("Set a MAF percentage to highlight in red the SNPs with a MAF that is +/- ")
    if use_maf:
        maf_threshold = st.number_input("MAF percentage (between 0 and 1)", min_value=0.0, max_value=1.0, step=0.01)
    else:
        maf_threshold = None

    # Ensure the required inputs are provided
    if not work_directory.strip() or not haploblock_file_path.strip():
        st.warning("Both output directory and Haploblock file path are mandatory.")
        return

    # Button to plot Haploblocks
    if st.button("Plot Haploblocks"):
        # Create a placeholder for live output
        st.write("### **Run logs :**")
        output_placeholder = st.empty()

        # Redirect stdout to capture print statements
        output_capture = StringIO()
        sys.stdout = output_capture

        try:
            with st.spinner('Generating the image...'):
                # Dictionary of possible options
                bloc_options = {
                    "min_size_snp": min_size_snp if min_size_snp else None,
                    "max_size_snp": max_size_snp if max_size_snp else None,
                    "min_size_bp": min_size_bp if min_size_bp else None,
                    "max_size_bp": max_size_bp if max_size_bp else None,
                    "bp_range_start": bp_range_start if bp_range_start else None,
                    "bp_range_end": bp_range_end if bp_range_end else None,
                    "maf_threshold": maf_threshold if maf_threshold else None
                }

                # Base file title
                file_title = f"{work_directory}/Haploblocks_plot"

                # Build the graph name
                plotting_parameters = {}
                for parameter, value in bloc_options.items():
                    if value:  # Include only non-empty or non-None values
                        plotting_parameters[parameter] = value
                        file_title += f".{parameter}_{value}"
                file_title += ".png"

                # Call the main function
                success = print_haploblocks(
                    work_directory=work_directory,
                    haploblock_file=haploblock_file_path,
                    snp_file=snp_file_path,
                    **plotting_parameters
                )
                if success == True:
                    st.success(f"✅ Haploblock plot created successfully!")
                else:
                    st.error("❌ Error: Haploblock plot creation failed.")
        except Exception as e:
            st.error(f"An error occurred: {e}")

        # Retrieve the captured output
        output = output_capture.getvalue()
        # Display the output in Streamlit
        output_placeholder.text(output)
        # Restore stdout
        sys.stdout = sys.__stdout__

# Main app routing based on sidebar choice
if option == "Home Page":
    run_homepage()
elif option == "Find Haploblocks for SNPs":
    run_find_Haploblocks_for_given_snps()
elif option == "List SNPs Carrying Given SNPs":
    run_list_snps_carrying_given_snps()
elif option == "Find Haploblocks in a Region":
    run_find_Haploblocks_in_a_region()
elif option == "Plot Haploblocks":
    run_print_Haploblocks()


# Footer content
st.markdown("""
    ---
    **HaploExplore** | Developed by Samuel Hiet, GBCM Team - CNAM, Paris""")
    
st.markdown("""           
    Version: 1.0 | Last Update: January 2025
    
    For questions or feedback, contact: samuel.hiet@lecnam.net
""")

st.markdown("""
    ---
    **Follow us for updates:**  
    [GitHub](link) | [LinkedIn] (link)
""")

st.markdown("""
    **Disclaimer:** This tool is intended for research purposes only.  
    Results are not guaranteed to be accurate for clinical decision-making.
""")

st.markdown("""
    ---
    <div style="text-align: center; font-size: 12px;">
        © 2025 HaploExplore. All rights reserved.
    </div>
""", unsafe_allow_html=True)

