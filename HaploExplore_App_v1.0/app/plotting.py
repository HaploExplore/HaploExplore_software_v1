from matplotlib import pyplot as plt
import numpy as np
import logging
from math import ceil

colormapp = plt.get_cmap('hsv')

def plot_haploblocks(save_path: str, snp_maf_dict: dict[str, tuple[float, int]], 
                      haploblock_composition: dict[str, tuple[str, tuple[str, ...]]], 
                      maf_threshold: float = 0.1) -> bool:
    """
    Plot haploblocks on a graphic.

    Parameters:
    - save_path: str, path to save the plot
    - snp_maf_dict: dict, mapping of SNP IDs to (MAF, position) tuples
    - haploblock_composition: dict, mapping of haploblock IDs to (Core SNP, list of SNP IDs)
    - maf_threshold: float, threshold for MAF comparison
    """
    try:
        if not save_path:
            raise ValueError("save_path must be provided")
        if not snp_maf_dict:
            raise ValueError("snp_maf_dict is empty or not provided")
        if not haploblock_composition:
            raise ValueError("haploblock_composition is empty or not provided")

        # Validate SNP data before sorting
        valid_haploblocks = {}
        for haploblock_id, (core_snp, snps) in haploblock_composition.items():
            if not isinstance(haploblock_id, str):
                raise TypeError(f"Haploblock ID must be a string, got {type(haploblock_id)}")

            if not snps:
                logging.warning(f"Warning: Empty SNP list for haploblock {haploblock_id}")
                continue

            # Ensure SNPs are valid and convert list to tuple
            valid_snps = tuple(snp for snp in snps if snp in snp_maf_dict and all(x is not None for x in snp_maf_dict[snp]))
            if not valid_snps:
                logging.warning(f"Warning: No valid SNPs found for haploblock {haploblock_id}")
                continue

            valid_haploblocks[haploblock_id] = valid_snps

        if not valid_haploblocks:
            raise ValueError("No valid haploblocks found after filtering")

        # Sorting haploblocks by the first SNP's position
        haploblocks_sorted = sorted(
            valid_haploblocks.items(),
            key=lambda x: snp_maf_dict.get(x[1][0], (0, float('inf')))[1] if x[1] else float('inf')
        )

        # Plot settings
        line_thickness = 3
        bloc_intervals = 62 * line_thickness
        starting_y = 150
        num_haploblocks = len(valid_haploblocks)
        fig_height = max(5, min(15, num_haploblocks * 0.5))
        plt.figure(figsize=(15, fig_height), dpi=300)

        # Find min/max SNP positions for x-axis limits
        all_positions = [snp_maf_dict[snp][1] for snps in valid_haploblocks.values() for snp in snps if snp in snp_maf_dict and snp_maf_dict[snp][1] is not None]
        if not all_positions:
            raise ValueError("No valid SNP positions found")

        min_pos, max_pos = min(all_positions), max(all_positions)
        x_padding = (max_pos - min_pos) * 0.05
        min_pos -= x_padding
        max_pos += x_padding

        plt.hlines(0, xmin=min_pos, xmax=max_pos, color='black', alpha=0.75, linewidth=1.5, zorder=2, label="_nolegend_")
        
        # Tick marks on x-axis
        num_ticks = 10
        tick_positions = [min_pos + (max_pos - min_pos) * i / (num_ticks - 1) for i in range(num_ticks)]
        for tick in tick_positions:
            plt.vlines(tick, ymin=-50, ymax=20, color='black', alpha=1, linewidth=0.8, zorder=1, label="_nolegend_")
            plt.text(tick, -60, f"{int(tick):,}", fontsize=8, ha='center', va='top', zorder=2)

        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        legend_labels = []

        for haploblock_idx, (haploblock_id, snps) in enumerate(haploblocks_sorted, start=1):
            tag_snp_id = snps[0]  # Core SNP is the first SNP in the list
            tag_snp_data = snp_maf_dict.get(tag_snp_id)
            if not tag_snp_data or None in tag_snp_data:
                continue
            tag_snp_maf, tag_snp_bp = tag_snp_data

            snp_positions = [snp_maf_dict[snp][1] for snp in snps if snp in snp_maf_dict and snp_maf_dict[snp][1] is not None]
            if not snp_positions:
                continue

            bloc_start_bp = min(snp_positions)
            bloc_end_bp = max(snp_positions)
            bloc_center_y = starting_y + (haploblock_idx - 1) * (line_thickness + bloc_intervals)
            plt.hlines(bloc_center_y, bloc_start_bp, bloc_end_bp, color="blue", alpha=0.75, linewidth=line_thickness, zorder=3)
            plt.scatter(tag_snp_bp, bloc_center_y, color='black', marker='*', alpha=0.75, linewidth=0.1, zorder=4, s=30, label="_nolegend_")
            plt.text(tag_snp_bp, bloc_center_y + 15, f"{haploblock_idx}", fontsize=6, color='black', alpha=0.75, zorder=4)
            legend_labels.append(f"Haploblock {haploblock_idx}: Tag SNP = {tag_snp_id}, start BP = {bloc_start_bp}, end BP = {bloc_end_bp}")

            for snp_id in snps[1:]:
                snp_data = snp_maf_dict.get(snp_id)
                if not snp_data or None in snp_data:
                    continue
                snp_maf, snp_pos = snp_data
                if maf_threshold is not None:
                    if abs(snp_maf - tag_snp_maf) <= maf_threshold * tag_snp_maf:
                        plt.hlines(bloc_center_y, snp_pos - 500, snp_pos + 500, color='red', linewidth=4, zorder=3, label="_nolegend_")

            plt.vlines(tag_snp_bp, ymin=0, ymax=bloc_center_y, colors='blue', alpha=0.5, linestyles='dashed', linewidth=0.5, zorder=2, label="_nolegend_")

        if legend_labels:
            plt.legend(labels=legend_labels, bbox_to_anchor=(0.5, -0.2), loc="upper center", ncol=3, fontsize=7, frameon=False)

        plt.yticks([])
        plt.ylabel("Haploblock Line")
        plt.xlabel("Position (bp)")
        plt.title("Haploblock Composition")
        plt.tight_layout(rect=[0, 0.25, 1, 1])
        plt.savefig(save_path, bbox_inches="tight")
        plt.close()
        return True

    except Exception as e:
        logging.error(f"Error: {str(e)}")
        plt.close()
        return False

def generate_haploblock_histograms(file_name : str, output1 : str, output2 : str):
    def generate_snp_histogram(file_name : str, output1 : str):
        try:
            # Initialize lists to store SNP counts
            snp_counts = []

            # Read the file line by line, skipping the header
            with open(file_name, 'r') as f:
                next(f)  # Skip the header row
                for line in f:
                    values = line.strip().split("\t")
                    snp_count = int(values[-1])  # Get the last column (size_in_snps)
                    snp_counts.append(snp_count)

            # Define the SNP count ranges for the histogram
            snp_count_ranges = [1, 2, 3, 4, 5, 10, 50, 100, 500, 1000, float("inf")]
            labels = ["1","2","3","4","5-10", "10-50", "50-100", "100-500", "500-1000", "> 1000"]

            # Count the number of haploblocks in each SNP count range
            haplo_counts, _ = np.histogram(snp_counts, bins=[x for x in snp_count_ranges])

            # Generate the histogram for haploblock SNP counts
            plt.figure(figsize=(12, 10))
            plt.bar(labels, haplo_counts)
            plt.xlabel("Haploblock SNP Count")
            plt.ylabel("Number of Haploblocks")
            plt.title("Haploblock SNP Count Distribution")
            plt.xticks(range(len(labels)), labels, rotation=60, ha='right')

            # Add the number of haploblocks on the bars
            for i, count in enumerate(haplo_counts):
                plt.text(i, count+3, str(count), ha='center', va='bottom')

            plt.savefig(output1)
            plt.clf()
            logging.info(f"Saved SNP count histogram to {output1}")
            return output1
        except Exception as e:
            logging.error(f"Error generating SNP histogram: {e}")
            return False

    def generate_size_histogram(file_name : str, output2 : str):
        try:
            # Initialize lists to store start and end positions
            start_bp = []
            end_bp = []

            # Read the file line by line, skipping the header
            with open(file_name, 'r') as f:
                next(f)  # Skip the header row
                for line in f:
                    # Split the line into columns
                    cols = line.strip().split("\t")
                    # Append the start and end positions to the lists
                    start_bp.append(int(cols[4]))
                    end_bp.append(int(cols[6]))

            # Calculate the size of haploblocks in kb
            sizes_in_kb = [(end - start) / 1000 for start, end in zip(start_bp, end_bp)]

            # Define the size ranges for the histogram
            size_ranges = [0, 100, 200, 500, 1e3, 2e3, 3e3, 4e3, 5e3, float("inf")]
            labels = ["1-100 kb", "100-200 kb", "200-500 kb", "500 kb-1 MB", "1-2 MB", "2-3 MB", "3-4 MB", "4-5MB", "> 5 MB"]

            # Count the number of haploblocks in each size range
            haplo_counts, _ = np.histogram(sizes_in_kb, bins=[x for x in size_ranges])

            # Generate the histogram for haploblock sizes
            plt.figure(figsize=(12, 10))
            plt.bar(labels, haplo_counts)
            plt.xlabel("Haploblock Size (kb)")
            plt.ylabel("Number of haploblocks")
            plt.title("Haploblock Size Distribution")
            plt.xticks(range(len(labels)), labels, rotation=60, ha='right')

            # Add the number of haploblocks on the bars
            for i, count in enumerate(haplo_counts):
                plt.text(i, count+3, str(count), ha='center', va='bottom')

            plt.savefig(output2)
            plt.clf()
            logging.info(f"Saved haploblock size histogram to {output2}")
            return output2
        except Exception as e:
            logging.error(f"Error generating size histogram: {e}")
            return False

    snp_success = generate_snp_histogram(file_name, output1)
    size_success = generate_size_histogram(file_name, output2)

    return snp_success, size_success


