import main_script as haploblocks
import json
import time

with open("config.json", 'r') as file:
    config = json.load(file)

def log_time(message):
    timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
    return f"[{timestamp}] {message}"

r_square_cut = config["r_square_cut"][0]
d_prime_cut = config["d_prime_cut"][0]
maf_percentage_cut = config["maf_percentage_cut"][0]
carrier_percentage_cut = config["carriere_percentage_cut"][0]
MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC = config["MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC"][0]
MAX_BP_SIZE_FOR_HAPLOBLOC = config["MAX_BP_SIZE_FOR_HAPLOBLOC"][0]
REGION_SIZE=config["REGION_SIZE"][0]
maf_threshold = config["maf_threshold"][0]

print(f"Using r_square_cut: {r_square_cut}, d_prime_cut: {d_prime_cut}, maf_percentage_cut: {maf_percentage_cut}, carriere_percentage_cut: {carrier_percentage_cut},  Max_empty_gap: {MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC}, Max_BP_size: {MAX_BP_SIZE_FOR_HAPLOBLOC}, REGION_SIZE: {REGION_SIZE}")
start = log_time("Starting...")
print(start)

haploblocks.find_haploblocks_in_a_region(r_square_cut, d_prime_cut, maf_percentage_cut, carrier_percentage_cut, REGION_SIZE, MAX_BP_SIZE_FOR_HAPLOBLOC, MAX_EMPTY_SNP_GAP_BETWEEN_TWO_SNPS_IN_A_BLOC, maf_threshold)
#haploblocks.find_haploblocks_for_given_snps()
#haploblocks.print_haploblocks()
#haploblocks.list_snps_carrying_given_snps(carrier_percentage_cut, maf_threshold)

end = log_time("Analyze done !")
print("start : ", start)
print("end : ",end)
