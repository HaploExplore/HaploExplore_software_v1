class SNP :
    """
    A Single Nucleotide Polymorphism defined by : its id, chromosome, position in bp, ref and alt alleles.
    """
    def __init__(self, rs_id : str, chromosome : int, position : int, ref_allele : str, alt_allele : str):
        self.__rs_id = rs_id
        self.__chromosome = chromosome
        self.__position = position
        self.__ref_allele = ref_allele
        self.__alt_allele = alt_allele

    def get_rs_id(self) -> str :
        return self.__rs_id

    def get_chromosome(self) -> int :
        return self.__chromosome

    def get_position(self) -> int :
        return self.__position

    def get_ref_allele(self) -> str :
        return self.__ref_allele

    def get_alt_allele(self) -> str :
        return self.__alt_allele

    def __str__(self) -> str :
        return f"SNP - rsID : {self.__rs_id}, chromosome : {self.__chromosome}, position : {self.__position} bp, ref : {self.__ref_allele}, alt : {self.__alt_allele}."

class SNPInformations :
    """
    A SNP within a context defined by : the SNP, its Minor Allele Frequency, data phasing, the list of individuals genotyped and of those that carry the minor allele on each chromosome.
    """
    def __init__(self, snp : SNP, maf : float, is_phased : bool, individuals : list[int], alt_carriers : list[int], alt_carriers_1 : list[int] | None, alt_carriers_2 : list[int] | None, ref_carriers : list[int], ref_carriers_1 : list[int] | None, ref_carriers_2 : list[int] | None):
        self.__snp = snp
        self.__maf = maf
        self.__phased = is_phased
        self.__individuals = individuals
        self.__alt_allele_carriers = alt_carriers
        self.__alt_allele_carriers_chr_1 = alt_carriers_1
        self.__alt_allele_carriers_chr_2 = alt_carriers_2
        self.__ref_allele_carriers = ref_carriers
        self.__ref_allele_carriers_chr_1 = ref_carriers_1
        self.__ref_allele_carriers_chr_2 = ref_carriers_2


    def get_snp(self) -> SNP :
        return self.__snp

    def get_maf(self) -> float :
        return self.__maf

    def is_phased(self) -> bool :
        return self.__phased

    def get_individuals_genotyped(self) -> list[int] :
        return self.__individuals

    def get_individuals_count(self) -> int :
        return len(self.__individuals)

    def get_alt_allele_carriers(self) -> list[int] :
        return self.__alt_allele_carriers

    def get_alt_allele_carriers_chr_1(self) -> list[int] :
        return self.__alt_allele_carriers_chr_1

    def get_alt_allele_carriers_chr_2(self) -> list[int] :
        return self.__alt_allele_carriers_chr_2

    def get_ref_allele_carriers(self) -> list[int] :
        return self.__ref_allele_carriers

    def get_ref_allele_carriers_chr_1(self) -> list[int] :
        return self.__ref_allele_carriers_chr_1

    def get_ref_allele_carriers_chr_2(self) -> list[int] :
        return self.__ref_allele_carriers_chr_2

class Haploblock :
    """
    A haploblock defined by its core SNP and its informations, as well as for all of the SNPs it contains.
    """
    def __init__(self, core_snp : SNPInformations, all_snps : list[SNPInformations], region:int):
        self.__core_snp = core_snp
        self.__all_snps = all_snps
        self.__region = region

    def get_core_snp(self) -> SNPInformations :
        return self.__core_snp

    def get_region(self) -> int :
        return self.__region
    
    def get_core_snp_index(self) -> int :
        for i, snp in enumerate(self.__all_snps) :
            if snp == self.get_core_snp():
                return i
        return -1

    def get_other_snps(self) -> list[SNPInformations] :
        all_snps = self.__all_snps[:]
        all_snps.pop(self.get_core_snp_index())
        return all_snps

    def get_all_snps(self) -> list[SNPInformations]:
        return self.__all_snps

    def get_start_bp(self) -> int :
        return min(snp.get_snp().get_position() for snp in self.__all_snps)

    def get_start_snp(self) -> SNP :
        def get_bp(snp : SNPInformations) :
            return snp.get_snp().get_position()
        all_snps = self.__all_snps[:]
        all_snps.sort(
            key = get_bp
        )
        return all_snps[0]

    def get_end_bp(self) -> int :
        return max(snp.get_snp().get_position() for snp in self.__all_snps)
    
    def get_end_snp(self) -> SNP :
        def get_bp(snp : SNPInformations) :
            return snp.get_snp().get_position()
        all_snps = self.__all_snps[:]
        all_snps.sort(
            key = get_bp,
            reverse = True)
        return all_snps[0]

    def get_snp_size(self) -> int :
        return len(self.__all_snps)

    def get_bp_size(self) -> int :
        return self.get_end_bp() - self.get_start_bp()

    def contains_snp(self, snp : str) -> bool :
        for snp_in_bloc in self.get_all_snps() :
            if snp_in_bloc.get_snp().get_rs_id() == snp or str(snp_in_bloc.get_snp().get_position()) == snp :
                return True
        return False
    
    def filter_snps_by_distance(self, max_snp_gap: int, dataset_snp_id: list[str]):
        """
        Filter SNPs within a haploblock by ensuring the gap between consecutive SNPs 
        (based on their index in the reference dataset) does not exceed `max_snp_gap`.
        """

        # Validate input: all SNP IDs must be strings
        if not all(isinstance(snp_id, str) for snp_id in dataset_snp_id):
            raise ValueError("All elements in dataset_snp_id must be of type str.")

        # Precompute a mapping from SNP ID -> position in the reference dataset
        snp_index_map = {snp_id: idx for idx, snp_id in enumerate(dataset_snp_id)}

        # Ensure SNPs are sorted by genomic position before filtering
        self.__all_snps.sort(key=lambda snp_info: snp_info.get_snp().get_position())

        # Identify the core SNP index (anchor SNP around which filtering is applied)
        core_index = self.get_core_snp_index()
        if core_index == -1:
            raise ValueError("The core SNP was not found in the list of SNPs.")

        core_snp_info = self.get_core_snp()

        def filter_direction(snps, reverse=False, prev_index=None):
            """
            Generator that yields SNPs in one direction (before or after the core SNP),
            stopping when the allowed SNP gap is exceeded.

            Args:
                snps (list): Subset of SNPs to traverse.
                reverse (bool): If True, traverse backwards (before core SNP).
                prev_index (int | None): Reference index of the previously accepted SNP.

            Yields:
                SNPs that satisfy the distance constraint relative to the previous SNP.
            """
            iterable = (reversed(snps) if reverse else snps)

            for snp_info in iterable:
                rs_id = snp_info.get_snp().get_rs_id()
                current_index = snp_index_map.get(rs_id, -1)

                # Compute the SNP gap relative to the previous SNP
                if prev_index is not None and prev_index != -1 and current_index != -1:
                    snp_gap = abs(current_index - prev_index) - 1
                    if snp_gap > max_snp_gap:
                        break  # Stop if the allowed gap is exceeded

                yield snp_info
                prev_index = current_index

        # Collect SNPs before the core SNP (reversed back into natural order)
        snps_before_gen = filter_direction(self.__all_snps[:core_index], reverse=True)
        snps_before = list(snps_before_gen)[::-1]

        # Collect SNPs after the core SNP
        snps_after_gen = filter_direction(
            self.__all_snps[core_index + 1:],
            reverse=False,
            prev_index=snp_index_map.get(core_snp_info.get_snp().get_rs_id(), -1)
        )
        snps_after = list(snps_after_gen)

        # Rebuild the SNP list: [before core] + [core SNP] + [after core]
        self.__all_snps = snps_before + [core_snp_info] + snps_after

