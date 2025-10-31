HAPLOEXPLORE

Motivation: Haploblocks in the genome are constitutive of evolution patterns and they play a pivotal role in shaping the genomic variability and susceptibility/resistance to diseases. Several software have been developed for haploblock detection, but they do not distinguish between the impacts of major and minor SNP alleles. In this study, we present a powerful haploblock detection software, specifically designed for identifying haploblocks associated with SNP minor allele haploblocks (MiA-haploblocks). The focus on minor alleles is essential since they represent the most recent evolutionary signatures and are closely linked to the various selection pressures experienced in recent history. 
Results: HaploExplore operates on vcf files containing phased data, exhibiting rapid processing times (a few minutes for analyzing the whole chromosome 22) and generating user-friendly outputs. Its results are convergent for populations starting from 100 individuals, and, as expected, shorter MiA-haploblocks are observed in populations of African descent compared to those of European descent. A comparative analysis of HaploExplore against other haploblock detection software revealed its superiority in terms of either simplicity, or flexibility, or speed, with the unique capability to target minor alleles. HaploExplore will be very useful for evolutionary genomics and for GWAS analysis in human diseases, given that the effects of genetic associations may accumulate within a specific haploblocks.
__________________________________________________________________________________________________________________________________________________________________________________________________

A Streamlit-based web application has been developed for an interactive and user-friendly experience. Users can easily upload input files, adjust parameters, and visualize results without running command-line scripts. The app is accessible via Docker or can be launched locally.

Running with Docker: 

Set the Volumes path :
	- set the accessible Folders for your Docker in the file .env

How to build the docker image (app) ?
$ docker compose -f 'docker-compose.yml' up -d --build
	or
Run the script "docker-compose.yml"

How to run the docker image (app) ?
$ docker run -p 8501:8501 haploexplore
	or
Run the docker image via Docker desktop

To run easily the application we recommend building it on Visual studio and using Docker Desktop.


Running without Docker (requires bcftools): 

To run HaploExplore only with Streamlit, ensure you have Python 3.10.12 installed along with bcftools, which is required for processing VCF files. 
To use the graphical interface of HaploExplore, install the required dependencies and launch the Streamlit app:

	1. Install dependencies:

		pip install -r requirements.txt

	2. Run the application:

		streamlit run app/application.py

	3. Open the displayed localhost link in a web browser.


Warning : If the App can't access to the provided folders : Go to 1. DockerDesktop
								  2. Settings
								  3. File Sharing
								  4. Add the directory path
__________________________________________________________________________________________________________________________________________________________________________________________________

A version of HaploExplore exists without the application (program-only) and can be run directly from the terminal with the same functionalities as the application.  
To run it, an "execution.py" file is provided to show how to use the different functions. Also to set the different parameters, a config.json file is provided.

This version requires Python 3.10.12 and bcftools, which is needed for processing VCF files.
__________________________________________________________________________________________________________________________________________________________________________________________________

Current functionalities:

	- Haploblock creation with LD (Find haploblocks in a region)
		- Building modes :
			- Standard mode : if a haploblock cannot be added to any haplobloc then a new one is generated.
			- Exhaustive mode : build a haploblock for each SNPs.
			- ListSNP mode : build a haploblock for SNPs of the provided list only. (possibility to integrate the SNPs of the list in haploblocks or not (useful for the study of genes for example) and also to generate haploblocks for other SNPs that are not in the list. (in addition of the SNPs of the list) such as the standard mode)
		- Removal of SNPs with a MAF <  0.01 (this parameter can by changed)

	- Extract the Haploblocs where one (or more) SNP is contained from a result file. (Find Haploblocks for Given SNPs)
	
	- Generate a table with the SNPs that are carried by the SNPs of the list. (List SNPs Carrying Given SNPs)

	- Generate a graphic to visualize haploblocks on a chromosome. Each haploblock is represented as a rectangle, with its core SNP marked by a star. (Plot haploblocks)

	- All parameters can be changed by the user.
__________________________________________________________________________________________________________________________________________________________________________________________________

Input files : 

	- Find haploblocks in a region function : 
		- A VCF file containing the region to analyze (no limit of individuals).
		- A list of SNPs ID (Chromosome:Position:Position:RefAllele:AltAllele format) in a text file - for ListSNP mode only.

	- Find Haploblocks for Given SNPs function :
		- A list of SNP IDs. (Chromosome:Position:Position:RefAllele:AltAllele) in a text file.
		- A composition file generated by "Find haploblocks in a region" function.

	- List SNPs Carrying Given SNPs function :
    		- VCF file. (phased/unphased & compressed/not compressed)
    		- List of SNP IDs. (Chromosome:Position:Position:RefAllele:AltAllele format)

	- Plot Haploblocks function :
		- A composition file generated by "Find haploblocks in a region" function.
		- A SNP information file generated by "Find haploblocks in a region" function.
__________________________________________________________________________________________________________________________________________________________________________________________________

Default Settings and Parameters:
HaploExplore includes several default parameters that can be modified through the app:

. LD Thresholds:
	r²: 0.1 (default)
	D': 0.7 (default)
. Carrier Percentage (CP): 80% (default)
. MAF percentage cut: 0.8 (default)
. Maximum SNP Gap Within a Block: 200 SNPs (default)
. Maximum Haploblock Size: 5,000,000 base pairs (default)
. Region Size for Splitting Datasets: 10,000,000 base pairs (default)
. Minimum MAF Threshold to Include a SNP: 1%(default) - when loading the VCF file
. Extend mode : 0.9 (default) - extend haploblocks (without region and size limits) when its size reach x% (per default 90%) of Maximum Haploblock Size
. Pruning mode : 0.95 (default) - delete haploblocks that have at least x% (per default 95%) of its SNPs that are contained in another haploblocks (of the same size or bigger)


Note: The Carrier Percentage (CP) is the proportion of individuals who carry the minor allele of a tag SNP and also carry the minor alleles of other SNPs within the same haploblock. This ensures that only SNPs with strong biological relevance and statistical correlation are grouped together.

__________________________________________________________________________________________________________________________________________________________________________________________________

Output files :

	- Find haploblocks in a region function : 
		- Tables (.txt):
    			- Basics haploblocks informations
			- Haploblocks SNP composition
			- Haploblocks sequence
			- SNPs informations (MAF & BP)
		- Histograms (.pdf):
			- Haploblocks distribution based on the size in bp
			- Haploblocks distribution based on the size in SNPs

	- Find Haploblocks for Given SNPs function :
		- Text file containing for each SNP the different haploblocks it can be found in.

	- List SNPs Carrying Given SNPs function :
		- Tables containing SNPs and their Carrier Percentage. (.txt)

	- Plot Haploblocks function :
		- Graphics (.png)


The software is built with Python 3.10.12 and utilizes Streamlit for an interactive graphical interface.
