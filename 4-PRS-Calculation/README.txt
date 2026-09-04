To use PGS-Calculator

1. Clone PGS Calculation Pipeline to this directory 
$ git clone https://github.com/mattlee821/PGS.git

2. Navigate inside the repository 
$ cd PGS

3. Run setup
$ bash PGS/src/setup.sh

4. Create params.yml based on example i.e.:
i) path to genetic data
ii) samplest name
iii) target_build
iv) directory to cache Apptainer container images
v) HPC module name for Apptainer
vi) min overlap for scores

5. Alterations:
i) change partition and account number in PGS/src/run.sh
ii) change account number in PGS/pipeline/pgsc_calc.config
iii) change "bfile" to "pfile" in PGS/pipeline/PGS.sh (line 150)

6. Run the calculator:
$ sbatch --job-name=pgs-[PGS_ID] PGS/src/run.sh --trait "[PGS_ID]" --dir_out PATH/TO/Data/Intermediate/4-PRS

X. For further guidance see the repository's readme