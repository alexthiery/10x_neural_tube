#SBATCH -t 6:00:00
#SBATCH --mem=0
#SBATCH -c 10
#SBATCH --job-name=10x_neural
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

ml purge
ml Anaconda2
ml R/3.6.0-foss-2019a
source activate 10x_new
cd ~/working/alexthiery/analysis/10x_scRNAseq_2019/repo
Rscript ./scripts/1_seurat_full.R CAMP