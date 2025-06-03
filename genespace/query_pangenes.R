##load the prerequisite modules: paste in terminal
#eval "$(conda shell.bash hook)"
#conda activate genespace

devtools::install_github("jtlovell/GENESPACE")
library(GENESPACE)

#change the working directory path to match your own
wd = "/scratch/ac05869/gs_pcu_cro/wd"
path2mcscanx = '/home/ac05869/miniconda/envs/genespace/bin'

load(file = "/scratch/ac05869/gs_pcu_cro/wd/results/gsParams.rda")

roi <- data.frame(
  genome = c("cro_v3", "pcu_v0_h1", "pcu_v0_h2"),
  chr = c("Chr3", "scf_4_1", "scf_4_2"),
  start = c(71588447, 87471718, 98404806),
  end = c(71691498, 87510312, 98461144))


query_pangenes(
  gsParam,
  bed = NULL,
  refGenome = NULL,
  transform = TRUE,
  showArrayMem = TRUE,
  showNSOrtho = TRUE,
  maxMem2Show = Inf,
  showUnPlacedPgs = FALSE
)

test <- query_pangenes(
  gsParam = gpar, bed = roi)