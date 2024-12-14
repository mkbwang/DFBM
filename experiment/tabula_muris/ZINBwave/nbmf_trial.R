
rm(list=ls())
count_mat <- read.table("experiment/tabula_muris/ZINBwave/simulated_counts.tsv",
                        header=TRUE, row.names=1, sep="\t")

example_slice <- count_mat > 11
example_mask <- count_mat > 10

begin <- proc.time()
result <- nbmf(Y = example_slice, Z = example_mask, k=3,
               tol=1e-5)
end <- proc.time()
