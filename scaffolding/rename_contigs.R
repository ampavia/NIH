setwd("/scratch/ac05869/gese_final_yahs/gese_v2.asm/100kb")
library(readr)
library(dplyr)
library(stringr)

gese_v2 <- read_tsv("./gese_v2.tsv", col_names = c("name", "seq"))
renamed <- gese_v2  %>%
  mutate(order = c(1:46)) %>%
  mutate(order2 = if_else(order > 8, cumsum(order > 8), NA_integer_)) %>%
  mutate(name = case_when(
    name == "scaffold_1" ~ ">Chr1",
    name == "scaffold_2" ~ ">Chr2",
    name == "scaffold_3" ~ ">Chr3",
    name == "scaffold_4" ~ ">Chr4",
    name == "scaffold_5" ~ ">Chr5",
    name == "scaffold_6" ~ ">Chr6",
    name == "scaffold_7" ~ ">Chr7",
    name == "scaffold_8" ~ ">Chr8",
    T ~ "NA")) %>%
  mutate(name = case_when(name == "NA" ~ paste(">scaffold", order2),
                                                     T ~ name)) %>%
  select(name, seq)

renamed$name <- gsub(pattern = " ", replacement = "_", renamed$name)
tail(renamed)
write_tsv(renamed, file = "./gese_v2_renamed.tsv", col_names = F)