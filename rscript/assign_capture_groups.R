library(tidyverse)

meta <- read_tsv("FA01p/Pigs_allTo20042016_shared.tsv", col_types = cols(.default = "c")) %>%
  select(sample="Extract No.", age="Age", period="Period", neolithic="Neolithic (Y/N/W)", status="Wild/Dom Status", internal="Internal", external="External") %>%
  unique()

groups <- read_tsv("FA01p/FA01p_NewGroupIDs_200217_EIP.txt", col_types = cols()) %>%
  select(sample=`Sample ID`, type=`Lib or Ext`, old_group=`New Group ID`)

ssDNA <- read_tsv("FA01p/FA01p_endo_ssDNA.txt", col_types = cols()) %>%
  select(sample=`Extraction No.`, mapped=`% Mapped`, mapped_q30=`% Mapped-Q30`, duplicates=Duplicates) %>%
  mutate(uniq_q30=mapped_q30*(1-duplicates))

# take the average of multiple library stats
dsDNA <- read_tsv("FA01p/Pig_allTo_13032017_clean.txt", col_types = cols()) %>%
  select(sample=`Extraction No.`, mapped=`% Mapped`, mapped_q30=`% Mapped-Q30`, duplicates=Duplicates) %>%
  group_by(sample) %>%
  summarise(mapped=mean(mapped), mapped_q30=mean(mapped_q30), duplicates=mean(duplicates)) %>%
  mutate(uniq_q30=mapped_q30*(1-duplicates))

# assign the new groups for the double-stranded DNA
dsDNA_libs <- groups %>%
  filter(type=="Lib") %>%
  inner_join(dsDNA, by="sample") %>%
  drop_na() %>%
  arrange(desc(uniq_q30)) %>%
  mutate(new_group=ceiling(row_number()/4)) %>%
  head(n=68)

last_group <- max(dsDNA_libs$new_group)

ssDNA_libs <- groups %>%
  filter(type=="Ext") %>%
  inner_join(ssDNA, by="sample") %>%
  drop_na() %>%
  arrange(desc(uniq_q30)) %>%
  mutate(new_group=last_group + ceiling(row_number()/4)) %>%
  head(n=28)

libraries <- bind_rows(dsDNA_libs, ssDNA_libs) %>%
  inner_join(meta, by="sample") 

# how many samples in each status grouping?
libraries %>%
  group_by(status) %>%
  tally()

# are there any index clashed in a capture group?
libraries %>%
  group_by(new_group, internal, external) %>%
  tally() %>%
  drop_na() %>%
  filter(n>1)

write_tsv(libraries, "FA01p/FA01p_capture_groups.tsv", na = "")

