library(tidyverse)

fa01p <- read_tsv("FA01p/FA01p_endo-difference.txt")

fa01p.long <-fa01p %>%
    mutate(`% MT_reads`=(MT_reads/`Mapped Reads`)*100) %>%
    select(-c("Mapped Reads", "Total Reads", "MT_reads")) %>%
    pivot_longer(cols=c("% Mapped", "% Mapped-Q30", "% MT_reads", "Duplicates", "ReadLen (Mapped)")) %>%
    pivot_wider(id_cols=c("Extraction No.", "name"), names_from="Library", values_from="value") %>%
    drop_na()

fa01p.limits <- fa01p.long %>%
    mutate(value=pmax(ssDNA, dsDNA)) %>%
    group_by(name) %>%
    summarise(max=max(value))

ggplot(data=fa01p.long, aes(x=ssDNA, y=dsDNA)) +
    facet_wrap("name", scales="free") +
    geom_point() +
    geom_smooth(method=lm) +
    geom_abline(color="red", linetype="dashed") +
    geom_blank(data = fa01p.limits, mapping = aes(x=max, y=max)) +
    theme(aspect.ratio = 1)
