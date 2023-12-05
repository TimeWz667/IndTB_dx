library(tidyverse)


targets <- read_csv(here::here("data", "targets.csv"))


targets %>% 
  filter(Year >= 2017) %>% 
  filter(Index %in% c("IncR", "CNR", "TestR", "PrevUt", "PrAsym", "PrSym", "PrExCS", 
                      "PrInc", "PrDR", "PrDetPub", "DrugTime", "TxI")) %>% 
  arrange(Index, Year) %>% 
  select(Index, Year, Tag, M, L, U, N, Error, Std) %>% 
  write_csv(here::here("data", "targets_sel.csv"))
