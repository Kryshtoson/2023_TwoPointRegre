library(vegan)
library(tidyverse)
library(readxl)

list.files(r'(C:\Users\chytryk\OneDrive - MUNI\Devin\Jirka)')
brbl2num <- function(x){
  factor(x,
         levels = c('r', '+', '1', '2a', '2m', '2b', '3', '4', '5'),
         labels = c(.01, 1, 3, 10, 5, 20, 38, 68, 88)) %>%
  as.character() %>%
  as.numeric()}

read_xlsx(r'(C:\Users\chytryk\OneDrive - MUNI\Devin\Jirka\Jirkaplots1993-2021.xlsx)', skip = 10) %>%
  pivot_longer(-1) %>%
  mutate(value = brbl2num(value)) %>%
  separate(name, c('year', 'site')) %>%
  drop_na() %>%
  pivot_wider(values_fill = 0, names_from = taxon) -> spe_jirka

spe_init <- read_xlsx(r'(D:\Dropbox\PC (4)\Documents\R-project\2023\2023_OakBeech\data\KRIVOKLATSKO_Forest plots_oakwood and beechwood.xlsx)') %>%
  select(Plot, `Forest type`, year = Year, `Acer campestre_4`:`Viola x dubia_6`) %>%
  mutate(Plot = paste0(str_sub(Plot,1,1), '_', `Forest type`)) %>%
  select(-`Forest type`) %>%
  pivot_longer(`Acer campestre_4`:`Viola x dubia_6`) %>%
  separate(name, c('name', 'layer'), sep = '_') %>%
    filter(layer %in% c(6,7)) %>%
    select(-layer) %>%
  mutate(value = brbl2num(value)) %>%
  drop_na() %>%
  pivot_wider(values_fill = 0) %>%
  separate(Plot, c('plot', 'dataset'))

spe_oak <- spe_init %>% filter(dataset == 'Oakwood') %>% select(-dataset) %>% rename('site' = plot) %>% relocate(year)
spe_beech <- spe_init %>% filter(dataset == 'Beechwood') %>% select(-dataset) %>% rename('site' = plot) %>% relocate(year)

getwhatweneed <- function(spe, dataset, pa = F){
  if(pa == T){
    spex <- spe[-c(1:2)]
    spe[-c(1:2)] <- (spex != 0)}
  mx <- as_tibble(as.matrix(vegdist(spe[-c(1:2)])))
  names(mx) <- paste0(spe$year, '-', spe$site)
  mx[lower.tri(mx)] <- NA
  bind_cols(spe[1:2], mx) %>%
    pivot_longer(cols = -c(1:2)) %>%
    separate(name, c('year2', 'site2')) %>%
    filter(site == site2) %>%
    mutate(dataset = dataset,
           year = as.numeric(year))
}

map2(list(spe_oak, spe_beech, spe_jirka), list('Oakwood', 'Beechwood', 'Dry grassland'), ~getwhatweneed(.x, .y)) %>%
  bind_rows() -> tocomp

map2(list(spe_oak, spe_beech, spe_jirka), list('Oakwood', 'Beechwood', 'Dry grassland'), ~getwhatweneed(.x, .y, pa = T)) %>%
  bind_rows() -> tocomp_pa

tocomp_pa %>%
  group_by(dataset) %>%
  filter(year > 1995 & year2 > 1995) %>%
  filter(year == min(year)) %>%
  #filter(year2 != min(year)) %>%
  group_by(year2, dataset) %>%
  summarise(value = mean(value)) %>%
  ggplot(aes(as.numeric(year2), value)) +
  geom_line(aes(group = dataset, colour = dataset)) +
  theme_bw() +
  scale_x_continuous(breaks = 1996:2021, expand = c(0,0,0.01,0)) +
  scale_y_continuous(expand = c(0, 0, .1, 0)) +
  scale_colour_discrete(name = 'Dataset', breaks = c('Dry grassland', 'Oakwood', 'Beechwood')) +
  labs(y = 'Bray Curtis dissimilarity',
       x = 'Year',
  caption = 'Average difference from the initial year (1993)') +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(1,0),
        axis.title.x = element_blank(),
        legend.justification = c(1,0),
        legend.background = element_blank())

ggsave('Permanent_plots_change-in-time.png', height = 6, width = 10)


######################################################################################
# spp richness

getwhatweneed2 <- function(spe, dataset){
  spe[1:2] %>%
  mutate(spe_rich = rowSums(spe[-c(1:2)] != 0)) %>%
  filter(year > 1995)-> df
  df %>%
    left_join(df %>%
                filter(year == 1996) %>%
                select(site, spe_rich_ini = spe_rich)) %>%
    mutate(diff = spe_rich - spe_rich_ini) %>%
    group_by(year)%>%
    summarise(diff = mean(diff)) %>%
    mutate(dataset = dataset,
           year = as.numeric(year))
}

map2(list(spe_oak, spe_beech, spe_jirka), list('Oakwood', 'Beechwood', 'Dry grassland'), ~getwhatweneed2(.x, .y)) %>%
  bind_rows() -> tocomp_richness

getconfi <- function(type){
  x <- tocomp_richness$diff[tocomp_richness$dataset == type]
  t.test(x - mean(x))$conf.int[2] }


confis <- tibble(
  dataset = c('Dry grassland', 'Oakwood', 'Beechwood'),
  confi = c(getconfi('Dry grassland'), getconfi('Oakwood'), getconfi('Beechwood')))

tocomp_richness %>%
  ggplot(aes(year, diff, colour = dataset)) +
  geom_line(aes(group = dataset)) +
  labs(y = 'Species Richness', caption = 'Lines indicate confidence intervals.') +
  geom_hline(data = confis, aes(yintercept = confi, colour = dataset), linetype = 2, alpha = .5) +
  geom_hline(data = confis, aes(yintercept = -confi, colour = dataset), linetype = 2, alpha = .5) +
  scale_colour_discrete(name = 'Dataset', breaks = c('Dry grassland', 'Oakwood', 'Beechwood')) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(1,0),
        axis.title.x = element_blank(),
        legend.justification = c(1,0),
        legend.background = element_blank())

ggsave('Permanent_plots_richness_change-in-time.png', height = 6, width = 10)
