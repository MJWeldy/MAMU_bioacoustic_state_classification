# read file path
all_paths <-
  list.files(path = "E://Scores//Combined",
             pattern = "*.csv",
             full.names = TRUE)

# read file content
all_content <-
  all_paths %>%
  lapply(fread)

# read file name
all_filenames <- all_paths %>%
  basename() %>%
  as.list()

# combine file content list and file name list
all_lists <- mapply(c, all_filenames,  all_content,SIMPLIFY = FALSE)

# unlist all lists and change column name
all_result <- rbindlist(all_lists, fill = T)
names(all_result)[1:2] <- c("site","site_rows")
tbl_plot <- all_result %>%
  dplyr::select(site, date, Offset, BRMA) %>% 
  mutate(week = week(date)) %>% 
  group_by(site,week) %>% 
  summarise(n = sum(BRMA>0.9))

#For testing with HMM
HMM_df <- tbl_plot %>% 
  drop_na(week) %>% 
  mutate(n_week = as.numeric(week)) %>% 
  #unique(week)
  #group_by(site) %>% 
  pivot_wider(site,names_from = n_week, values_from = n, names_sort = TRUE) 
#write.csv(HMM_df, "./data/clean/full_df.csv")

site_states <- read.csv("./data/raw/site_state_observations.csv")
site_states$response <- ifelse(site_states$True_state < 1, NA,
                               ifelse(site_states$True_state>3,3,site_states$True_state))