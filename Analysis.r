####  Moisture response curve experiment - script - Juan F. Dueñas juanfduenas@proton.me

# Load Packages we will work with.
pkgs <- c("tidyverse", "mgcv", "gratia", "ggpubr", "vegan", "phyloseq", "ggordiplots", "knitr")

vapply(pkgs, FUN = library, FUN.VALUE = logical(1L), logical.return = TRUE, character.only = TRUE)

# Figure 1 ------ Control soil

## load data, pivot and store on a list
rs <- read_delim("rs_db", delim = ",") %>% 
       filter(treatment%in%"Ctrl")%>%
       pivot_longer(
       cols = 7:10,
       names_to = "var",
       values_to = "y")%>%
       select(potID, dl, whc, no_poll, var, y)%>%
       rename(x=dl, z=whc)%>%
       group_split(var)%>% # split into several databases by treatment and store on a list
       set_names(., nm = c("cf", "mwd", "imp", "wsa")) #give names to each set  

# community data
f <- readRDS("RS_all.rds")%>% 
     subset_taxa(., Kingdom=="Fungi")%>% # eliminate non target organisms
     filter_taxa(., function(x) sum(x) > 0, TRUE)%>% # keep variants with more than 0 reads
     subset_samples(., treatment%in%c("ctrl")) # ignore warning

# alpha diversity
alf <- cbind(sample_data(f)[,c(1)], estimate_richness(f, measures=c("InvSimpson")))%>%
       mutate(var=rep("invsimp",44))%>%
       right_join(.,rs[["cf"]], by="potID")%>%
       select(-c(y, var.y))%>%
       rename(var=var.x, y=InvSimpson)%>%
       select(1,4,5,6,3,2) %>% as_tibble(.)

rs$alf <- alf # append to db list
rm(alf)

# beta diversity (will not go on the list)
betf <- vegdist(otu_table(f), method = "jac", binary = T)%>%
        metaMDS(., parallel=2, trace=F, weakties=T)

f_tr <- sample_data(f)[,c(1,3,4)] %>% data.frame(.) %>% rename(dl=moistureL, z=whc) %>%
        mutate_at('dl', as.factor)
        
ford <- gg_ordiplot(betf, groups = f_tr$dl, hull = F, label = F,
                    spiders = F, ellipse = F, plot = F, choices = c(1, 2), scaling=1)%>%
         pluck('df_ord')%>%
         rownames_to_column("potID")%>%
         mutate_at('potID', as.numeric)%>%
         left_join(.,f_tr, by="potID") %>% select(-Group)

rm(betf)
## Plots

# function for plots
Myplotsc=function(df, fml, laby){
           plot=list(ggplot()+
                     geom_point(data=df, aes(x=x, y=y, colour=z))+
                     geom_smooth(data=df, aes(x=x, y=y), method = "gam", formula = y ~ s(x, k=3), 
                            method.args=list(family=fml),
                            se=TRUE, color="Black")+
                     scale_x_continuous(breaks = c(1, 3,5,7,9,11,13,15))+
                     scale_colour_viridis_c(option="E", direction = -1)+
                     labs(y=laby, x="Dryness level", colour="WHC %")+
                     theme(text = element_text(size = 12),
                           panel.grid.minor = element_blank()))
}

fml <- list('gaussian', 'gaussian', 'poisson',
         'betar', 'gaussian')
laby <- list('Log of copy number per g of soil', 'Mean weight diameter (mm)', 'Water drop penetration time (sec.)',
          'Fraction of water stable aggregates', 'Inverse Simpson index (ASVs)')

pl <- mapply(Myplotsc, rs, fml, laby)

ord <- ggplot() + 
            geom_vline(xintercept=0.0, color="White", linewidth=1, linetype=1)+
            geom_hline(yintercept=0.0, color="White", linewidth=1, linetype=1)+    
            geom_point(data=ford, aes(x=x, y=y, color=z), alpha=0.8,  size=2, show.legend = T) +
            scale_colour_viridis_c(option = "E", direction = -1)+
            labs(x="Axis 1", y="Axis 2", color="WHC %")+
            theme(axis.title = element_text(size = 12),
                  legend.key = element_blank(),  #removes the box around each legend item
                  legend.position = "right", #legend at the bottom
                  legend.text = element_text(size=12),
                  panel.border = element_rect(colour = "Black", fill = F),
                  panel.grid = element_blank())        

pl$ord <- ord # append ordination plot to list

(fig1 <- ggarrange(pl$mwd, pl$wsa, pl$imp, pl$cf, pl$alf, pl$ord, ncol = 3, nrow = 2,
                  common.legend = T, legend = "right", labels = c('a', 'b', 'c', 'd', 'e', 'f')))

# function for GAMS

Mymodsc=function(df, fml){
  plot=list(mod=gam(data = df, formula = y~s(x, k=3), family = fml, method="REML"))
}

mds <- mapply(Mymodsc, rs, fml)

sm <- lapply(mds, summary)%>% # get model summary
      lapply(., '[[', 's.table')%>% # extract s.table
      imap_dfr(., ~ bind_cols(mod = .y,  edf = .x[,1], ref.df = .x[,2], 
                              Fval = round(.x[,3], 3), pval = round(.x[,4], 3))) # format stats

# check diagnostics 
#overview(mds$alf.mod)
check <- function(b, k.sample = 5000, k.rep = 200) {
  mgcv:::k.check(b, subsample = k.sample, n.rep = k.rep)
}

chs <- lapply(mds, check)%>%
       imap_dfr(., ~ bind_cols(mod = .y,  k = .x[,1], edf = .x[,2], 
                          kindex = round(.x[,3], 3), pval = round(.x[,4], 3))) # format stats

appraise(mds$cf.mod, method = "simulate", n_simulate = 1000)
appraise(mds$mwd.mod, method = "simulate", n_simulate = 1000)
appraise(mds$imp.mod, method = "simulate", n_simulate = 1000)
appraise(mds$wsa.mod, method = "simulate", n_simulate = 1000)
appraise(mds$alf.mod, method = "simulate", n_simulate = 1000)

#kable(sm, format = "simple", caption = "Table 1", digits = 3)
write_csv(sm, "table1")

#Multivariate
adonis2(otu_table(f)~f_tr$dl, permutations = 999, method = 'jac', binary=T, by="margin")
