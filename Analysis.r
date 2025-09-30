####  Moisture response curve experiment - script - Juan F. Dueñas juanfduenas@proton.me

# Load Packages we will work with.
pkgs <- c("tidyverse", "mgcv", "gratia", "ggpubr", "vegan", "phyloseq", "ggordiplots", "knitr")

vapply(pkgs, FUN = library, FUN.VALUE = logical(1L), logical.return = TRUE, character.only = TRUE)
#

# Figure 1 and table ---------- response curves control soil ####

## load data, pivot and store on a list 
rs <- read_delim("rs_db.csv", delim = ",") %>% 
       filter(treatment%in%"Ctrl")%>%
       pivot_longer(
       cols = 7:10,
       names_to = "var",
       values_to = "y")%>%
       select(potID, dl, whc, no_poll, var, y)%>%
       rename(x=dl, z=whc)%>%
       filter(!is.na(y))%>%
       group_split(var)%>% # split into several databases by treatment and store on a list
       set_names(., nm = c("cf", "mwd", "imp", "wsa"))%>% #give names to each set  
       as.list(.)

# community data
f <- readRDS("RS_all.rds")%>% 
     subset_taxa(., Kingdom=="Fungi")%>% # eliminate non target organisms
     subset_samples(., treatment%in%c("ctrl")) # ignore warning

# alpha diversity
alf <- cbind(sample_data(f)[,c(1)], estimate_richness(f, measures=c("InvSimpson")))%>%
       mutate(var=rep("invsimp",44))%>%
       right_join(.,rs[["cf"]], by="potID")%>% #fixed 
       select(-c(y, var.y))%>%
       rename(var=var.x, y=InvSimpson)%>%
       select(1,4,5,6,3,2) %>% 
       filter(!is.na(y)) %>% as_tibble(.) %>%
       mutate_at('potID', as.integer) 

rs$alf <- alf # append to db list
rm(alf)

# beta diversity
betf <- vegdist(otu_table(f), method = "jac", binary = T)%>%
  metaMDS(., parallel=2, trace=F, weakties=T)

f_tr <- sample_data(f)[,c(1,3,4)] %>% data.frame(.) %>% rename(dl=moistureL, z=whc) %>%
  mutate_at('dl', as.factor)

# Model selection routine

Mymodselc <- function(df, fml) {
  # function to set two pairs of models one linear one gam  
  set1 <- list(
                 gm     = gam(y ~ s(x, bs = "tp", k = 3), data = df, family = fml),
                 linear = lm(y ~ x, data = df)
                )
}
fml <- list('gaussian', 'gaussian', 'poisson',
            'betar', 'gaussian')
models <- Map(Mymodselc, rs, fml) # iteration of model fit

# function to extract AIC and rank models
extract_model_selection_table <- function(model_list) {
        out <- lapply(names(model_list), function(ds) {
        sublist <- model_list[[ds]]
    
         # extract AICs
          aics <- sapply(sublist, AIC)
    
         # ΔAIC
          delta <- aics - min(aics)
    
         # Akaike weights
         weights <- exp(-0.5 * delta)
         weights <- weights / sum(weights)
    
        # build selection table
        tab <- data.frame(
               Dataset = ds,
               Model = names(aics),
               AIC = aics,
               DeltaAIC = delta,
               Weight = weights,
               Rank = rank(aics, ties.method = "first")
              )
    
        tab[order(tab$AIC), ]
  })
  
  # combine into one big table
  do.call(rbind, out)
}


myAICs <- extract_model_selection_table(models)

# function for Models (GAMS)
Mygamsc=function(df, fml){
  # fits a GAM to each dataset on a list in df with max 3 base functions
  # returns a list of named gams.
  # argument fml is a string list stating the family of the conditional prob. and must be provided separately
  md=list(gm=gam(data = df, formula = y~s(x, bs="tp", k=3), family = fml, method="REML"))
}

fml <- list('gaussian', 'gaussian', 'poisson',
            'betar', 'gaussian')

mds <- mapply(Mygamsc, rs, fml) # iteratively apply Mygamsc

sm <- lapply(mds, summary)%>% # get each model summary
      lapply(., '[[', 's.table')%>% # extract each s.table
      imap_dfr(., ~ bind_cols(mod = .y,  edf = .x[,1], ref.df = .x[,2], 
                          Fval = round(.x[,3], 3), pval = round(.x[,4], 3))) # format stats

#tables
kable(sm, format = "simple", caption = "Table 1", digits = 3) # pretty simple table
write_csv(sm, "table1") # good old way to save results
# check model fit with diagnostics 
#overview(mds$alf.mod)
check <- function(b, k.sample = 5000, k.rep = 200) {
  mgcv:::k.check(b, subsample = k.sample, n.rep = k.rep)
}

# function to check if k is correctly specified
chs <- lapply(mds, check)%>%
  imap_dfr(., ~ bind_cols(mod = .y,  k = .x[,1], edf = .x[,2], 
                          kindex = round(.x[,3], 3), pval = round(.x[,4], 3))) # format stats

appraise(mds$cf.gm, method = "simulate", n_simulate = 1000)
appraise(mds$mwd.gm, method = "simulate", n_simulate = 1000)
appraise(mds$imp.gm, method = "simulate", n_simulate = 1000)
appraise(mds$wsa.gm, method = "simulate", n_simulate = 1000)
appraise(mds$alf.gm, method = "simulate", n_simulate = 1000)

# Multivariate Model (PERMANOVA)
adonis2(otu_table(f)~f_tr$dl, permutations = 9999, method = 'jac', binary=T, by="term")

## Plots
# estimate uncertainty around fits
Myse <- function(mds, df, i){
         # function to extract the fit and standard error in the response scale
         # models are the list of models defined by Mymodsc
         # df is the list of data.frames that contain the predictor needed to pass to function predict
         # i is the number of iterations of this routine = number of datasets/models
         nd <- lapply(df, function(y) with(y, expand.grid(x = seq(min(x),max(x), length.out=100))))
         fit <- lapply(mds, function(x){do.call(cbind, predict(x, nd[[i]], type='response', se.fit=T, unconditional=T))})
         fit2 <- lapply(fit, function(x){cbind(x, nd[[i]])})
}
se <- Myse(mds=mds, df=rs, i=5)

# function for plots
Myplotsc=function(df, fml, laby, se){
           plot=list(ggplot()+
                     geom_point(data=df, aes(x=x, y=y, colour=z))+
                     geom_line(data=se, aes(x=x, y=fit), inherit.aes = F)+ # use fit estimated by predict
                     geom_ribbon(data=se,                                  # estimate alpha=95% Credible Intervals
                                   aes(x=x, ymin = fit- 1.96*(se.fit), ymax = fit + 1.96*(se.fit)), 
                                   alpha = 0.5, fill = "Grey", inherit.aes = F)+ # this constant comes from qnorm((1-0.95)/2, lower.tail = F)
                     scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15))+
                     scale_colour_viridis_c(option="E", direction = -1)+
                     labs(y=laby, x="Dryness level", colour="WHC %")+
                     theme_bw()+
                     theme(text = element_text(size = 12),
                           panel.grid.minor = element_blank()))
}


laby <- list('Log of copy number per g of soil', 'Mean weight diameter (mm)', 'Water drop penetration time (sec.)',
          'Fraction of water stable aggregates', 'Inverse Simpson index (ASVs)')

pl <- mapply(Myplotsc, rs, fml, laby, se)


#plot ord

# beta diversity (will not go on the list)
ford <- gg_ordiplot(betf, groups = f_tr$dl, hull = F, label = F,
                    spiders = F, ellipse = F, plot = F, choices = c(1, 2), scaling=1)%>%
        pluck('df_ord')%>%
        rownames_to_column("potID")%>%
        mutate_at('potID', as.numeric)%>%
        left_join(.,f_tr, by="potID") %>% select(-Group)

rm(betf)

ord <- ggplot() + 
            geom_vline(xintercept=0.0, color="Grey", linewidth=1, linetype=1)+
            geom_hline(yintercept=0.0, color="Grey", linewidth=1, linetype=1)+    
            geom_point(data=ford, aes(x=x, y=y, color=z), alpha=0.8,  size=2, show.legend = T) +
            scale_colour_viridis_c(option = "E", direction = -1)+
            labs(x="Axis 1", y="Axis 2", color="WHC %")+
            theme_bw()+
            theme(axis.title = element_text(size = 12),
                  legend.key = element_blank(),  #removes the box around each legend item
                  legend.position = "right", #legend at the bottom
                  legend.text = element_text(size=12),
                  panel.border = element_rect(colour = "Black", fill = F),
                  panel.grid = element_blank())        

pl$ord <- ord # append ordination plot to list

(fig1 <- ggarrange(pl$mwd, pl$wsa, pl$imp, pl$cf, pl$alf, pl$ord, ncol = 3, nrow = 2,
                  common.legend = T, legend = "right", labels = c('a', 'b', 'c', 'd', 'e', 'f')))

# Figure 2 -------------------- Thresholds control soil  ####

# load dataset and select only variables that are responding to the water gradient
rs <- read_delim("rs_db.csv", delim = ",") %>% 
      filter(treatment%in%"Ctrl")%>%
      select(c(1:7,9))%>%
      pivot_longer(
      cols = 7:8,
      names_to = "var",
      values_to = "y")%>%
      select(potID, dl, whc, no_poll, var, y)%>%
      rename(x=dl, z=whc, n=no_poll)%>%
      filter(!is.na(y))%>%
      group_split(var)%>% # split into several databases by treatment and store on a list
      set_names(., nm = c("mwd", "imp"))%>% #give names to each set  
      as.list(.)


#1)Find fits and calculate AICs - the conditional dist. can only be changed to binomial if the response variable requires
# otherwise one has to transform
lbs <- c("chngpt", "boot")
vapply(lbs, FUN = library, FUN.VALUE = logical(1L), logical.return = TRUE, character.only = TRUE)
      
Mymodels=function(df, fml="gaussian"){
  models=list(
    linear=gam(data = df, formula=y ~ x, family = fml, method="REML"),
    gm=gam(data = df, formula = y ~ s(x, k=3), family = fml, method="REML"), 
    stegmented=chngptm(formula.1 = y~1, formula.2 = ~x, data = df, type="stegmented", family = fml, REML = T), # removed argument: #, chngpt.init = 0.7#
    stepm=chngptm(formula.1 = y~1, formula.2 = ~x, data = df, type="step", family = fml, REML = T),
    segmented=chngptm(formula.1 = y~1, formula.2 = ~x, data = df, type="segmented", family = fml, REML = T))
}

#apply function to all the sets within our data list
ths <- lapply(rs, Mymodels)

MyAICs=function(models){
  AICs=with(models,c(AICLin=AIC(linear),
                     AICgam=gm$aic,
                     AICSteg=stegmented$best.fit$aic,
                     AICStep=stepm$best.fit$aic,
                     AIC.segmented=segmented$best.fit$aic))
}

aics <- t(round(do.call(rbind, lapply(ths, MyAICs)),2))%>% as_tibble(rownames=NA) %>% 
        rownames_to_column() %>%
        mutate(Delta_mwd = mwd - min(mwd),
               Delta_imp = imp - min(imp, na.rm = TRUE)) # Stegmented is the best simplification in all cases

# visualize best models
plot(ths$mwd$stegmented)
plot(ths$imp$stegmented)

# Function to bootstrap thresholds - get an idea of uncertainty in threshold determination
bs=function(data, indices, fml="gaussian"){
            d=data[indices,] #allows boot to select sample
           models=list(stegmented=chngptm(formula.1 = y~1,formula.2 = ~x, data = d, type="stegmented", family = fml, 
                                 REML = T))
           OUT= c(thresSteg=models$stegmented$chngpt)
} #If data are qualitative, do family="binomial"

# It takes a few seconds with multicores 
results=lapply(rs, boot, statistic = bs, R=200, parallel = "multicore", ncpus = 4)

# Get the vector of thresholds estimated in each case 
dfboot=lapply(results, "[[", "t")%>%
       set_names(., nm = c('mwd', 'imp')) #give names to each database
      new_colnames=c("thresSteg")
      dfboot <- lapply(dfboot, 
                 function(x) {colnames(x) <- new_colnames; x})

#3) Validate thresholds bootstrap of slope and intcp changes
# bootthres=function(data,indices,formula=y~x,thres=thres){
#   d=data[indices,] #allows boot to select sample
#   mdl=lm(data=d,formula = formula)
#   slp=coef(mdl)[2]
#   intcp=coef(mdl)[1]
#   thy=predict.lm(mdl,newdata = data.frame(x=thres))
#   return(c(slp,intcp,thy))
# }
# 
# funcdiff=function(dff,thres,variable,bootthres=bootthres){
#   dfs=list(before=dff[dff$x<=thres,],
#            after=dff[dff$x>thres,])
#   
#   mdls=lapply(dfs, function(x) return(boot(data = x,statistic = bootthres,R=200,thres=thres)))
#   
#   lboot=lapply(mdls, function(x) {
#     res=as.data.frame(x$t)
#     colnames(res)=c("slope","intcp","value")
#     return(res)
#   })
# }  
# 
# #Validating thresholds - in Berdugo et al these values are compared with Mann-Whitney U tests
# resdf.ctrl=funcdiff(rs[["mwd"]],median(dfboot[["mwd"]][,1]),variable,bootthres = bootthres) 

# #test if median difference between slopes or intercept before of after the threshold is 0 (H0)
# # use Mann-Whitney U test - unpaired
# wilcox.test(resdf.ctrl$before$slope, resdf.ctrl$after$slope, paired = F) # significant differences in slope were found
# wilcox.test(resdf.ctrl$before$intcp, resdf.ctrl$after$intcp, paired = F) # idem for intercepts
# wilcox.test(resdf.ctrl$before$value, resdf.ctrl$after$value, paired = F) # idem for predicted values

# Plots
th.mwd <- data.frame(cbind(dfboot[["mwd"]][,1],c(rep("3", 200))))%>% # in almost all treatments stegmented was the best solution
                           rename(th=X1, var=X2)%>%
                           mutate_at('th', as.numeric)%>%
                           mutate_at('var', as.integer)

th.imp <- data.frame(cbind(dfboot[["imp"]][,1],c(rep("3", 200))))%>% # in almost all treatments stegmented was the best solution
                          rename(th=X1, var=X2)%>%
                          mutate_at('th', as.numeric)%>%
                          mutate_at('var', as.integer)


# function for Models (GAMS)
Mygamsc=function(df, fml){
  # fits a GAM to each dataset on a list in df with max 3 base functions
  # returns a list of named gams.
  # argument fml is a string list stating the family of the conditional prob. and must be provided separately
  md=list(gm=gam(data = df, formula = y~s(x, bs="tp", k=3), family = fml, method="REML"))
}

fml <- list('gaussian', 'poisson')

mds <- mapply(Mygamsc, rs, fml) # iteratively apply Mygamsc

Myse <- function(mds, df, i){
  # function to extract the fit and standard error in the response scale
  # models are the list of models defined by Mymodsc
  # df is the list of data.frames that contain the predictor needed to pass to function predict
  # i is the number of iterations of this routine = number of datasets/models
  nd <- lapply(df, function(y) with(y, expand.grid(x = seq(min(x),max(x), length.out=100))))
  fit <- lapply(mds, function(x){do.call(cbind, predict(x, nd[[i]], type='response', se.fit=T, unconditional=T))})
  fit2 <- lapply(fit, function(x){cbind(x, nd[[i]])})
}
se <- Myse(mds=mds, df=rs, i=2)

# Plot both trends in same plot
colrs <- RColorBrewer::brewer.pal(3, "Set2")

(p <-        ggplot()+
             geom_density(data=th.mwd, aes(x=th), colour=colrs[3], fill="Grey", alpha=0.4)+
             geom_segment(data=th.mwd, aes(x = median(th), xend = median(th), y = 0, yend = 0.8), colour=colrs[3], linetype=2, linewidth = 1)+
             geom_density(data=th.imp, aes(x=th), colour=colrs[2], fill="Grey", alpha=0.4)+
             geom_segment(data=th.imp, aes(x = median(th), xend = median(th), y = 0, yend = 0.8), colour=colrs[2], linetype=2, linewidth = 1)+
             geom_line(data=se[[1]], aes(x=x, y=fit), inherit.aes = F)+ # use fit estimated by predict
             geom_ribbon(data=se[[1]],                                  # estimate alpha=95% Credible Intervals
                          aes(x=x, ymin = fit- 1.96*(se.fit), ymax = fit + 1.96*(se.fit)), 
                          alpha = 0.5, fill = "Grey", inherit.aes = F)+ # this constant comes from qnorm((1-0.95)/2, lower.tail = F)
             geom_line(data=se[[2]], aes(x=x, y=fit/10), inherit.aes = F)+ # use fit estimated by predict
             geom_ribbon(data=se[[2]],                                  # estimate alpha=95% Credible Intervals
                aes(x=x, ymin = fit/10- 1.96*(se.fit/10), ymax = fit/10 + 1.96*(se.fit/10)), 
                alpha = 0.5, fill = "Grey", inherit.aes = F)+ # t
             scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15))+
             scale_y_continuous(name = 'Mean weight diameter (mm)', sec.axis = sec_axis(transform=~./10, name='Water drop penetration time (sec./10)'))+
             labs(x="Dryness level")+
             theme_bw()+
             theme(text = element_text(size = 12),
                    panel.grid.minor = element_blank()))
# Figure S1, and tables ------- Fungal communities in control ####

## load data, pivot and store on a list 
rs <- read_delim("rs_db.csv", delim = ",") %>% 
  filter(treatment%in%"Ctrl")%>%
  pivot_longer(
    cols = 7:10,
    names_to = "var",
    values_to = "y")%>%
  select(potID, dl, whc, no_poll, var, y)%>%
  rename(x=dl, z=whc)%>%
  filter(!is.na(y))%>%
  group_split(var)%>% # split into several databases by treatment and store on a list
  set_names(., nm = c("cf", "mwd", "imp", "wsa"))%>% #give names to each set  
  as.list(.)

# community data
f <- readRDS("RS_all.rds")%>% 
     subset_taxa(., Kingdom=="Fungi")%>% # eliminate non target organisms
     filter_taxa(., function(x) sum(x) > 0, TRUE)%>% # keep variants with more than 0 reads
     subset_samples(., treatment%in%c("ctrl")) # ignore warning

#subset per Phyla
a=subset_taxa(f, Phylum=="Ascomycota")
b=subset_taxa(f, Phylum=="Basidiomycota")
c=subset_taxa(f, Phylum=="Mortierellomycota")
d=subset_taxa(f, Phylum=="Mucoromycota")

# alpha diversity per phyla

alf.a <- cbind(sample_data(a)[,c(1)], estimate_richness(a, measures=c("InvSimpson")))%>%rename(Ascomycota=InvSimpson)
alf.b <- cbind(sample_data(b)[,c(1)], estimate_richness(b, measures=c("InvSimpson")))%>%rename(Basidiomycota=InvSimpson)
alf.c <- cbind(sample_data(c)[,c(1)], estimate_richness(c, measures=c("InvSimpson")))%>%rename(Mortierellomycota=InvSimpson)
alf.d <- cbind(sample_data(d)[,c(1)], estimate_richness(d, measures=c("InvSimpson")))%>%rename(Mucoromycota=InvSimpson)
tl <- list(alf.a, alf.b, alf.c, alf.d) 

alf_p <- reduce(tl, full_join, by="potID")%>%
  right_join(.,rs[["cf"]], by="potID")%>% 
  select(1:8)%>%
  pivot_longer(
    cols = 2:5,
    names_to = "p",
    values_to = "y")%>%
  filter(!is.na(y)) %>%
  mutate_at('p', as.factor)%>%
  as_tibble(.) 

rm(alf.a, alf.b, alf.c, alf.d, tl)

# beta diversity - relative abundance phylum
ra_p = psmelt(f) %>% select(potID, Abundance, moistureL, whc, Phylum) %>% #melt and select columns that I am interested in
  group_by(potID, moistureL, whc, Phylum) %>%
  summarize(abund = sum(Abundance)) %>%
  mutate(ra = abund / sum(abund)) %>%
  filter(ra>0.05) %>%
  ungroup(.)%>% rename(x=moistureL, z=whc, p=Phylum, y=ra) %>% 
  mutate_at('p', as.factor) %>% as_tibble(.)

# beta diversity - relative abundance guild
# ra_g = psmelt(f) %>% select(potID, Abundance, moistureL, whc, Guild) %>% #melt and select columns that I am interested in
#                      group_by(potID, moistureL, whc, Guild) %>%
#                      summarize(abund = sum(Abundance)) %>%
#                      mutate(ra = abund / sum(abund)) %>%
#                      ungroup()%>% as_tibble()%>%
#                      rename(x=moistureL, z=whc, p=Guild, y=ra) %>%
#                      mutate_at('p', as.factor)
cd <- list(alf=alf_p, bet=ra_p)
rm(ra_p, alf_p, a, b, c, d)           
# function for Models (GAMS)

Mygamsc=function(df, fml, k){
  # fits a GAM to each dataset on a list in df with max 3 base functions
  # returns a list of named gams.
  # argument fml is a string list stating the family of the conditional prob. and must be provided separately
  md=list(gm=gam(data = df, formula = y ~ p + s(x, by=p, k=k), family = fml, method="REML", select = T))
}

fml <- list('Gamma', 'betar')

gmsp <- mapply(Mygamsc, cd, fml, 3) # iteratively apply Mygamsc

sm <- lapply(gmsp, summary)%>% # get each model summary
  lapply(., '[[', 's.table')%>% # extract each s.table
  imap_dfr(., ~ bind_cols(mod = .y,  edf = .x[,1], ref.df = .x[,2], 
                          Fval = round(.x[,3], 3), pval = round(.x[,4], 3))) # format stats

pm <- lapply(gmsp, summary)%>% # get each model summary
  lapply(., '[[', 'p.table')%>%
  imap_dfr(., ~ bind_cols(mod = .y, term = rownames(.x), bet = round(.x[,1], 3), StdE = round(.x[,2],3), 
                          tval = round(.x[,3], 3), pval = round(.x[,4], 3)))#%>% # format stats

# Tables
kable(sm, format = "simple", caption = "Table 2", digits = 3) # pretty simple table
kable(pm, format = "simple", caption = "Table 3", digits = 3) # pretty simple table
write_csv(sm, "table2") # good old way to save results
write_csv(pm, "table3") # good old way to save results

# check model fit with diagnostics 
appraise(gmsp$alf.gm, method = "simulate", n_simulate = 1000)
appraise(gmsp$bet.gm, method = "simulate", n_simulate = 1000)
#appraise(gmsp$gld.gm, method = "simulate", n_simulate = 1000)

# Plots
# estimate uncertainty around fits
Myse <- function(mds, df, i){
  # function to extract the fit and standard error in the response scale
  # models are the list of models defined by Mymodsc
  # df is the list of data.frames that contain the predictor needed to pass to function predict
  # i is the number of iterations of this routine = number of datasets/models
  nd <- lapply(df, function(y) with(y, expand.grid(x = evenly(x, n=100), p=levels(p))))
  fit <- lapply(mds, function(x){do.call(tibble, predict(x, nd[[i]], type='response', se.fit=T, unconditional=T))})
  fit2 <- lapply(fit, function(x){cbind(x, nd[[i]])})
}
se <- Myse(mds=gmsp, df=cd, i=2)

# function for plots
Myplotsc=function(df, fml, laby, se){
  clrs <- viridis::viridis(3, option = "E", direction = -1)  
  plot=list(ggplot(data=df, aes(x=x, y=y, colour = z))+
              annotate("rect", xmin=0.9,xmax=5.9,ymin=0,ymax=Inf, fill=alpha(clrs[3], 0.2), colour=NA)+
              annotate("rect", xmin=5.9,xmax=9.1,ymin=0,ymax=Inf, fill=alpha(clrs[2], 0.2), colour=NA)+
              annotate("rect", xmin=9.1,xmax=15.2,ymin=0,ymax=Inf, fill=alpha(clrs[1], 0.2), colour=NA)+
              geom_point(alpha=0.5)+
              geom_line(data=se, aes(x=x, y=fit, linetype = p), linewidth=0.5, inherit.aes = F)+ # use fit estimated by predict
              geom_ribbon(data=se,                                  # estimate alpha=95% Credible Intervals
                          aes(x=x, ymin = fit- 1.96*(se.fit), ymax = fit + 1.96*(se.fit), group = p), 
                          alpha = 0.5, fill = "Grey", inherit.aes = F)+ # this constant comes from qnorm((1-0.95)/2, lower.tail = F)
              scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15))+
              #scale_color_manual(values=colrs)+
              scale_colour_viridis_c(option = "E", direction = -1)+
              labs(y=laby, x="Dryness level", colour="WHC %", linetype="Group ID")+
              theme_bw()+
              theme(text = element_text(size = 12),
                    panel.grid.minor = element_blank()))
}

laby <- list('Inverse Simpson index (ASVs)', 'Relative abundance of reads')

pl <- mapply(Myplotsc, cd, fml, laby, se)

(figs1 <- ggarrange(pl$alf, pl$bet, ncol = 2, nrow = 1,
                    common.legend = T, legend = "right", labels = c('a', 'b')))


# Figure 3, figs2 and tables -- Single factor effects ####

'%notin%' <- Negate('%in%') # useful custom function to filter out

# load data and select all treatments except for control

rs <- read_delim("rs_db.csv", delim = ",") %>% 
      filter(treatment %notin% "All_in") %>% 
      mutate(pha=case_when(dl < 6 ~ "Wet",
                       dl < 10  ~ "Medium",
                       dl >= 10 ~ "Dry" )) %>%
      mutate(pha=fct_relevel(pha, c("Wet", "Medium", "Dry")))%>%
      mutate_at(c('no_poll', "treatment"), as.factor)%>%
      mutate(treatment=fct_relevel(treatment, c("Ctrl", "Copper","µPlastic", "Nitrogen", "Salinity", "Surfactant")))%>%   
      pivot_longer(
      cols = 7:10,
      names_to = "var",
      values_to = "y")%>%
      select(potID, dl, whc, treatment, no_poll, pha, var, y)%>%
      rename(x=dl, z=whc, n=no_poll, t=treatment, p=pha)%>%
      filter(!is.na(y))%>%
      group_split(var)%>% # split into several databases by treatment and store on a list
      set_names(., nm = c("cf", "mwd", "imp", "wsa"))%>% #give names to each set  
      as.list(.)

# community data
f <- readRDS("RS_all.rds")%>% 
     subset_taxa(., Kingdom=="Fungi")%>% # eliminate non target organisms
     subset_samples(., treatment%notin%c("all"))

f <-  subset_samples(f, sample_names(f) %notin% c("294","295", "285", "304", "313"))%>% # Eliminate samples with ultra low read numbers (<1000)
     filter_taxa(., function(x) sum(x) > 0, TRUE) # keep variants with more than 0 reads    

# alpha diversity
alf <- cbind(sample_data(f)[,c(1)], estimate_richness(f, measures=c("InvSimpson")))%>%
       mutate(var=rep("invsimp",250))%>%
       right_join(.,rs[["cf"]], by="potID")%>% #fixed 
       select(-c(y, var.y))%>%
       rename(var=var.x, y=InvSimpson)%>%
       select(1,4,5,6,7,8,3,2) %>%
       filter(!is.na(y)) %>% as_tibble(.) %>%
       mutate_at('potID', as.integer) 

rs$alf <- alf # append to db list
rm(alf)

# beta diversity
betf <- vegdist(otu_table(f), method = "jac", binary = T)%>%
        metaMDS(., parallel=2, trace=F, weakties=T)

f_tr <- sample_data(f)[,c(1,2,3,4)] %>% data.frame(.) %>% rename(t=treatment, m=moistureL, z=whc) %>%
        mutate_at('t', as.factor)%>%
        mutate(p=case_when(m < 6 ~ "Wet",
                           m < 10  ~ "Medium",
                           m >= 10 ~ "Dry" )) %>%
        mutate(p=fct_relevel(p, c("Wet", "Medium", "Dry")))%>%
        pivot_wider(names_from = t, values_from = t, values_fn = ~1, values_fill = 0 )

# Model selection

Mymodselc <- function(df, fml) {
  # function to set two pairs of models one linear one gam  
  gms <- list(
    gm1 = gam(data = df, formula = y ~ t + s(x, by=t, bs = 'ts', k=3), family = fml, method="REML", select = T),
    gm2 = gam(data = df, formula = y ~ s(x, by=t, bs = 'ts', k=3), family = fml, method="REML", select = T)
  )
}

fml <- list('gaussian', 'Gamma', 'poisson',
            'betar', 'gaussian')
gmssel <- Map(Mymodselc, rs, fml) # iteration of model fit

# function to extract AIC and rank models
extract_model_selection_table <- function(model_list) {
  out <- lapply(names(model_list), function(ds) {
    sublist <- model_list[[ds]]
    
    # extract AICs
    aics <- sapply(sublist, AIC)
    
    # ΔAIC
    delta <- aics - min(aics)
    
    # Akaike weights
    weights <- exp(-0.5 * delta)
    weights <- weights / sum(weights)
    
    # build selection table
    tab <- data.frame(
      Dataset = ds,
      Model = names(aics),
      AIC = aics,
      DeltaAIC = delta,
      Weight = weights,
      Rank = rank(aics, ties.method = "first")
    )
    
    tab[order(tab$AIC), ]
  })
  
  # combine into one big table
  do.call(rbind, out)
}

myAICs <- extract_model_selection_table(gmssel)
# differences are minimal, chose GAM formula 1

# function for Models with interactions (GAMS)

Mygamsc=function(df, fml, k){
  # fits a GAM to each dataset on a list in df with max 3 base functions
  # returns a list of named gams.
  # argument fml is a string list stating the family of the conditional prob. and must be provided separately
  md=list(gm=gam(data = df, formula = y ~ t + s(x, by=t, k=k), family = fml, method="REML", select = T))
}

fml <- list('gaussian', 'Gamma', 'poisson',
            'betar', 'gaussian')

gms <- mapply(Mygamsc, rs, fml, 3) # iteratively apply Mygamsc

sm <- lapply(gms, summary)%>% # get each model summary
      lapply(., '[[', 's.table')%>% # extract each s.table
      imap_dfr(., ~ bind_cols(mod = .y,  edf = .x[,1], ref.df = .x[,2], 
                          Fval = round(.x[,3], 3), pval = round(.x[,4], 3)))%>% # format stats
      mutate(term=rep(c("Ctrl", "Copper","µPlastic", "Nitrogen", "Salinity", "Surfactant"),5)) 

pm <- lapply(gms, summary)%>% # get each model summary parametric part
      lapply(., '[[', 'p.table')%>%
      imap_dfr(., ~ bind_cols(mod = .y, term = rownames(.x), bet = round(.x[,1], 3), StdE = round(.x[,2],3), 
                          tval = round(.x[,3], 3), pval = round(.x[,4], 3)))#%>% # format stats

# Tables
kable(sm, format = "simple", caption = "Table 4", digits = 3) # pretty simple table
kable(pm, format = "simple", caption = "Table 5", digits = 3) # pretty simple table
write_csv(sm, "table4") # good old way to save results
write_csv(pm, "table5") # good old way to save results
# check model fit with diagnostics 
#overview(mds$alf.mod)
check <- function(b, k.sample = 5000, k.rep = 200) {
  mgcv:::k.check(b, subsample = k.sample, n.rep = k.rep)
}

# function to check if k is correctly specified
chs <- lapply(gms, check)%>%
  imap_dfr(., ~ bind_cols(mod = .y,  k = .x[,1], edf = .x[,2], 
                          kindex = round(.x[,3], 3), pval = round(.x[,4], 3))) # format stats

appraise(gms$cf.gm, method = "simulate", n_simulate = 1000)
appraise(gms$mwd.gm, method = "simulate", n_simulate = 1000)
appraise(gms$imp.gm, method = "simulate", n_simulate = 1000)
appraise(gms$wsa.gm, method = "simulate", n_simulate = 1000)
appraise(gms$alf.gm, method = "simulate", n_simulate = 1000)

# Permanovas
# Multivariate Model (PERMANOVA)
adonis2(otu_table(f)~f_tr$ctrl+f_tr$mp+f_tr$sdbs+f_tr$cu+f_tr$n+f_tr$nacl, permutations = 9999, 
        method = 'jac', binary=T, by="term", parallel=2) # only surfactant
adonis2(otu_table(f)~f_tr$m*f_tr$sdbs, permutations = 9999, 
        method = 'jac', binary=T, by="term", parallel=2) # test interaction

# Plots
# estimate uncertainty around fits
Myse <- function(mds, df, i){
  # function to extract the fit and standard error in the response scale
  # models are the list of models defined by Mymodsc
  # df is the list of data.frames that contain the predictor needed to pass to function predict
  # i is the number of iterations of this routine = number of datasets/models
  nd <- lapply(df, function(y) with(y, expand.grid(x = evenly(x, n=100), t=levels(t))))
  fit <- lapply(mds, function(x){do.call(tibble, predict(x, nd[[i]], type='response', se.fit=T, unconditional=T))})
  fit2 <- lapply(fit, function(x){cbind(x, nd[[i]])})
}
se <- Myse(mds=gms, df=rs, i=5)

# function for plots
Myplotsc=function(df, fml, laby, se){
  
  clrs <- viridis::viridis(3, option = "E", direction = -1)
  plot=list(ggplot(data=df, aes(x=x, y=y, colour = z))+
              annotate("rect", xmin=0.9,xmax=5.9,ymin=0,ymax=Inf, fill=alpha(clrs[3], 0.2), colour=NA)+
              annotate("rect", xmin=5.9,xmax=9.1,ymin=0,ymax=Inf, fill=alpha(clrs[2], 0.2), colour=NA)+
              annotate("rect", xmin=9.1,xmax=15.2,ymin=0,ymax=Inf, fill=alpha(clrs[1], 0.2), colour=NA)+
              geom_point(alpha=0.4)+
              geom_line(data=se, aes(x=x, y=fit, linetype = t), linewidth=0.5, inherit.aes = F)+ # use fit estimated by predict
              scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15))+
              #scale_color_manual(values=colrs)+
              scale_colour_viridis_c(option = "E", direction = -1)+
              labs(y=laby, x="Dryness level", colour="WHC %", linetype="ID of Pollutants")+
              theme_bw()+
              theme(text = element_text(size = 12),
                    panel.grid.minor = element_blank()))
}

laby <- list('Log(cop. g-1)', 'MWD (mm)', 'WDPT (sec.)',
             'WSA', '1/D (ASVs)')

pl <- mapply(Myplotsc, rs, fml, laby, se)

pl$cf <- pl$cf + coord_cartesian(ylim=c(23, 24.26)) # zooms in without losing the background of thresholds

# Plot ordinations
# beta diversity 
# points for NMDS
fs <- subset_samples(f, treatment%in%c('ctrl','sdbs'))
betf <- vegdist(otu_table(fs), method = "jac", binary = T)%>%
  metaMDS(., parallel=2, trace=F, weakties=T)

f_trs <- sample_data(fs)[,c(1,2,3,4)] %>% data.frame(.) %>% rename(t=treatment, m=moistureL, z=whc) %>%
        mutate_at('t', as.factor)%>%
        mutate(t=fct_recode(t, Ctrl='ctrl', Surfactant='sdbs'))%>%
        mutate(p=case_when(m < 6 ~ "Wet",
                     m < 10  ~ "Medium",
                     m >= 10 ~ "Dry" )) %>%
        mutate(p=fct_relevel(p, c("Wet", "Medium", "Dry")))

ford <- gg_ordiplot(betf, groups = f_trs$m, hull = F, label = F,
                    spiders = F, ellipse = F, plot = F, choices = c(1, 2), scaling=1)%>%
        pluck('df_ord')%>%
        rownames_to_column("potID")%>%
        mutate_at('potID', as.numeric)%>%
        left_join(.,f_trs, by="potID") %>% select(-Group)

# Spider for NMDS
sord <- gg_ordiplot(betf, groups = f_trs$t, hull = F, label = F,
                    spiders = T, ellipse = F, plot = F, choices = c(1, 2), scaling=1)%>%
        pluck('df_spiders')%>%
        rownames_to_column("potID")%>%
        mutate_at('potID', as.numeric)%>%
        left_join(.,f_tr, by="potID") %>% select(-Group)#%>%filter(p%notin%c('Medium'))

#ordination
clrs <- viridis::viridis(3, option = "E", direction = 1)
ord <- ggplot() + 
       geom_vline(xintercept=0.0, color="Grey", linewidth=1, linetype=1)+
       geom_hline(yintercept=0.0, color="Grey", linewidth=1, linetype=1)+    
       geom_point(data=ford, aes(x=x, y=y, colour=p, shape=t), alpha=0.8,  size=2, show.legend = T) +
       #geom_segment(data = sord, mapping = aes(x = cntr.x, y = cntr.x, xend= x, yend = y), alpha=0.2)+
       scale_color_manual(values=clrs[c(1,2,3)])+
       #scale_colour_viridis_c(option = "E", direction = -1)+
       facet_wrap(.~t)+
       labs(x="Axis 1", y="Axis 2", color="Phase", shape="")+
       theme_bw()+
       theme(axis.title = element_text(size = 12),
             legend.key = element_blank(),  #removes the box around each legend item
             legend.position = "bottom", #legend position
             legend.text = element_text(size=12),
             panel.border = element_rect(colour = "Black", fill = F),
             panel.grid = element_blank())    

(fig3 <- ggarrange(pl$mwd, pl$imp, pl$cf, pl$alf, ncol = 2, nrow = 2,
                 common.legend = T, legend = "right", labels = c('a', 'b', 'c', 'd')))
(figs2 <- ggarrange(pl$wsa, ord, ncol = 1, nrow = 2, 
                common.legend = F, legend = "right", labels = c('a', 'b')))  
# Figure 4, figs3 and tables -- Multiple factor effects  ####
'%notin%' <- Negate('%in%') # useful custom function to filter out
# load data and select all treatments except for control
rs <- read_delim("rs_db.csv", delim = ",") %>% 
      mutate(pha=case_when(dl < 6 ~ "Wet",
                          dl < 10  ~ "Medium",
                          dl >= 10 ~ "Dry" )) %>%
      mutate(pha=fct_relevel(pha, c("Wet", "Medium", "Dry")))%>%
      select(-c(8,10))%>%
      select(1:6,9,7,8)%>%
      mutate_at(c('no_poll', "treatment"), as.factor)%>%
      mutate(treatment=fct_relevel(treatment, c("Ctrl", "Copper","µPlastic", "Nitrogen", "Salinity", "Surfactant", "All_in")))%>%   
      pivot_longer(
        cols = 8:9,
      names_to = "var",
      values_to = "y")%>%
      select(potID, dl, whc, treatment, no_poll, pha, var, y)%>%
      rename(x=dl, z=whc, n=no_poll, t=treatment, p=pha)%>%
      filter(!is.na(y))%>%
      group_split(var)%>% # split into several databases by treatment and store on a list
      set_names(., nm = c("mwd", "imp"))%>% #give names to each set  
      as.list(.)

# community data
f <- readRDS("RS_all.rds")%>% subset_taxa(., Kingdom=="Fungi")# eliminate non target organisms
     
f <- subset_samples(f, sample_names(f) %notin% c("294","295", "285", "304", "313"))%>% # Eliminate samples with ultra low read numbers (<1000)
     filter_taxa(., function(x) sum(x) > 0, TRUE) # keep variants with more than 0 reads

# alpha diversity
alf <- cbind(sample_data(f)[,c(1)], estimate_richness(f, measures=c("InvSimpson")))%>%
       mutate(var=rep("invsimp",294))%>%
       right_join(.,rs[["mwd"]], by="potID")%>% #fixed 
       select(-c(y, var.y))%>%
       rename(var=var.x, y=InvSimpson)%>%
       select(1,4,5,6,7,8,3,2) %>%
       filter(!is.na(y)) %>% as_tibble(.) %>%
       mutate_at('potID', as.integer) 

rs$alf <- alf # append to db list
rm(alf)

# beta diversity
betf <- vegdist(otu_table(f), method = "jac", binary = T)%>%
  metaMDS(., parallel=2, trace=F, weakties=T)

f_tr <- sample_data(f)[,c(1,2,3,4,6)] %>% data.frame(.) %>% rename(t=treatment, m=moistureL, z=whc, n=no_poll) %>%
        mutate_at('t', as.factor)%>%
        mutate_at('n', as.factor)%>%
        mutate(p=case_when(m < 6 ~ "Wet",
                           m < 10  ~ "Medium",
                           m >= 10 ~ "Dry" )) %>%
        mutate(p=fct_relevel(p, c("Wet", "Medium", "Dry")))

# function for Models (GAMS)
Mygamsc=function(df, fml){
  # fits a GAM to each dataset on a list in df with max 3 base functions
  # returns a list of named gams.
  # argument fml is a string list stating the family of the conditional prob. and must be provided separately
  md=list(gm=gam(data = df, formula = y ~ n + s(x, by=n, k=3), family = fml, method="REML", select = T))
}

fml <- list('gaussian', 'poisson', 'gaussian')

gms <- mapply(Mygamsc, rs, fml) # iteratively apply Mygamsc

sm <- lapply(gms, summary)%>% # get each model summary
  lapply(., '[[', 's.table')%>% # extract each s.table
  imap_dfr(., ~ bind_cols(mod = .y,  edf = .x[,1], ref.df = .x[,2], 
                          Fval = round(.x[,3], 3), pval = round(.x[,4], 3))) # format stats

pm <- lapply(gms, summary)%>% # get each model summary
  lapply(., '[[', 'p.table')%>%
  imap_dfr(., ~ bind_cols(mod = .y, term = rownames(.x), bet = round(.x[,1], 3), StdE = round(.x[,2],3), 
                          tval = round(.x[,3], 3), pval = round(.x[,4], 3)))#%>% # format stats

# Tables
kable(sm, format = "simple", caption = "Table 6", digits = 3) # pretty simple table
kable(pm, format = "simple", caption = "Table 7", digits = 3) # pretty simple table
write_csv(sm, "table6") # good old way to save results
write_csv(pm, "table7") # good old way to save results
# check model fit with diagnostics 
#overview(mds$alf.mod)
check <- function(b, k.sample = 5000, k.rep = 200) {
  mgcv:::k.check(b, subsample = k.sample, n.rep = k.rep)
}

# function to check if k is correctly specified
chs <- lapply(gms, check)%>%
  imap_dfr(., ~ bind_cols(mod = .y,  k = .x[,1], edf = .x[,2], 
                          kindex = round(.x[,3], 3), pval = round(.x[,4], 3))) # format stats

appraise(gms$mwd.gm, method = "simulate", n_simulate = 1000)
appraise(gms$imp.gm, method = "simulate", n_simulate = 1000)
appraise(gms$alf.gm, method = "simulate", n_simulate = 1000)

# Permanovas
# Multivariate Model (PERMANOVA)
adonis2(otu_table(f)~f_tr$m*f_tr$n, permutations = 999, method = 'jac', binary=T, by="term")

# Plots
# estimate uncertainty around fits
Myse <- function(mds, df, i){
  # function to extract the fit and standard error in the response scale
  # models are the list of models defined by Mymodsc
  # df is the list of data.frames that contain the predictor needed to pass to function predict
  # i is the number of iterations of this routine = number of datasets/models
  nd <- lapply(df, function(y) with(y, expand.grid(x = evenly(x, n=100), n=levels(n))))
  fit <- lapply(mds, function(x){do.call(tibble, predict(x, nd[[i]], type='response', se.fit=T, unconditional=T))})
  fit2 <- lapply(fit, function(x){cbind(x, nd[[i]])})
}
se <- Myse(mds=gms, df=rs, i=3)

# function for plots
Myplotsc=function(df, fml, laby, se){
  
clrs <- viridis::viridis(3, option = "E", direction = -1)
plot=list(ggplot(data=df, aes(x=x, y=y, colour=z))+
              annotate("rect", xmin=0.9,xmax=5.9,ymin=0,ymax=Inf, fill=alpha(clrs[3], 0.2), colour=NA)+
              annotate("rect", xmin=5.9,xmax=9.1,ymin=0,ymax=Inf, fill=alpha(clrs[2], 0.2), colour=NA)+
              annotate("rect", xmin=9.1,xmax=15.2,ymin=0,ymax=Inf, fill=alpha(clrs[1], 0.2), colour=NA)+
              geom_point(alpha=0.4)+
              geom_line(data=se, aes(x=x, y=fit, linetype = n), linewidth=0.5, inherit.aes = F)+ # use fit estimated by predict
              scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15))+
              #scale_color_manual(values=colrs)+
              scale_colour_viridis_c(option = "E", direction = -1)+
              labs(y=laby, x="Dryness level", linetype="No of Pollutants", colour="WHC %")+
              theme_bw()+
              theme(text = element_text(size = 12),
                    panel.grid.minor = element_blank()))
}

laby <- list('Mean weight diameter (mm)', 'Water drop penetration time (sec.)','Inverse Simpson index (ASVs)')

pl <- mapply(Myplotsc, rs, fml, laby, se)

# Plot ordinations
# beta diversity 
ford <- gg_ordiplot(betf, groups = f_tr$n, hull = F, label = F,
                    spiders = F, ellipse = F, plot = F, choices = c(1, 2), scaling=1)%>%
        pluck('df_ord')%>%
        rownames_to_column("potID")%>%
        mutate_at('potID', as.numeric)%>%
        left_join(.,f_tr, by="potID") %>% select(-Group)%>%filter(n%notin%c('0', '1'))%>%filter(p%notin%c('Medium'))


sord <- gg_ordiplot(betf, groups = f_tr$p, hull = F, label = F,
                    spiders = T, ellipse = F, plot = F, choices = c(1, 2), scaling=1)%>%
        pluck('df_spiders')%>%
        rownames_to_column("potID")%>%
        mutate_at('potID', as.numeric)%>%
        left_join(.,f_tr, by="potID") %>% select(-Group)%>%filter(n%notin%c('0', '1'))%>%filter(p%notin%c('Medium'))
        
# ford <- ford %>% filter(potID%notin%c(285, 294, 295, 296, 304, 313)) # this are suspect of being contaminated
# sord <- sord %>% filter(potID%notin%c(285, 294, 295, 296, 304, 313))

#ordination
clrs <- viridis::viridis(3, option = "E", direction = 1)
ord <- ggplot() + 
       geom_vline(xintercept=0.0, color="Grey", linewidth=1, linetype=1)+
       geom_hline(yintercept=0.0, color="Grey", linewidth=1, linetype=1)+    
       geom_point(data=ford, aes(x=x, y=y, colour=p), alpha=0.8,  size=2, show.legend = T) +
       geom_segment(data = sord, mapping = aes(x = cntr.x, y = cntr.x, xend= x, yend = y), alpha=0.2)+
       scale_color_manual(values=clrs[c(1,3)])+
       #scale_colour_viridis_c(option = "E", direction = -1)+
       labs(x="Axis 1", y="Axis 2", color="Phase")+
       #facet_wrap(.~n, scales='free')+
       theme_bw()+
       theme(axis.title = element_text(size = 12),
             legend.key = element_blank(),  #removes the box around each legend item
              legend.position = "right", #legend position
              legend.text = element_text(size=12),
              panel.border = element_rect(colour = "Black", fill = F),
              panel.grid = element_blank())        

pl$ord <- ord # append ordination plot to list

(p1 <- ggarrange(pl$mwd, pl$imp, ncol = 2, nrow = 1,
                    common.legend = T, legend = "top", labels = c('a', 'b')))

(fig4 <- ggarrange(p1, ord, ncol = 1, nrow = 2, common.legend = F, legend = "bottom", labels = c('', 'c')))
#

