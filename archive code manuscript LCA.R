
# last edit: 05/03/2020

# accompanying code to the manuscript:
# "Heterogeneity in quality of life of long-term colon cancer survivors: 
# a latent class analysis of the population-based PROFILES registry"

# F.J. Clouth, MSc.
# Department of Methodology and Statistics, Tilburg University
# f.j.clouth@tilburguniversity.edu


# 1 data processing       line 45
# 2 run latentGOLD        line 258
# 3 data processing       line 324
# 4 run latentGOLD        line 380
# Table 1                 line 466
# Table 2                 line 497
# Table 3                 line 555
# Table 4                 line 608
# Figure 1                line 613
# Figure 2                line 677
# Figure 3                line 1535


library(dplyr)
library(alluvial)
library(data.table)
library(survival)
library(readstata13)
library(ggplot2)
library(gridExtra)
library(gtable)
library(foreign)
library(rio)
library(zoo)
library(car)
library(survminer)
library(tidyr)
library(grid)
library(Gmisc)
library(hrbrthemes)
library(RColorBrewer)
library(viridis)

setwd("adjust wd")

kern_crc <- read.csv("kern_crc.csv")

profiles.data <- read.dta13("profiles10000dataset_v6_17012018.dta")

crcw1 <- read.spss("crcw1final_pv28052014.sav", to.data.frame = T)

crcw2 <- read.spss("crcw2final_pv28052014.sav", to.data.frame = T)

crcw3 <- read.spss("crcw3final_pv28052014.sav", to.data.frame = T)

newmerge <- merge(kern_crc, profiles.data, by = "eid")

newmerge <- newmerge %>%
  filter(studie == 7) %>%
  filter(topo == "C18") %>%
  filter(topog != "C181") %>%
  filter(!morf %in% c("8002", "8013", "8041", "8042", "8043", "8044", "8045", "8150", 
                      "8151", "8152", "8153", "8154", "8155", "8156", "8157", "8243", 
                      "8244", "8245", "8246")) %>%
  filter(response == 1) %>%
  filter(stadium %in% c("1", "2", "2A", "2B", "2C", "3", "3A", "3B", "3C"))

newmerge <- newmerge %>%
  select(rn.x, eid, diffgr, topo, topog, morf, lyond, lypos, vitstat, vitfup)

newmerge$RN <- newmerge$rn.x

w1 <- crcw1 %>%
  filter(STUDYSTATUS == "responder") %>%
  filter(local == "colon") %>%
  filter(stage %in% c("Stage 1", "Stage 2", "Stage 3"))

w1 <- merge(w1, newmerge, by = "RN")

w1 <- w1 %>%
  filter(topog %in% c("C180", "C182", "C183", "C184", "C185", "C186", "C187", "C188", "C189")) %>%
  filter(!morf %in% c("8002", "8013", "8041", "8042", "8043", "8044", "8045", "8150", 
                      "8151", "8152", "8153", "8154", "8155", "8156", "8157", "8243", 
                      "8244", "8245", "8246"))

w1 <- w1 %>%
  select(RN, eid, diffgr, topog, morf, vitstat, vitfup, stage, stoma_kr, CT, surgery, RT,
         T_SES3, lyond, lypos, LEEFTIJD, TIMESINCEDIAG, GESLACHT, 
         P_HART, P_BER, P_HBLO, P_LONG, P_SUI, P_MAAG, P_NIER, P_LEVER, P_ANEM, P_SCHILD,
         P_DEPRES, P_ARTROS, P_RUG, P_REUMA, AANT_COMORB, OPLEID,
         QL, PF, RF, EF, CF, SF, FA, NV, PA, DY, SL, AP, CO, DI, FI)

w1 <- w1 %>%
  rename(chemo = CT,
         radio = RT,
         ses = T_SES3,
         edu = OPLEID,
         age = LEEFTIJD,
         t.s.diag = TIMESINCEDIAG,
         sex = GESLACHT,
         comorb = AANT_COMORB,
         QL.1 = QL, PF.1 = PF, RF.1 = RF, EF.1 = EF, CF.1 = CF,
         SF.1 = SF, FA.1 = FA, NV.1 = NV, PA.1 = PA, DY.1 = DY,
         SL.1 = SL, AP.1 = AP, CO.1 = CO, DI.1 = DI, FI.1 = FI)

w1$pf.1 <- car::recode(w1$PF.1, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w1$rf.1 <- car::recode(w1$RF.1, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w1$ef.1 <- car::recode(w1$EF.1, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w1$cf.1 <- car::recode(w1$CF.1, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w1$sf.1 <- car::recode(w1$SF.1, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w1$ql.1 <- car::recode(w1$QL.1, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w1$fa.1 <- car::recode(w1$FA.1, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w1$nv.1 <- car::recode(w1$NV.1, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w1$pa.1 <- car::recode(w1$PA.1, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w1$dy.1 <- car::recode(w1$DY.1, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w1$sl.1 <- car::recode(w1$SL.1, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w1$ap.1 <- car::recode(w1$AP.1, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w1$co.1 <- car::recode(w1$CO.1, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w1$di.1 <- car::recode(w1$DI.1, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w1$fi.1 <- car::recode(w1$FI.1, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")

w1$diffgr1 <- car::recode(w1$diffgr, "1 = 'Well differentiated'; 2 = 'Moderately differentiated';
                         3 = 'Poorly differentiated'; 4:9 = NA")

w1$diffgr2 <- car::recode(w1$diffgr, "1 = 'Well differentiated'; 2 = 'Moderately differentiated';
                         3 = 'Poorly differentiated'; 4:9 = 'unknown'")

w1$topog <- car::recode(w1$topog, "c('C180', 'C182', 'C183', 'C184') = 'proximal'; 
                        c('C185', 'C186', 'C187') = 'distal';
                        c('C188', 'C189') = 'other'")

w1$morf <- car::recode(w1$morf, "c(8470, 8480, 8481) = 'mucinous'; 
                       8490 = 'signet ring cell'; else = 'adenocarcinoma'")

w2 <- crcw2 %>%
  filter(STUDYSTATUS == "responder") %>%
  filter(local == "colon") %>%
  filter(stage %in% c("Stage 1", "Stage 2", "Stage 3"))

w2 <- merge(w2, newmerge, by = "RN")

w2 <- w2 %>%
  filter(topog != "C181") %>%
  filter(!morf %in% c("8002", "8013", "8041", "8042", "8043", "8044", "8045", "8150", 
                      "8151", "8152", "8153", "8154", "8155", "8156", "8157", "8243", 
                      "8244", "8245", "8246"))

w2 <- w2 %>%
  select(RN, QL, PF, RF, EF, CF, SF, FA, NV, PA, DY, SL, AP, CO, DI, FI)

w2 <- w2 %>%
  rename(QL.2 = QL, PF.2 = PF, RF.2 = RF, EF.2 = EF, CF.2 = CF,
         SF.2 = SF, FA.2 = FA, NV.2 = NV, PA.2 = PA, DY.2 = DY,
         SL.2 = SL, AP.2 = AP, CO.2 = CO, DI.2 = DI, FI.2 = FI)

w2$pf.2 <- car::recode(w2$PF.2, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w2$rf.2 <- car::recode(w2$RF.2, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w2$ef.2 <- car::recode(w2$EF.2, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w2$cf.2 <- car::recode(w2$CF.2, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w2$sf.2 <- car::recode(w2$SF.2, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w2$ql.2 <- car::recode(w2$QL.2, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w2$fa.2 <- car::recode(w2$FA.2, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w2$nv.2 <- car::recode(w2$NV.2, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w2$pa.2 <- car::recode(w2$PA.2, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w2$dy.2 <- car::recode(w2$DY.2, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w2$sl.2 <- car::recode(w2$SL.2, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w2$ap.2 <- car::recode(w2$AP.2, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w2$co.2 <- car::recode(w2$CO.2, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w2$di.2 <- car::recode(w2$DI.2, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w2$fi.2 <- car::recode(w2$FI.2, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")

w3 <- crcw3 %>%
  filter(STUDYSTATUS == "responder") %>%
  filter(local == "colon") %>%
  filter(stage %in% c("Stage 1", "Stage 2", "Stage 3"))

w3 <- merge(w3, newmerge, by = "RN")

w3 <- w3 %>%
  filter(topog != "C181") %>%
  filter(!morf %in% c("8002", "8013", "8041", "8042", "8043", "8044", "8045", "8150", 
                      "8151", "8152", "8153", "8154", "8155", "8156", "8157", "8243", 
                      "8244", "8245", "8246"))

w3 <- w3 %>%
  select(RN, QL, PF, RF, EF, CF, SF, FA, NV, DY, SL, AP, CO, DI, FI)

w3$PA <- NA

w3 <- w3 %>%
  rename(QL.3 = QL, PF.3 = PF, RF.3 = RF, EF.3 = EF, CF.3 = CF,
         SF.3 = SF, FA.3 = FA, NV.3 = NV, PA.3 = PA, DY.3 = DY,
         SL.3 = SL, AP.3 = AP, CO.3 = CO, DI.3 = DI, FI.3 = FI)

w3$pf.3 <- car::recode(w3$PF.3, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w3$rf.3 <- car::recode(w3$RF.3, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w3$ef.3 <- car::recode(w3$EF.3, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w3$cf.3 <- car::recode(w3$CF.3, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w3$sf.3 <- car::recode(w3$SF.3, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w3$ql.3 <- car::recode(w3$QL.3, "0:25 = 1; 26:50 = 2; 51:75 = 3; 76:100 = 4")
w3$fa.3 <- car::recode(w3$FA.3, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w3$nv.3 <- car::recode(w3$NV.3, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w3$pa.3 <- car::recode(w3$PA.3, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w3$dy.3 <- car::recode(w3$DY.3, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w3$sl.3 <- car::recode(w3$SL.3, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w3$ap.3 <- car::recode(w3$AP.3, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w3$co.3 <- car::recode(w3$CO.3, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w3$di.3 <- car::recode(w3$DI.3, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")
w3$fi.3 <- car::recode(w3$FI.3, "0:25 = 4; 26:50 = 3; 51:75 = 2; 76:100 = 1")


crcw.all <- merge(w1, w2, by = "RN", all.x = T)
crcw.all <- merge(crcw.all, w3, by = "RN", all.x = T)


crcw.long1 <- crcw.all %>%
  select(RN:edu, pf.1:fi.1, diffgr1, diffgr2)

crcw.long1$wave <- 1

crcw.long1 <- crcw.long1 %>%
  rename(pf = pf.1, rf = rf.1, ef = ef.1, cf = cf.1, sf = sf.1,
         ql = ql.1, fa = fa.1, nv = nv.1, pa = pa.1, dy = dy.1,
         sl = sl.1, ap = ap.1, co = co.1, di = di.1, fi = fi.1)

crcw.long2 <- crcw.all %>%
  select(RN:edu, pf.2:fi.2, diffgr1, diffgr2)

crcw.long2$wave <- 2

crcw.long2 <- crcw.long2 %>%
  rename(pf = pf.2, rf = rf.2, ef = ef.2, cf = cf.2, sf = sf.2,
         ql = ql.2, fa = fa.2, nv = nv.2, pa = pa.2, dy = dy.2,
         sl = sl.2, ap = ap.2, co = co.2, di = di.2, fi = fi.2)

crcw.long3 <- crcw.all %>%
  select(RN:edu, pf.3:fi.3, diffgr1, diffgr2)

crcw.long3$wave <- 3

crcw.long3 <- crcw.long3 %>%
  rename(pf = pf.3, rf = rf.3, ef = ef.3, cf = cf.3, sf = sf.3,
         ql = ql.3, fa = fa.3, nv = nv.3, pa = pa.3, dy = dy.3,
         sl = sl.3, ap = ap.3, co = co.3, di = di.3, fi = fi.3)

crcw.long <- rbind(crcw.long1, crcw.long2)
crcw.long <- rbind(crcw.long, crcw.long3)

crcw.long <- crcw.long %>%
  select(RN, wave, diffgr1, diffgr2, diffgr:edu, pf:fi)

export(crcw.long, "CRCW_all_long.sav")




# run latentGOLD
# Random Seed	185925						
# Best Start Seed	456866

makeNewSyntax <- function(syntaxName){
  
  newSyntaxToBe <- utils::capture.output(cat(paste("
//LG5.1//
version = 5.1
infile 'C:\\adjust wd\\CRCW_all_long.sav'
model
options
  maxthreads=8;
  algorithm 
    tolerance=1e-008 emtolerance=0,01 emiterations=250 nriterations=50;
  startvalues
    seed=0 sets=16 tolerance=1e-005 iterations=50;
  bayes
    categorical=1 variances=1 latent=1 poisson=1;
  montecarlo
    seed=0 sets=0 replicates=500 tolerance=1e-008;
  quadrature  nodes=10;
  missing  includeall;
  output      
    parameters=effect  betaopts=wl standarderrors profile probmeans=posterior
    loadings bivariateresiduals estimatedvalues=model reorderclasses;
  outfile  'C:\\adjust wd\\5class_all_long.sav'
    classification keep RN, diffgr, topog, morf, vitstat, vitfup, stage, stoma_kr,
       chemo, surgery, radio, ses, lyond, lypos, age, t.s.diag, sex, comorb;
variables
  holdout cases wave !=  1;
  dependent pf, rf, ef, cf, sf, ql, fa, nv, pa, dy, sl, ap, co, di, fi;
  latent
    Cluster nominal 5;
equations
  Cluster <- 1;
  pf <- 1 + Cluster;
  rf <- 1 + Cluster;
  ef <- 1 + Cluster;
  cf <- 1 + Cluster;
  sf <- 1 + Cluster;
  ql <- 1 + Cluster;
  fa <- 1 + Cluster;
  nv <- 1 + Cluster;
  pa <- 1 + Cluster;
  dy <- 1 + Cluster;
  sl <- 1 + Cluster;
  ap <- 1 + Cluster;
  co <- 1 + Cluster;
  di <- 1 + Cluster;
  fi <- 1 + Cluster;
end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}

LG <- "adjust wd/LatentGOLD5.1/lg51.exe"

makeNewSyntax(syntaxName = "5classSyntax")

setwd("adjust wd")

shell(paste(LG, "5classSyntax.lgs", "/b"))




class5_all_long <- read.spss("5class_all_long.sav", to.data.frame = T)

class5_all_long$MV <- 0

class5_all_long1 <- class5_all_long %>%
  filter(wave == 1)
for(i in 1:1489) {
  class5_all_long1$MV[i][all(is.na(c(class5_all_long1$pf[i], class5_all_long1$rf[i], class5_all_long1$ef[i], 
                                     class5_all_long1$cf[i], class5_all_long1$sf[i], class5_all_long1$ql[i], 
                                     class5_all_long1$fa[i], class5_all_long1$nv[i], class5_all_long1$pa[i], 
                                     class5_all_long1$dy[i], class5_all_long1$sl[i], class5_all_long1$ap[i], 
                                     class5_all_long1$co[i], class5_all_long1$di[i], class5_all_long1$fi[i])))] <- 1
}

class5_all_long2 <- class5_all_long %>%
  filter(wave == 2)
for(i in 1:1489) {
  class5_all_long2$MV[i][all(is.na(c(class5_all_long2$pf[i], class5_all_long2$rf[i], class5_all_long2$ef[i], 
                                  class5_all_long2$cf[i], class5_all_long2$sf[i], class5_all_long2$ql[i], 
                                  class5_all_long2$fa[i], class5_all_long2$nv[i], class5_all_long2$pa[i], 
                                  class5_all_long2$dy[i], class5_all_long2$sl[i], class5_all_long2$ap[i], 
                                  class5_all_long2$co[i], class5_all_long2$di[i], class5_all_long2$fi[i])))] <- 1
}

class5_all_long3 <- class5_all_long %>%
  filter(wave == 3)
for(i in 1:1489) {
  class5_all_long3$MV[i][all(is.na(c(class5_all_long3$pf[i], class5_all_long3$rf[i], class5_all_long3$ef[i], 
                                     class5_all_long3$cf[i], class5_all_long3$sf[i], class5_all_long3$ql[i], 
                                     class5_all_long3$fa[i], class5_all_long3$nv[i], class5_all_long3$pa[i], 
                                     class5_all_long3$dy[i], class5_all_long3$sl[i], class5_all_long3$ap[i], 
                                     class5_all_long3$co[i], class5_all_long3$di[i], class5_all_long3$fi[i])))] <- 1
}


class5_all_long <- as.data.frame(rbind(class5_all_long1, rbind(class5_all_long2, class5_all_long3)))

class5_all_long$Cluster.1[class5_all_long$MV == 1] <- NA
class5_all_long$Cluster.2[class5_all_long$MV == 1] <- NA
class5_all_long$Cluster.3[class5_all_long$MV == 1] <- NA
class5_all_long$Cluster.4[class5_all_long$MV == 1] <- NA
class5_all_long$Cluster.5[class5_all_long$MV == 1] <- NA

class5_all_long2 <- class5_all_long %>%
  filter(morf != "signet ring cell") %>%
  filter(ses != 4) 

class5_all_long2$diffgr <- car::recode(class5_all_long2$diffgr, "9 = NA")

export(class5_all_long2, "5class_all_long2.sav")

export(class5_all_long, "5class_all_long.sav")



# Random Seed	486060						
# Best Start Seed	486060
makeNewSyntax2 <- function(syntaxName){
  
  newSyntaxToBe <- utils::capture.output(cat(paste("
//LG5.1//
version = 5.1
infile 'C:\\adjust wd\\5class_all_long2.sav'
model
  options
    maxthreads=8;
    algorithm 
      tolerance=1e-008 emtolerance=0,01 emiterations=250 nriterations=50 ;
    startvalues
      seed=0 sets=16 tolerance=1e-005 iterations=50;
    bayes
      categorical=1 variances=1 latent=1 poisson=1;
    montecarlo
      seed=202022 sets=0 replicates=500 tolerance=1e-008;
    quadrature  nodes=10;
    missing  includeall;
    step3 modal ml;
    output      
      parameters=first  betaopts=swl standarderrors=robust profile=posterior
      probmeans=posterior estimatedvalues=model reorderclasses;
  variables
    select wave =  1;
    independent sex nominal, age, ses nominal, comorb nominal, t.s.diag,
      stage nominal, diffgr nominal, topog nominal, morf nominal, chemo nominal;
    latent Cluster nominal posterior = ( Cluster.1 Cluster.2 Cluster.3 Cluster.4 Cluster.5 ) ;
  equations
    Cluster <- 1 + sex + age + ses + comorb + t.s.diag + stage + diffgr + topog + morf + chemo;
end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


makeNewSyntax2(syntaxName = "step3MNL")

shell(paste(LG, "step3MNL.lgs", "/b"))



#Random Seed	155347						
#Best Start Seed	155347

makeNewSyntax3 <- function(syntaxName){
  
  newSyntaxToBe <- utils::capture.output(cat(paste("
//LG5.1//
version = 5.1
infile 'C:\\adjust wd\\5class_all_long.sav'
model
  options
    maxthreads=8;
    algorithm 
      tolerance=1e-008 emtolerance=0,01 emiterations=250 nriterations=50;
    startvalues
      seed=0 sets=16 tolerance=1e-005 iterations=50;
    bayes
      categorical=1 variances=1 latent=0 poisson=1;
    montecarlo
      seed=0 sets=0 replicates=500 tolerance=1e-008;
    quadrature  nodes=10;
    missing  includeall;
    step3 modal ml;
    output      
      parameters=effect  betaopts=wl standarderrors=robust profile=posterior
      probmeans=posterior estimatedvalues=model reorderclasses;
    variables
      caseid RN;
      latent Cluster dynamic nominal posterior = ( Cluster.1 Cluster.2 Cluster.3 Cluster.4 Cluster.5 ) ;
    equations
      Cluster[=0] <- 1;
      Cluster <- (~tra) 1 | Cluster[-1];
end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


makeNewSyntax3(syntaxName = "step3LTA")

shell(paste(LG, "step3LTA.lgs", "/b"))



# Table 1.
class5_all_long_w1 <- class5_all_long %>%
  filter(wave == 1)


mean(class5_all_long_w1$age, na.rm = T)
sd(class5_all_long_w1$age, na.rm = T)
class5_all_long_w1 %>%
  filter(is.na(age)) %>%
  count()
table(class5_all_long_w1$sex, useNA = "always")
prop.table(table(class5_all_long_w1$sex, useNA = "always"))
table(class5_all_long_w1$ses, useNA = "always")
prop.table(table(class5_all_long_w1$ses, useNA = "always"))
mean(class5_all_long_w1$t.s.diag, na.rm = T)
sd(class5_all_long_w1$t.s.diag, na.rm = T)
class5_all_long_w1 %>%
  filter(is.na(t.s.diag)) %>%
  count()
table(class5_all_long_w1$stage, useNA = "always")
prop.table(table(class5_all_long_w1$stage, useNA = "always"))
table(class5_all_long_w1$diffgr, useNA = "always")
prop.table(table(class5_all_long_w1$diffgr, useNA = "always"))
table(class5_all_long_w1$topog, useNA = "always")
prop.table(table(class5_all_long_w1$topog, useNA = "always"))
table(class5_all_long_w1$morf, useNA = "always")
prop.table(table(class5_all_long_w1$morf, useNA = "always"))
table(class5_all_long_w1$chemo, useNA = "always")
prop.table(table(class5_all_long_w1$chemo, useNA = "always"))


# Table 2.
for(i in 1:7) {
makeNewSyntax <- function(syntaxName){
  
  newSyntaxToBe <- utils::capture.output(cat(paste("
//LG5.1//
version = 5.1
infile 'C:\\adjust wd\\CRCW_all_long.sav'
model
options
maxthreads=8;
algorithm 
tolerance=1e-008 emtolerance=0,01 emiterations=250 nriterations=50;
startvalues
seed=0 sets=16 tolerance=1e-005 iterations=50;
bayes
categorical=1 variances=1 latent=1 poisson=1;
montecarlo
seed=0 sets=0 replicates=500 tolerance=1e-008;
quadrature  nodes=10;
missing  includeall;
output      
parameters=effect  betaopts=wl standarderrors profile probmeans=posterior
loadings bivariateresiduals estimatedvalues=model reorderclasses;
variables
holdout cases wave !=  1;
dependent pf, rf, ef, cf, sf, ql, fa, nv, pa, dy, sl, ap, co, di, fi;
latent
Cluster nominal ", i,";
equations
Cluster <- 1;
pf <- 1 + Cluster;
rf <- 1 + Cluster;
ef <- 1 + Cluster;
cf <- 1 + Cluster;
sf <- 1 + Cluster;
ql <- 1 + Cluster;
fa <- 1 + Cluster;
nv <- 1 + Cluster;
pa <- 1 + Cluster;
dy <- 1 + Cluster;
sl <- 1 + Cluster;
ap <- 1 + Cluster;
co <- 1 + Cluster;
di <- 1 + Cluster;
fi <- 1 + Cluster;
end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}

makeNewSyntax(syntaxName = paste(i, "class", sep = ""))

shell(paste(LG, paste(i, "class.lgs", sep = ""), "/b"))
}




# Table 3.

eortc <- merge(crcw1, class5_all_long_w1, by="RN", all.y = T)


eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(QL, na.rm = T), sd(QL, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(PF, na.rm = T), sd(PF, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(RF, na.rm = T), sd(RF, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(EF, na.rm = T), sd(EF, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(CF, na.rm = T), sd(CF, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(SF, na.rm = T), sd(SF, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(FA, na.rm = T), sd(FA, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(NV, na.rm = T), sd(NV, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(PA, na.rm = T), sd(PA, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(DY, na.rm = T), sd(DY, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(SL, na.rm = T), sd(SL, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(AP, na.rm = T), sd(AP, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(CO, na.rm = T), sd(CO, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(DI, na.rm = T), sd(DI, na.rm = T))
eortc %>%
  group_by(Cluster.) %>%
  summarise(mean(FI, na.rm = T), sd(FI, na.rm = T))



# Table 4.
# see step3MNL.lst



# Figure 1.

png("rflow.png", width = 2900, height = 2100, res = 300)
grid.newpage()
# set some parameters to use repeatedly
leftx <- .25
midx <- .3
rightx <- .7
width <- .28
height <- .12
gp <- gpar(fill = "white")
# create boxes
(total <- boxGrob("Total\n N = 3875", 
                  x=midx, y=.9, box_gp = gp, width = width, height = height))

(resw1 <- boxGrob("Wave 1\n N = 1489", 
                  x=midx, y=.5, box_gp = gp, width = width, height = height))

(resw2 <- boxGrob("Wave 2\n N = 912", 
                  x=midx, y=.3, box_gp = gp, width = width, height = height))

(resw3 <- boxGrob("Wave 3\n N = 809", 
                  x=midx, y=.1, box_gp = gp, width = width, height = height))

(exclu <- boxGrob("Exclusions\n\nNon-respondents\nRectal\nMetastatic\nAppendix and NET", 
                  x=rightx, y=.7, box_gp = gp, width = width, height = .25, just = "left"))

(excluno <- boxGrob("2386\n\n1250\n1020\n97\n19", 
                    x=.9, y=.7, box_gp = gpar(alpha=0), width = width, height = .25, just = "left"))

(nores2 <- boxGrob("Non-respondents\n N = 577", 
                   x=rightx, y=.4, box_gp = gp, width = width, height = height))

(nores3 <- boxGrob("Non-respondents\n N = 103", 
                   x=rightx, y=.2, box_gp = gp, width = width, height = height))


connectGrob(total, resw1, "v", arrow_obj = 
              getOption("connectGrobArrow", default = arrow(angle = 25, length = unit(0.2, "inches"),
                                                            ends = "last", type = "closed")))

connectGrob(resw1, resw2, "v", arrow_obj = 
              getOption("connectGrobArrow", default = arrow(angle = 25, length = unit(0.2, "inches"),
                                                            ends = "last", type = "closed")))

connectGrob(resw2, resw3, "v", arrow_obj = 
              getOption("connectGrobArrow", default = arrow(angle = 25, length = unit(0.2, "inches"),
                                                            ends = "last", type = "closed")))

connectGrob(total, exclu, "-", arrow_obj = 
              getOption("connectGrobArrow", default = arrow(angle = 25, length = unit(0.2, "inches"),
                                                            ends = "last", type = "closed")))

connectGrob(resw1, nores2, "-", arrow_obj = 
              getOption("connectGrobArrow", default = arrow(angle = 25, length = unit(0.2, "inches"),
                                                            ends = "last", type = "closed")))

connectGrob(resw2, nores3, "-", arrow_obj = 
              getOption("connectGrobArrow", default = arrow(angle = 25, length = unit(0.2, "inches"),
                                                            ends = "last", type = "closed")))
dev.off()



# Figure 2.
w1 <- eortc %>%
  select(RN, Cluster., QL, PF, RF, EF, CF,
         SF, FA, NV, PA, DY,
         SL, AP, CO, DI, FI)

w1 <- w1 %>%
  rename(clu. = Cluster.,
         ql = QL, pf = PF, rf = RF, ef = EF, cf = CF,
         sf = SF, fa = FA, nv = NV, pa = PA, dy = DY,
         sl = SL, ap = AP, co = CO, di = DI, fi = FI)

table(w1$clu.)

c1_mean <- w1 %>%
  filter(clu. == 1) %>%
  summarise(., mean(ql, na.rm = T), mean(ql, na.rm = T) - 1.96*sd(ql, na.rm = T)/sqrt(587), mean(ql, na.rm = T) + 1.96*sd(ql, na.rm = T)/sqrt(587),
            mean(pf, na.rm = T), mean(pf, na.rm = T) - 1.96*sd(pf, na.rm = T)/sqrt(587), mean(pf, na.rm = T) + 1.96*sd(pf, na.rm = T)/sqrt(587),
            mean(rf, na.rm = T), mean(rf, na.rm = T) - 1.96*sd(rf, na.rm = T)/sqrt(587), mean(rf, na.rm = T) + 1.96*sd(rf, na.rm = T)/sqrt(587),
            mean(ef, na.rm = T), mean(ef, na.rm = T) - 1.96*sd(ef, na.rm = T)/sqrt(587), mean(ef, na.rm = T) + 1.96*sd(ef, na.rm = T)/sqrt(587),
            mean(cf, na.rm = T), mean(cf, na.rm = T) - 1.96*sd(cf, na.rm = T)/sqrt(587), mean(cf, na.rm = T) + 1.96*sd(cf, na.rm = T)/sqrt(587),
            mean(sf, na.rm = T), mean(sf, na.rm = T) - 1.96*sd(sf, na.rm = T)/sqrt(587), mean(sf, na.rm = T) + 1.96*sd(sf, na.rm = T)/sqrt(587),
            mean(fa, na.rm = T), mean(fa, na.rm = T) - 1.96*sd(fa, na.rm = T)/sqrt(587), mean(fa, na.rm = T) + 1.96*sd(fa, na.rm = T)/sqrt(587),
            mean(nv, na.rm = T), mean(nv, na.rm = T) - 1.96*sd(nv, na.rm = T)/sqrt(587), mean(nv, na.rm = T) + 1.96*sd(nv, na.rm = T)/sqrt(587),
            mean(pa, na.rm = T), mean(pa, na.rm = T) - 1.96*sd(pa, na.rm = T)/sqrt(587), mean(pa, na.rm = T) + 1.96*sd(pa, na.rm = T)/sqrt(587),
            mean(dy, na.rm = T), mean(dy, na.rm = T) - 1.96*sd(dy, na.rm = T)/sqrt(587), mean(dy, na.rm = T) + 1.96*sd(dy, na.rm = T)/sqrt(587),
            mean(sl, na.rm = T), mean(sl, na.rm = T) - 1.96*sd(sl, na.rm = T)/sqrt(587), mean(sl, na.rm = T) + 1.96*sd(sl, na.rm = T)/sqrt(587),
            mean(ap, na.rm = T), mean(ap, na.rm = T) - 1.96*sd(ap, na.rm = T)/sqrt(587), mean(ap, na.rm = T) + 1.96*sd(ap, na.rm = T)/sqrt(587),
            mean(co, na.rm = T), mean(co, na.rm = T) - 1.96*sd(co, na.rm = T)/sqrt(587), mean(co, na.rm = T) + 1.96*sd(co, na.rm = T)/sqrt(587),
            mean(di, na.rm = T), mean(di, na.rm = T) - 1.96*sd(di, na.rm = T)/sqrt(587), mean(di, na.rm = T) + 1.96*sd(di, na.rm = T)/sqrt(587),
            mean(fi, na.rm = T), mean(fi, na.rm = T) - 1.96*sd(fi, na.rm = T)/sqrt(587), mean(fi, na.rm = T) + 1.96*sd(fi, na.rm = T)/sqrt(587))

c2_mean <- w1 %>%
  filter(clu. == 2) %>%
  summarise(., mean(ql, na.rm = T), mean(ql, na.rm = T) - 1.96*sd(ql, na.rm = T)/sqrt(454), mean(ql, na.rm = T) + 1.96*sd(ql, na.rm = T)/sqrt(454),
            mean(pf, na.rm = T), mean(pf, na.rm = T) - 1.96*sd(pf, na.rm = T)/sqrt(454), mean(pf, na.rm = T) + 1.96*sd(pf, na.rm = T)/sqrt(454),
            mean(rf, na.rm = T), mean(rf, na.rm = T) - 1.96*sd(rf, na.rm = T)/sqrt(454), mean(rf, na.rm = T) + 1.96*sd(rf, na.rm = T)/sqrt(454),
            mean(ef, na.rm = T), mean(ef, na.rm = T) - 1.96*sd(ef, na.rm = T)/sqrt(454), mean(ef, na.rm = T) + 1.96*sd(ef, na.rm = T)/sqrt(454),
            mean(cf, na.rm = T), mean(cf, na.rm = T) - 1.96*sd(cf, na.rm = T)/sqrt(454), mean(cf, na.rm = T) + 1.96*sd(cf, na.rm = T)/sqrt(454),
            mean(sf, na.rm = T), mean(sf, na.rm = T) - 1.96*sd(sf, na.rm = T)/sqrt(454), mean(sf, na.rm = T) + 1.96*sd(sf, na.rm = T)/sqrt(454),
            mean(fa, na.rm = T), mean(fa, na.rm = T) - 1.96*sd(fa, na.rm = T)/sqrt(454), mean(fa, na.rm = T) + 1.96*sd(fa, na.rm = T)/sqrt(454),
            mean(nv, na.rm = T), mean(nv, na.rm = T) - 1.96*sd(nv, na.rm = T)/sqrt(454), mean(nv, na.rm = T) + 1.96*sd(nv, na.rm = T)/sqrt(454),
            mean(pa, na.rm = T), mean(pa, na.rm = T) - 1.96*sd(pa, na.rm = T)/sqrt(454), mean(pa, na.rm = T) + 1.96*sd(pa, na.rm = T)/sqrt(454),
            mean(dy, na.rm = T), mean(dy, na.rm = T) - 1.96*sd(dy, na.rm = T)/sqrt(454), mean(dy, na.rm = T) + 1.96*sd(dy, na.rm = T)/sqrt(454),
            mean(sl, na.rm = T), mean(sl, na.rm = T) - 1.96*sd(sl, na.rm = T)/sqrt(454), mean(sl, na.rm = T) + 1.96*sd(sl, na.rm = T)/sqrt(454),
            mean(ap, na.rm = T), mean(ap, na.rm = T) - 1.96*sd(ap, na.rm = T)/sqrt(454), mean(ap, na.rm = T) + 1.96*sd(ap, na.rm = T)/sqrt(454),
            mean(co, na.rm = T), mean(co, na.rm = T) - 1.96*sd(co, na.rm = T)/sqrt(454), mean(co, na.rm = T) + 1.96*sd(co, na.rm = T)/sqrt(454),
            mean(di, na.rm = T), mean(di, na.rm = T) - 1.96*sd(di, na.rm = T)/sqrt(454), mean(di, na.rm = T) + 1.96*sd(di, na.rm = T)/sqrt(454),
            mean(fi, na.rm = T), mean(fi, na.rm = T) - 1.96*sd(fi, na.rm = T)/sqrt(454), mean(fi, na.rm = T) + 1.96*sd(fi, na.rm = T)/sqrt(454))

c3_mean <- w1 %>%
  filter(clu. == 3) %>%
  summarise(., mean(ql, na.rm = T), mean(ql, na.rm = T) - 1.96*sd(ql, na.rm = T)/sqrt(208), mean(ql, na.rm = T) + 1.96*sd(ql, na.rm = T)/sqrt(208),
            mean(pf, na.rm = T), mean(pf, na.rm = T) - 1.96*sd(pf, na.rm = T)/sqrt(208), mean(pf, na.rm = T) + 1.96*sd(pf, na.rm = T)/sqrt(208),
            mean(rf, na.rm = T), mean(rf, na.rm = T) - 1.96*sd(rf, na.rm = T)/sqrt(208), mean(rf, na.rm = T) + 1.96*sd(rf, na.rm = T)/sqrt(208),
            mean(ef, na.rm = T), mean(ef, na.rm = T) - 1.96*sd(ef, na.rm = T)/sqrt(208), mean(ef, na.rm = T) + 1.96*sd(ef, na.rm = T)/sqrt(208),
            mean(cf, na.rm = T), mean(cf, na.rm = T) - 1.96*sd(cf, na.rm = T)/sqrt(208), mean(cf, na.rm = T) + 1.96*sd(cf, na.rm = T)/sqrt(208),
            mean(sf, na.rm = T), mean(sf, na.rm = T) - 1.96*sd(sf, na.rm = T)/sqrt(208), mean(sf, na.rm = T) + 1.96*sd(sf, na.rm = T)/sqrt(208),
            mean(fa, na.rm = T), mean(fa, na.rm = T) - 1.96*sd(fa, na.rm = T)/sqrt(208), mean(fa, na.rm = T) + 1.96*sd(fa, na.rm = T)/sqrt(208),
            mean(nv, na.rm = T), mean(nv, na.rm = T) - 1.96*sd(nv, na.rm = T)/sqrt(208), mean(nv, na.rm = T) + 1.96*sd(nv, na.rm = T)/sqrt(208),
            mean(pa, na.rm = T), mean(pa, na.rm = T) - 1.96*sd(pa, na.rm = T)/sqrt(208), mean(pa, na.rm = T) + 1.96*sd(pa, na.rm = T)/sqrt(208),
            mean(dy, na.rm = T), mean(dy, na.rm = T) - 1.96*sd(dy, na.rm = T)/sqrt(208), mean(dy, na.rm = T) + 1.96*sd(dy, na.rm = T)/sqrt(208),
            mean(sl, na.rm = T), mean(sl, na.rm = T) - 1.96*sd(sl, na.rm = T)/sqrt(208), mean(sl, na.rm = T) + 1.96*sd(sl, na.rm = T)/sqrt(208),
            mean(ap, na.rm = T), mean(ap, na.rm = T) - 1.96*sd(ap, na.rm = T)/sqrt(208), mean(ap, na.rm = T) + 1.96*sd(ap, na.rm = T)/sqrt(208),
            mean(co, na.rm = T), mean(co, na.rm = T) - 1.96*sd(co, na.rm = T)/sqrt(208), mean(co, na.rm = T) + 1.96*sd(co, na.rm = T)/sqrt(208),
            mean(di, na.rm = T), mean(di, na.rm = T) - 1.96*sd(di, na.rm = T)/sqrt(208), mean(di, na.rm = T) + 1.96*sd(di, na.rm = T)/sqrt(208),
            mean(fi, na.rm = T), mean(fi, na.rm = T) - 1.96*sd(fi, na.rm = T)/sqrt(208), mean(fi, na.rm = T) + 1.96*sd(fi, na.rm = T)/sqrt(208))

c4_mean <- w1 %>%
  filter(clu. == 4) %>%
  summarise(., mean(ql, na.rm = T), mean(ql, na.rm = T) - 1.96*sd(ql, na.rm = T)/sqrt(119), mean(ql, na.rm = T) + 1.96*sd(ql, na.rm = T)/sqrt(119),
            mean(pf, na.rm = T), mean(pf, na.rm = T) - 1.96*sd(pf, na.rm = T)/sqrt(119), mean(pf, na.rm = T) + 1.96*sd(pf, na.rm = T)/sqrt(119),
            mean(rf, na.rm = T), mean(rf, na.rm = T) - 1.96*sd(rf, na.rm = T)/sqrt(119), mean(rf, na.rm = T) + 1.96*sd(rf, na.rm = T)/sqrt(119),
            mean(ef, na.rm = T), mean(ef, na.rm = T) - 1.96*sd(ef, na.rm = T)/sqrt(119), mean(ef, na.rm = T) + 1.96*sd(ef, na.rm = T)/sqrt(119),
            mean(cf, na.rm = T), mean(cf, na.rm = T) - 1.96*sd(cf, na.rm = T)/sqrt(119), mean(cf, na.rm = T) + 1.96*sd(cf, na.rm = T)/sqrt(119),
            mean(sf, na.rm = T), mean(sf, na.rm = T) - 1.96*sd(sf, na.rm = T)/sqrt(119), mean(sf, na.rm = T) + 1.96*sd(sf, na.rm = T)/sqrt(119),
            mean(fa, na.rm = T), mean(fa, na.rm = T) - 1.96*sd(fa, na.rm = T)/sqrt(119), mean(fa, na.rm = T) + 1.96*sd(fa, na.rm = T)/sqrt(119),
            mean(nv, na.rm = T), mean(nv, na.rm = T) - 1.96*sd(nv, na.rm = T)/sqrt(119), mean(nv, na.rm = T) + 1.96*sd(nv, na.rm = T)/sqrt(119),
            mean(pa, na.rm = T), mean(pa, na.rm = T) - 1.96*sd(pa, na.rm = T)/sqrt(119), mean(pa, na.rm = T) + 1.96*sd(pa, na.rm = T)/sqrt(119),
            mean(dy, na.rm = T), mean(dy, na.rm = T) - 1.96*sd(dy, na.rm = T)/sqrt(119), mean(dy, na.rm = T) + 1.96*sd(dy, na.rm = T)/sqrt(119),
            mean(sl, na.rm = T), mean(sl, na.rm = T) - 1.96*sd(sl, na.rm = T)/sqrt(119), mean(sl, na.rm = T) + 1.96*sd(sl, na.rm = T)/sqrt(119),
            mean(ap, na.rm = T), mean(ap, na.rm = T) - 1.96*sd(ap, na.rm = T)/sqrt(119), mean(ap, na.rm = T) + 1.96*sd(ap, na.rm = T)/sqrt(119),
            mean(co, na.rm = T), mean(co, na.rm = T) - 1.96*sd(co, na.rm = T)/sqrt(119), mean(co, na.rm = T) + 1.96*sd(co, na.rm = T)/sqrt(119),
            mean(di, na.rm = T), mean(di, na.rm = T) - 1.96*sd(di, na.rm = T)/sqrt(119), mean(di, na.rm = T) + 1.96*sd(di, na.rm = T)/sqrt(119),
            mean(fi, na.rm = T), mean(fi, na.rm = T) - 1.96*sd(fi, na.rm = T)/sqrt(119), mean(fi, na.rm = T) + 1.96*sd(fi, na.rm = T)/sqrt(119))

c5_mean <- w1 %>%
  filter(clu. == 5) %>%
  summarise(., mean(ql, na.rm = T), mean(ql, na.rm = T) - 1.96*sd(ql, na.rm = T)/sqrt(121), mean(ql, na.rm = T) + 1.96*sd(ql, na.rm = T)/sqrt(121),
            mean(pf, na.rm = T), mean(pf, na.rm = T) - 1.96*sd(pf, na.rm = T)/sqrt(121), mean(pf, na.rm = T) + 1.96*sd(pf, na.rm = T)/sqrt(121),
            mean(rf, na.rm = T), mean(rf, na.rm = T) - 1.96*sd(rf, na.rm = T)/sqrt(121), mean(rf, na.rm = T) + 1.96*sd(rf, na.rm = T)/sqrt(121),
            mean(ef, na.rm = T), mean(ef, na.rm = T) - 1.96*sd(ef, na.rm = T)/sqrt(121), mean(ef, na.rm = T) + 1.96*sd(ef, na.rm = T)/sqrt(121),
            mean(cf, na.rm = T), mean(cf, na.rm = T) - 1.96*sd(cf, na.rm = T)/sqrt(121), mean(cf, na.rm = T) + 1.96*sd(cf, na.rm = T)/sqrt(121),
            mean(sf, na.rm = T), mean(sf, na.rm = T) - 1.96*sd(sf, na.rm = T)/sqrt(121), mean(sf, na.rm = T) + 1.96*sd(sf, na.rm = T)/sqrt(121),
            mean(fa, na.rm = T), mean(fa, na.rm = T) - 1.96*sd(fa, na.rm = T)/sqrt(121), mean(fa, na.rm = T) + 1.96*sd(fa, na.rm = T)/sqrt(121),
            mean(nv, na.rm = T), mean(nv, na.rm = T) - 1.96*sd(nv, na.rm = T)/sqrt(121), mean(nv, na.rm = T) + 1.96*sd(nv, na.rm = T)/sqrt(121),
            mean(pa, na.rm = T), mean(pa, na.rm = T) - 1.96*sd(pa, na.rm = T)/sqrt(121), mean(pa, na.rm = T) + 1.96*sd(pa, na.rm = T)/sqrt(121),
            mean(dy, na.rm = T), mean(dy, na.rm = T) - 1.96*sd(dy, na.rm = T)/sqrt(121), mean(dy, na.rm = T) + 1.96*sd(dy, na.rm = T)/sqrt(121),
            mean(sl, na.rm = T), mean(sl, na.rm = T) - 1.96*sd(sl, na.rm = T)/sqrt(121), mean(sl, na.rm = T) + 1.96*sd(sl, na.rm = T)/sqrt(121),
            mean(ap, na.rm = T), mean(ap, na.rm = T) - 1.96*sd(ap, na.rm = T)/sqrt(121), mean(ap, na.rm = T) + 1.96*sd(ap, na.rm = T)/sqrt(121),
            mean(co, na.rm = T), mean(co, na.rm = T) - 1.96*sd(co, na.rm = T)/sqrt(121), mean(co, na.rm = T) + 1.96*sd(co, na.rm = T)/sqrt(121),
            mean(di, na.rm = T), mean(di, na.rm = T) - 1.96*sd(di, na.rm = T)/sqrt(121), mean(di, na.rm = T) + 1.96*sd(di, na.rm = T)/sqrt(121),
            mean(fi, na.rm = T), mean(fi, na.rm = T) - 1.96*sd(fi, na.rm = T)/sqrt(121), mean(fi, na.rm = T) + 1.96*sd(fi, na.rm = T)/sqrt(121))

c1.spider <- as.data.frame(c1_mean)
c2.spider <- as.data.frame(c2_mean)
c3.spider <- as.data.frame(c3_mean)
c4.spider <- as.data.frame(c4_mean)
c5.spider <- as.data.frame(c5_mean)


c1.spider <- c1.spider %>%
  rename(ql = `mean(ql, na.rm = T)`,
         ql_l = `mean(ql, na.rm = T) - 1.96 * sd(ql, na.rm = T)/sqrt(587)`,
         ql_h = `mean(ql, na.rm = T) + 1.96 * sd(ql, na.rm = T)/sqrt(587)`,
         ef = `mean(ef, na.rm = T)`,
         ef_l = `mean(ef, na.rm = T) - 1.96 * sd(ef, na.rm = T)/sqrt(587)`,
         ef_h = `mean(ef, na.rm = T) + 1.96 * sd(ef, na.rm = T)/sqrt(587)`,
         cf = `mean(cf, na.rm = T)`,
         cf_l = `mean(cf, na.rm = T) - 1.96 * sd(cf, na.rm = T)/sqrt(587)`,
         cf_h = `mean(cf, na.rm = T) + 1.96 * sd(cf, na.rm = T)/sqrt(587)`,
         sf = `mean(sf, na.rm = T)`,
         sf_l = `mean(sf, na.rm = T) - 1.96 * sd(sf, na.rm = T)/sqrt(587)`,
         sf_h = `mean(sf, na.rm = T) + 1.96 * sd(sf, na.rm = T)/sqrt(587)`,
         sl = `mean(sl, na.rm = T)`,
         sl_l = `mean(sl, na.rm = T) - 1.96 * sd(sl, na.rm = T)/sqrt(587)`,
         sl_h = `mean(sl, na.rm = T) + 1.96 * sd(sl, na.rm = T)/sqrt(587)`,
         nv = `mean(nv, na.rm = T)`,
         nv_l = `mean(nv, na.rm = T) - 1.96 * sd(nv, na.rm = T)/sqrt(587)`,
         nv_h = `mean(nv, na.rm = T) + 1.96 * sd(nv, na.rm = T)/sqrt(587)`,
         ap = `mean(ap, na.rm = T)`,
         ap_l = `mean(ap, na.rm = T) - 1.96 * sd(ap, na.rm = T)/sqrt(587)`,
         ap_h = `mean(ap, na.rm = T) + 1.96 * sd(ap, na.rm = T)/sqrt(587)`,
         co = `mean(co, na.rm = T)`,
         co_l = `mean(co, na.rm = T) - 1.96 * sd(co, na.rm = T)/sqrt(587)`,
         co_h = `mean(co, na.rm = T) + 1.96 * sd(co, na.rm = T)/sqrt(587)`,
         di = `mean(di, na.rm = T)`,
         di_l = `mean(di, na.rm = T) - 1.96 * sd(di, na.rm = T)/sqrt(587)`,
         di_h = `mean(di, na.rm = T) + 1.96 * sd(di, na.rm = T)/sqrt(587)`,
         fi = `mean(fi, na.rm = T)`,
         fi_l = `mean(fi, na.rm = T) - 1.96 * sd(fi, na.rm = T)/sqrt(587)`,
         fi_h = `mean(fi, na.rm = T) + 1.96 * sd(fi, na.rm = T)/sqrt(587)`,
         dy = `mean(dy, na.rm = T)`,
         dy_l = `mean(dy, na.rm = T) - 1.96 * sd(dy, na.rm = T)/sqrt(587)`,
         dy_h = `mean(dy, na.rm = T) + 1.96 * sd(dy, na.rm = T)/sqrt(587)`,
         pa = `mean(pa, na.rm = T)`,
         pa_l = `mean(pa, na.rm = T) - 1.96 * sd(pa, na.rm = T)/sqrt(587)`,
         pa_h = `mean(pa, na.rm = T) + 1.96 * sd(pa, na.rm = T)/sqrt(587)`,
         fa = `mean(fa, na.rm = T)`,
         fa_l = `mean(fa, na.rm = T) - 1.96 * sd(fa, na.rm = T)/sqrt(587)`,
         fa_h = `mean(fa, na.rm = T) + 1.96 * sd(fa, na.rm = T)/sqrt(587)`,
         rf = `mean(rf, na.rm = T)`,
         rf_l = `mean(rf, na.rm = T) - 1.96 * sd(rf, na.rm = T)/sqrt(587)`,
         rf_h = `mean(rf, na.rm = T) + 1.96 * sd(rf, na.rm = T)/sqrt(587)`,
         pf = `mean(pf, na.rm = T)`,
         pf_l = `mean(pf, na.rm = T) - 1.96 * sd(pf, na.rm = T)/sqrt(587)`,
         pf_h = `mean(pf, na.rm = T) + 1.96 * sd(pf, na.rm = T)/sqrt(587)`)

c2.spider <- c2.spider %>%
  rename(ql = `mean(ql, na.rm = T)`,
         ql_l = `mean(ql, na.rm = T) - 1.96 * sd(ql, na.rm = T)/sqrt(454)`,
         ql_h = `mean(ql, na.rm = T) + 1.96 * sd(ql, na.rm = T)/sqrt(454)`,
         ef = `mean(ef, na.rm = T)`,
         ef_l = `mean(ef, na.rm = T) - 1.96 * sd(ef, na.rm = T)/sqrt(454)`,
         ef_h = `mean(ef, na.rm = T) + 1.96 * sd(ef, na.rm = T)/sqrt(454)`,
         cf = `mean(cf, na.rm = T)`,
         cf_l = `mean(cf, na.rm = T) - 1.96 * sd(cf, na.rm = T)/sqrt(454)`,
         cf_h = `mean(cf, na.rm = T) + 1.96 * sd(cf, na.rm = T)/sqrt(454)`,
         sf = `mean(sf, na.rm = T)`,
         sf_l = `mean(sf, na.rm = T) - 1.96 * sd(sf, na.rm = T)/sqrt(454)`,
         sf_h = `mean(sf, na.rm = T) + 1.96 * sd(sf, na.rm = T)/sqrt(454)`,
         sl = `mean(sl, na.rm = T)`,
         sl_l = `mean(sl, na.rm = T) - 1.96 * sd(sl, na.rm = T)/sqrt(454)`,
         sl_h = `mean(sl, na.rm = T) + 1.96 * sd(sl, na.rm = T)/sqrt(454)`,
         nv = `mean(nv, na.rm = T)`,
         nv_l = `mean(nv, na.rm = T) - 1.96 * sd(nv, na.rm = T)/sqrt(454)`,
         nv_h = `mean(nv, na.rm = T) + 1.96 * sd(nv, na.rm = T)/sqrt(454)`,
         ap = `mean(ap, na.rm = T)`,
         ap_l = `mean(ap, na.rm = T) - 1.96 * sd(ap, na.rm = T)/sqrt(454)`,
         ap_h = `mean(ap, na.rm = T) + 1.96 * sd(ap, na.rm = T)/sqrt(454)`,
         co = `mean(co, na.rm = T)`,
         co_l = `mean(co, na.rm = T) - 1.96 * sd(co, na.rm = T)/sqrt(454)`,
         co_h = `mean(co, na.rm = T) + 1.96 * sd(co, na.rm = T)/sqrt(454)`,
         di = `mean(di, na.rm = T)`,
         di_l = `mean(di, na.rm = T) - 1.96 * sd(di, na.rm = T)/sqrt(454)`,
         di_h = `mean(di, na.rm = T) + 1.96 * sd(di, na.rm = T)/sqrt(454)`,
         fi = `mean(fi, na.rm = T)`,
         fi_l = `mean(fi, na.rm = T) - 1.96 * sd(fi, na.rm = T)/sqrt(454)`,
         fi_h = `mean(fi, na.rm = T) + 1.96 * sd(fi, na.rm = T)/sqrt(454)`,
         dy = `mean(dy, na.rm = T)`,
         dy_l = `mean(dy, na.rm = T) - 1.96 * sd(dy, na.rm = T)/sqrt(454)`,
         dy_h = `mean(dy, na.rm = T) + 1.96 * sd(dy, na.rm = T)/sqrt(454)`,
         pa = `mean(pa, na.rm = T)`,
         pa_l = `mean(pa, na.rm = T) - 1.96 * sd(pa, na.rm = T)/sqrt(454)`,
         pa_h = `mean(pa, na.rm = T) + 1.96 * sd(pa, na.rm = T)/sqrt(454)`,
         fa = `mean(fa, na.rm = T)`,
         fa_l = `mean(fa, na.rm = T) - 1.96 * sd(fa, na.rm = T)/sqrt(454)`,
         fa_h = `mean(fa, na.rm = T) + 1.96 * sd(fa, na.rm = T)/sqrt(454)`,
         rf = `mean(rf, na.rm = T)`,
         rf_l = `mean(rf, na.rm = T) - 1.96 * sd(rf, na.rm = T)/sqrt(454)`,
         rf_h = `mean(rf, na.rm = T) + 1.96 * sd(rf, na.rm = T)/sqrt(454)`,
         pf = `mean(pf, na.rm = T)`,
         pf_l = `mean(pf, na.rm = T) - 1.96 * sd(pf, na.rm = T)/sqrt(454)`,
         pf_h = `mean(pf, na.rm = T) + 1.96 * sd(pf, na.rm = T)/sqrt(454)`)

c3.spider <- c3.spider %>%
  rename(ql = `mean(ql, na.rm = T)`,
         ql_l = `mean(ql, na.rm = T) - 1.96 * sd(ql, na.rm = T)/sqrt(208)`,
         ql_h = `mean(ql, na.rm = T) + 1.96 * sd(ql, na.rm = T)/sqrt(208)`,
         ef = `mean(ef, na.rm = T)`,
         ef_l = `mean(ef, na.rm = T) - 1.96 * sd(ef, na.rm = T)/sqrt(208)`,
         ef_h = `mean(ef, na.rm = T) + 1.96 * sd(ef, na.rm = T)/sqrt(208)`,
         cf = `mean(cf, na.rm = T)`,
         cf_l = `mean(cf, na.rm = T) - 1.96 * sd(cf, na.rm = T)/sqrt(208)`,
         cf_h = `mean(cf, na.rm = T) + 1.96 * sd(cf, na.rm = T)/sqrt(208)`,
         sf = `mean(sf, na.rm = T)`,
         sf_l = `mean(sf, na.rm = T) - 1.96 * sd(sf, na.rm = T)/sqrt(208)`,
         sf_h = `mean(sf, na.rm = T) + 1.96 * sd(sf, na.rm = T)/sqrt(208)`,
         sl = `mean(sl, na.rm = T)`,
         sl_l = `mean(sl, na.rm = T) - 1.96 * sd(sl, na.rm = T)/sqrt(208)`,
         sl_h = `mean(sl, na.rm = T) + 1.96 * sd(sl, na.rm = T)/sqrt(208)`,
         nv = `mean(nv, na.rm = T)`,
         nv_l = `mean(nv, na.rm = T) - 1.96 * sd(nv, na.rm = T)/sqrt(208)`,
         nv_h = `mean(nv, na.rm = T) + 1.96 * sd(nv, na.rm = T)/sqrt(208)`,
         ap = `mean(ap, na.rm = T)`,
         ap_l = `mean(ap, na.rm = T) - 1.96 * sd(ap, na.rm = T)/sqrt(208)`,
         ap_h = `mean(ap, na.rm = T) + 1.96 * sd(ap, na.rm = T)/sqrt(208)`,
         co = `mean(co, na.rm = T)`,
         co_l = `mean(co, na.rm = T) - 1.96 * sd(co, na.rm = T)/sqrt(208)`,
         co_h = `mean(co, na.rm = T) + 1.96 * sd(co, na.rm = T)/sqrt(208)`,
         di = `mean(di, na.rm = T)`,
         di_l = `mean(di, na.rm = T) - 1.96 * sd(di, na.rm = T)/sqrt(208)`,
         di_h = `mean(di, na.rm = T) + 1.96 * sd(di, na.rm = T)/sqrt(208)`,
         fi = `mean(fi, na.rm = T)`,
         fi_l = `mean(fi, na.rm = T) - 1.96 * sd(fi, na.rm = T)/sqrt(208)`,
         fi_h = `mean(fi, na.rm = T) + 1.96 * sd(fi, na.rm = T)/sqrt(208)`,
         dy = `mean(dy, na.rm = T)`,
         dy_l = `mean(dy, na.rm = T) - 1.96 * sd(dy, na.rm = T)/sqrt(208)`,
         dy_h = `mean(dy, na.rm = T) + 1.96 * sd(dy, na.rm = T)/sqrt(208)`,
         pa = `mean(pa, na.rm = T)`,
         pa_l = `mean(pa, na.rm = T) - 1.96 * sd(pa, na.rm = T)/sqrt(208)`,
         pa_h = `mean(pa, na.rm = T) + 1.96 * sd(pa, na.rm = T)/sqrt(208)`,
         fa = `mean(fa, na.rm = T)`,
         fa_l = `mean(fa, na.rm = T) - 1.96 * sd(fa, na.rm = T)/sqrt(208)`,
         fa_h = `mean(fa, na.rm = T) + 1.96 * sd(fa, na.rm = T)/sqrt(208)`,
         rf = `mean(rf, na.rm = T)`,
         rf_l = `mean(rf, na.rm = T) - 1.96 * sd(rf, na.rm = T)/sqrt(208)`,
         rf_h = `mean(rf, na.rm = T) + 1.96 * sd(rf, na.rm = T)/sqrt(208)`,
         pf = `mean(pf, na.rm = T)`,
         pf_l = `mean(pf, na.rm = T) - 1.96 * sd(pf, na.rm = T)/sqrt(208)`,
         pf_h = `mean(pf, na.rm = T) + 1.96 * sd(pf, na.rm = T)/sqrt(208)`)

c4.spider <- c4.spider %>%
  rename(ql = `mean(ql, na.rm = T)`,
         ql_l = `mean(ql, na.rm = T) - 1.96 * sd(ql, na.rm = T)/sqrt(119)`,
         ql_h = `mean(ql, na.rm = T) + 1.96 * sd(ql, na.rm = T)/sqrt(119)`,
         ef = `mean(ef, na.rm = T)`,
         ef_l = `mean(ef, na.rm = T) - 1.96 * sd(ef, na.rm = T)/sqrt(119)`,
         ef_h = `mean(ef, na.rm = T) + 1.96 * sd(ef, na.rm = T)/sqrt(119)`,
         cf = `mean(cf, na.rm = T)`,
         cf_l = `mean(cf, na.rm = T) - 1.96 * sd(cf, na.rm = T)/sqrt(119)`,
         cf_h = `mean(cf, na.rm = T) + 1.96 * sd(cf, na.rm = T)/sqrt(119)`,
         sf = `mean(sf, na.rm = T)`,
         sf_l = `mean(sf, na.rm = T) - 1.96 * sd(sf, na.rm = T)/sqrt(119)`,
         sf_h = `mean(sf, na.rm = T) + 1.96 * sd(sf, na.rm = T)/sqrt(119)`,
         sl = `mean(sl, na.rm = T)`,
         sl_l = `mean(sl, na.rm = T) - 1.96 * sd(sl, na.rm = T)/sqrt(119)`,
         sl_h = `mean(sl, na.rm = T) + 1.96 * sd(sl, na.rm = T)/sqrt(119)`,
         nv = `mean(nv, na.rm = T)`,
         nv_l = `mean(nv, na.rm = T) - 1.96 * sd(nv, na.rm = T)/sqrt(119)`,
         nv_h = `mean(nv, na.rm = T) + 1.96 * sd(nv, na.rm = T)/sqrt(119)`,
         ap = `mean(ap, na.rm = T)`,
         ap_l = `mean(ap, na.rm = T) - 1.96 * sd(ap, na.rm = T)/sqrt(119)`,
         ap_h = `mean(ap, na.rm = T) + 1.96 * sd(ap, na.rm = T)/sqrt(119)`,
         co = `mean(co, na.rm = T)`,
         co_l = `mean(co, na.rm = T) - 1.96 * sd(co, na.rm = T)/sqrt(119)`,
         co_h = `mean(co, na.rm = T) + 1.96 * sd(co, na.rm = T)/sqrt(119)`,
         di = `mean(di, na.rm = T)`,
         di_l = `mean(di, na.rm = T) - 1.96 * sd(di, na.rm = T)/sqrt(119)`,
         di_h = `mean(di, na.rm = T) + 1.96 * sd(di, na.rm = T)/sqrt(119)`,
         fi = `mean(fi, na.rm = T)`,
         fi_l = `mean(fi, na.rm = T) - 1.96 * sd(fi, na.rm = T)/sqrt(119)`,
         fi_h = `mean(fi, na.rm = T) + 1.96 * sd(fi, na.rm = T)/sqrt(119)`,
         dy = `mean(dy, na.rm = T)`,
         dy_l = `mean(dy, na.rm = T) - 1.96 * sd(dy, na.rm = T)/sqrt(119)`,
         dy_h = `mean(dy, na.rm = T) + 1.96 * sd(dy, na.rm = T)/sqrt(119)`,
         pa = `mean(pa, na.rm = T)`,
         pa_l = `mean(pa, na.rm = T) - 1.96 * sd(pa, na.rm = T)/sqrt(119)`,
         pa_h = `mean(pa, na.rm = T) + 1.96 * sd(pa, na.rm = T)/sqrt(119)`,
         fa = `mean(fa, na.rm = T)`,
         fa_l = `mean(fa, na.rm = T) - 1.96 * sd(fa, na.rm = T)/sqrt(119)`,
         fa_h = `mean(fa, na.rm = T) + 1.96 * sd(fa, na.rm = T)/sqrt(119)`,
         rf = `mean(rf, na.rm = T)`,
         rf_l = `mean(rf, na.rm = T) - 1.96 * sd(rf, na.rm = T)/sqrt(119)`,
         rf_h = `mean(rf, na.rm = T) + 1.96 * sd(rf, na.rm = T)/sqrt(119)`,
         pf = `mean(pf, na.rm = T)`,
         pf_l = `mean(pf, na.rm = T) - 1.96 * sd(pf, na.rm = T)/sqrt(119)`,
         pf_h = `mean(pf, na.rm = T) + 1.96 * sd(pf, na.rm = T)/sqrt(119)`)

c5.spider <- c5.spider %>%
  rename(ql = `mean(ql, na.rm = T)`,
         ql_l = `mean(ql, na.rm = T) - 1.96 * sd(ql, na.rm = T)/sqrt(121)`,
         ql_h = `mean(ql, na.rm = T) + 1.96 * sd(ql, na.rm = T)/sqrt(121)`,
         ef = `mean(ef, na.rm = T)`,
         ef_l = `mean(ef, na.rm = T) - 1.96 * sd(ef, na.rm = T)/sqrt(121)`,
         ef_h = `mean(ef, na.rm = T) + 1.96 * sd(ef, na.rm = T)/sqrt(121)`,
         cf = `mean(cf, na.rm = T)`,
         cf_l = `mean(cf, na.rm = T) - 1.96 * sd(cf, na.rm = T)/sqrt(121)`,
         cf_h = `mean(cf, na.rm = T) + 1.96 * sd(cf, na.rm = T)/sqrt(121)`,
         sf = `mean(sf, na.rm = T)`,
         sf_l = `mean(sf, na.rm = T) - 1.96 * sd(sf, na.rm = T)/sqrt(121)`,
         sf_h = `mean(sf, na.rm = T) + 1.96 * sd(sf, na.rm = T)/sqrt(121)`,
         sl = `mean(sl, na.rm = T)`,
         sl_l = `mean(sl, na.rm = T) - 1.96 * sd(sl, na.rm = T)/sqrt(121)`,
         sl_h = `mean(sl, na.rm = T) + 1.96 * sd(sl, na.rm = T)/sqrt(121)`,
         nv = `mean(nv, na.rm = T)`,
         nv_l = `mean(nv, na.rm = T) - 1.96 * sd(nv, na.rm = T)/sqrt(121)`,
         nv_h = `mean(nv, na.rm = T) + 1.96 * sd(nv, na.rm = T)/sqrt(121)`,
         ap = `mean(ap, na.rm = T)`,
         ap_l = `mean(ap, na.rm = T) - 1.96 * sd(ap, na.rm = T)/sqrt(121)`,
         ap_h = `mean(ap, na.rm = T) + 1.96 * sd(ap, na.rm = T)/sqrt(121)`,
         co = `mean(co, na.rm = T)`,
         co_l = `mean(co, na.rm = T) - 1.96 * sd(co, na.rm = T)/sqrt(121)`,
         co_h = `mean(co, na.rm = T) + 1.96 * sd(co, na.rm = T)/sqrt(121)`,
         di = `mean(di, na.rm = T)`,
         di_l = `mean(di, na.rm = T) - 1.96 * sd(di, na.rm = T)/sqrt(121)`,
         di_h = `mean(di, na.rm = T) + 1.96 * sd(di, na.rm = T)/sqrt(121)`,
         fi = `mean(fi, na.rm = T)`,
         fi_l = `mean(fi, na.rm = T) - 1.96 * sd(fi, na.rm = T)/sqrt(121)`,
         fi_h = `mean(fi, na.rm = T) + 1.96 * sd(fi, na.rm = T)/sqrt(121)`,
         dy = `mean(dy, na.rm = T)`,
         dy_l = `mean(dy, na.rm = T) - 1.96 * sd(dy, na.rm = T)/sqrt(121)`,
         dy_h = `mean(dy, na.rm = T) + 1.96 * sd(dy, na.rm = T)/sqrt(121)`,
         pa = `mean(pa, na.rm = T)`,
         pa_l = `mean(pa, na.rm = T) - 1.96 * sd(pa, na.rm = T)/sqrt(121)`,
         pa_h = `mean(pa, na.rm = T) + 1.96 * sd(pa, na.rm = T)/sqrt(121)`,
         fa = `mean(fa, na.rm = T)`,
         fa_l = `mean(fa, na.rm = T) - 1.96 * sd(fa, na.rm = T)/sqrt(121)`,
         fa_h = `mean(fa, na.rm = T) + 1.96 * sd(fa, na.rm = T)/sqrt(121)`,
         rf = `mean(rf, na.rm = T)`,
         rf_l = `mean(rf, na.rm = T) - 1.96 * sd(rf, na.rm = T)/sqrt(121)`,
         rf_h = `mean(rf, na.rm = T) + 1.96 * sd(rf, na.rm = T)/sqrt(121)`,
         pf = `mean(pf, na.rm = T)`,
         pf_l = `mean(pf, na.rm = T) - 1.96 * sd(pf, na.rm = T)/sqrt(121)`,
         pf_h = `mean(pf, na.rm = T) + 1.96 * sd(pf, na.rm = T)/sqrt(121)`)


c1 <- c1.spider %>%
  select(ql, ef, cf, sf, sl, nv, ap, co, di, fi, dy, pa, fa, rf, pf)
c1_l <- c1.spider %>%
  select(ql_l, ef_l, cf_l, sf_l, sl_l, nv_l, ap_l, co_l, di_l, fi_l, dy_l, pa_l, fa_l, rf_l, pf_l)
c1_l <- c1_l %>%
  rename(ql=ql_l, ef=ef_l, cf=cf_l, sf=sf_l, sl=sl_l, nv=nv_l, ap=ap_l, co=co_l, di=di_l, 
         fi=fi_l, dy=dy_l, pa=pa_l, fa=fa_l, rf=rf_l, pf=pf_l)
c1_h <- c1.spider %>%
  select(ql_h, ef_h, cf_h, sf_h, sl_h, nv_h, ap_h, co_h, di_h, fi_h, dy_h, pa_h, fa_h, rf_h, pf_h)
c1_h <- c1_h %>%
  rename(ql=ql_h, ef=ef_h, cf=cf_h, sf=sf_h, sl=sl_h, nv=nv_h, ap=ap_h, co=co_h, di=di_h, 
         fi=fi_h, dy=dy_h, pa=pa_h, fa=fa_h, rf=rf_h, pf=pf_h)

c1_final <- rbind(c1, rbind(c1_l, c1_h))


c2 <- c2.spider %>%
  select(ql, ef, cf, sf, sl, nv, ap, co, di, fi, dy, pa, fa, rf, pf)
c2_l <- c2.spider %>%
  select(ql_l, ef_l, cf_l, sf_l, sl_l, nv_l, ap_l, co_l, di_l, fi_l, dy_l, pa_l, fa_l, rf_l, pf_l)
c2_l <- c2_l %>%
  rename(ql=ql_l, ef=ef_l, cf=cf_l, sf=sf_l, sl=sl_l, nv=nv_l, ap=ap_l, co=co_l, di=di_l, 
         fi=fi_l, dy=dy_l, pa=pa_l, fa=fa_l, rf=rf_l, pf=pf_l)
c2_h <- c2.spider %>%
  select(ql_h, ef_h, cf_h, sf_h, sl_h, nv_h, ap_h, co_h, di_h, fi_h, dy_h, pa_h, fa_h, rf_h, pf_h)
c2_h <- c2_h %>%
  rename(ql=ql_h, ef=ef_h, cf=cf_h, sf=sf_h, sl=sl_h, nv=nv_h, ap=ap_h, co=co_h, di=di_h, 
         fi=fi_h, dy=dy_h, pa=pa_h, fa=fa_h, rf=rf_h, pf=pf_h)

c2_final <- rbind(c2, rbind(c2_l, c2_h))


c3 <- c3.spider %>%
  select(ql, ef, cf, sf, sl, nv, ap, co, di, fi, dy, pa, fa, rf, pf)
c3_l <- c3.spider %>%
  select(ql_l, ef_l, cf_l, sf_l, sl_l, nv_l, ap_l, co_l, di_l, fi_l, dy_l, pa_l, fa_l, rf_l, pf_l)
c3_l <- c3_l %>%
  rename(ql=ql_l, ef=ef_l, cf=cf_l, sf=sf_l, sl=sl_l, nv=nv_l, ap=ap_l, co=co_l, di=di_l, 
         fi=fi_l, dy=dy_l, pa=pa_l, fa=fa_l, rf=rf_l, pf=pf_l)
c3_h <- c3.spider %>%
  select(ql_h, ef_h, cf_h, sf_h, sl_h, nv_h, ap_h, co_h, di_h, fi_h, dy_h, pa_h, fa_h, rf_h, pf_h)
c3_h <- c3_h %>%
  rename(ql=ql_h, ef=ef_h, cf=cf_h, sf=sf_h, sl=sl_h, nv=nv_h, ap=ap_h, co=co_h, di=di_h, 
         fi=fi_h, dy=dy_h, pa=pa_h, fa=fa_h, rf=rf_h, pf=pf_h)

c3_final <- rbind(c3, rbind(c3_l, c3_h))


c4 <- c4.spider %>%
  select(ql, ef, cf, sf, sl, nv, ap, co, di, fi, dy, pa, fa, rf, pf)
c4_l <- c4.spider %>%
  select(ql_l, ef_l, cf_l, sf_l, sl_l, nv_l, ap_l, co_l, di_l, fi_l, dy_l, pa_l, fa_l, rf_l, pf_l)
c4_l <- c4_l %>%
  rename(ql=ql_l, ef=ef_l, cf=cf_l, sf=sf_l, sl=sl_l, nv=nv_l, ap=ap_l, co=co_l, di=di_l, 
         fi=fi_l, dy=dy_l, pa=pa_l, fa=fa_l, rf=rf_l, pf=pf_l)
c4_h <- c4.spider %>%
  select(ql_h, ef_h, cf_h, sf_h, sl_h, nv_h, ap_h, co_h, di_h, fi_h, dy_h, pa_h, fa_h, rf_h, pf_h)
c4_h <- c4_h %>%
  rename(ql=ql_h, ef=ef_h, cf=cf_h, sf=sf_h, sl=sl_h, nv=nv_h, ap=ap_h, co=co_h, di=di_h, 
         fi=fi_h, dy=dy_h, pa=pa_h, fa=fa_h, rf=rf_h, pf=pf_h)

c4_final <- rbind(c4, rbind(c4_l, c4_h))


c5 <- c5.spider %>%
  select(ql, ef, cf, sf, sl, nv, ap, co, di, fi, dy, pa, fa, rf, pf)
c5_l <- c5.spider %>%
  select(ql_l, ef_l, cf_l, sf_l, sl_l, nv_l, ap_l, co_l, di_l, fi_l, dy_l, pa_l, fa_l, rf_l, pf_l)
c5_l <- c5_l %>%
  rename(ql=ql_l, ef=ef_l, cf=cf_l, sf=sf_l, sl=sl_l, nv=nv_l, ap=ap_l, co=co_l, di=di_l, 
         fi=fi_l, dy=dy_l, pa=pa_l, fa=fa_l, rf=rf_l, pf=pf_l)
c5_h <- c5.spider %>%
  select(ql_h, ef_h, cf_h, sf_h, sl_h, nv_h, ap_h, co_h, di_h, fi_h, dy_h, pa_h, fa_h, rf_h, pf_h)
c5_h <- c5_h %>%
  rename(ql=ql_h, ef=ef_h, cf=cf_h, sf=sf_h, sl=sl_h, nv=nv_h, ap=ap_h, co=co_h, di=di_h, 
         fi=fi_h, dy=dy_h, pa=pa_h, fa=fa_h, rf=rf_h, pf=pf_h)

c5_final <- rbind(c5, rbind(c5_l, c5_h))


final.spider <- rbind(c1_final, rbind(c2_final, rbind(c3_final, rbind(c4_final, c5_final))))

final.spider <- as.data.frame(final.spider)
final.spider$ap[2] <- 0

final.spider <- final.spider %>%
  rename(sl1=sl, nv1=nv, ap1=ap, co1=co, di1=di, fi1=fi, dy1=dy, pa1=pa, fa1=fa)

final.spider$sl[1] <- 100-final.spider$sl1[1]
final.spider$sl[2] <- 100-final.spider$sl1[3]
final.spider$sl[3] <- 100-final.spider$sl1[2]
final.spider$sl[4] <- 100-final.spider$sl1[4]
final.spider$sl[5] <- 100-final.spider$sl1[6]
final.spider$sl[6] <- 100-final.spider$sl1[5]
final.spider$sl[7] <- 100-final.spider$sl1[7]
final.spider$sl[8] <- 100-final.spider$sl1[9]
final.spider$sl[9] <- 100-final.spider$sl1[8]
final.spider$sl[10] <- 100-final.spider$sl1[10]
final.spider$sl[11] <- 100-final.spider$sl1[12]
final.spider$sl[12] <- 100-final.spider$sl1[11]
final.spider$sl[13] <- 100-final.spider$sl1[13]
final.spider$sl[14] <- 100-final.spider$sl1[15]
final.spider$sl[15] <- 100-final.spider$sl1[14]

final.spider$nv[1] <- 100-final.spider$nv1[1]
final.spider$nv[2] <- 100-final.spider$nv1[3]
final.spider$nv[3] <- 100-final.spider$nv1[2]
final.spider$nv[4] <- 100-final.spider$nv1[4]
final.spider$nv[5] <- 100-final.spider$nv1[6]
final.spider$nv[6] <- 100-final.spider$nv1[5]
final.spider$nv[7] <- 100-final.spider$nv1[7]
final.spider$nv[8] <- 100-final.spider$nv1[9]
final.spider$nv[9] <- 100-final.spider$nv1[8]
final.spider$nv[10] <- 100-final.spider$nv1[10]
final.spider$nv[11] <- 100-final.spider$nv1[12]
final.spider$nv[12] <- 100-final.spider$nv1[11]
final.spider$nv[13] <- 100-final.spider$nv1[13]
final.spider$nv[14] <- 100-final.spider$nv1[15]
final.spider$nv[15] <- 100-final.spider$nv1[14]

final.spider$ap[1] <- 100-final.spider$ap1[1]
final.spider$ap[2] <- 100-final.spider$ap1[3]
final.spider$ap[3] <- 100-final.spider$ap1[2]
final.spider$ap[4] <- 100-final.spider$ap1[4]
final.spider$ap[5] <- 100-final.spider$ap1[6]
final.spider$ap[6] <- 100-final.spider$ap1[5]
final.spider$ap[7] <- 100-final.spider$ap1[7]
final.spider$ap[8] <- 100-final.spider$ap1[9]
final.spider$ap[9] <- 100-final.spider$ap1[8]
final.spider$ap[10] <- 100-final.spider$ap1[10]
final.spider$ap[11] <- 100-final.spider$ap1[12]
final.spider$ap[12] <- 100-final.spider$ap1[11]
final.spider$ap[13] <- 100-final.spider$ap1[13]
final.spider$ap[14] <- 100-final.spider$ap1[15]
final.spider$ap[15] <- 100-final.spider$ap1[14]

final.spider$co[1] <- 100-final.spider$co1[1]
final.spider$co[2] <- 100-final.spider$co1[3]
final.spider$co[3] <- 100-final.spider$co1[2]
final.spider$co[4] <- 100-final.spider$co1[4]
final.spider$co[5] <- 100-final.spider$co1[6]
final.spider$co[6] <- 100-final.spider$co1[5]
final.spider$co[7] <- 100-final.spider$co1[7]
final.spider$co[8] <- 100-final.spider$co1[9]
final.spider$co[9] <- 100-final.spider$co1[8]
final.spider$co[10] <- 100-final.spider$co1[10]
final.spider$co[11] <- 100-final.spider$co1[12]
final.spider$co[12] <- 100-final.spider$co1[11]
final.spider$co[13] <- 100-final.spider$co1[13]
final.spider$co[14] <- 100-final.spider$co1[15]
final.spider$co[15] <- 100-final.spider$co1[14]

final.spider$di[1] <- 100-final.spider$di1[1]
final.spider$di[2] <- 100-final.spider$di1[3]
final.spider$di[3] <- 100-final.spider$di1[2]
final.spider$di[4] <- 100-final.spider$di1[4]
final.spider$di[5] <- 100-final.spider$di1[6]
final.spider$di[6] <- 100-final.spider$di1[5]
final.spider$di[7] <- 100-final.spider$di1[7]
final.spider$di[8] <- 100-final.spider$di1[9]
final.spider$di[9] <- 100-final.spider$di1[8]
final.spider$di[10] <- 100-final.spider$di1[10]
final.spider$di[11] <- 100-final.spider$di1[12]
final.spider$di[12] <- 100-final.spider$di1[11]
final.spider$di[13] <- 100-final.spider$di1[13]
final.spider$di[14] <- 100-final.spider$di1[15]
final.spider$di[15] <- 100-final.spider$di1[14]

final.spider$fi[1] <- 100-final.spider$fi1[1]
final.spider$fi[2] <- 100-final.spider$fi1[3]
final.spider$fi[3] <- 100-final.spider$fi1[2]
final.spider$fi[4] <- 100-final.spider$fi1[4]
final.spider$fi[5] <- 100-final.spider$fi1[6]
final.spider$fi[6] <- 100-final.spider$fi1[5]
final.spider$fi[7] <- 100-final.spider$fi1[7]
final.spider$fi[8] <- 100-final.spider$fi1[9]
final.spider$fi[9] <- 100-final.spider$fi1[8]
final.spider$fi[10] <- 100-final.spider$fi1[10]
final.spider$fi[11] <- 100-final.spider$fi1[12]
final.spider$fi[12] <- 100-final.spider$fi1[11]
final.spider$fi[13] <- 100-final.spider$fi1[13]
final.spider$fi[14] <- 100-final.spider$fi1[15]
final.spider$fi[15] <- 100-final.spider$fi1[14]

final.spider$dy[1] <- 100-final.spider$dy1[1]
final.spider$dy[2] <- 100-final.spider$dy1[3]
final.spider$dy[3] <- 100-final.spider$dy1[2]
final.spider$dy[4] <- 100-final.spider$dy1[4]
final.spider$dy[5] <- 100-final.spider$dy1[6]
final.spider$dy[6] <- 100-final.spider$dy1[5]
final.spider$dy[7] <- 100-final.spider$dy1[7]
final.spider$dy[8] <- 100-final.spider$dy1[9]
final.spider$dy[9] <- 100-final.spider$dy1[8]
final.spider$dy[10] <- 100-final.spider$dy1[10]
final.spider$dy[11] <- 100-final.spider$dy1[12]
final.spider$dy[12] <- 100-final.spider$dy1[11]
final.spider$dy[13] <- 100-final.spider$dy1[13]
final.spider$dy[14] <- 100-final.spider$dy1[15]
final.spider$dy[15] <- 100-final.spider$dy1[14]

final.spider$pa[1] <- 100-final.spider$pa1[1]
final.spider$pa[2] <- 100-final.spider$pa1[3]
final.spider$pa[3] <- 100-final.spider$pa1[2]
final.spider$pa[4] <- 100-final.spider$pa1[4]
final.spider$pa[5] <- 100-final.spider$pa1[6]
final.spider$pa[6] <- 100-final.spider$pa1[5]
final.spider$pa[7] <- 100-final.spider$pa1[7]
final.spider$pa[8] <- 100-final.spider$pa1[9]
final.spider$pa[9] <- 100-final.spider$pa1[8]
final.spider$pa[10] <- 100-final.spider$pa1[10]
final.spider$pa[11] <- 100-final.spider$pa1[12]
final.spider$pa[12] <- 100-final.spider$pa1[11]
final.spider$pa[13] <- 100-final.spider$pa1[13]
final.spider$pa[14] <- 100-final.spider$pa1[15]
final.spider$pa[15] <- 100-final.spider$pa1[14]

final.spider$fa[1] <- 100-final.spider$fa1[1]
final.spider$fa[2] <- 100-final.spider$fa1[3]
final.spider$fa[3] <- 100-final.spider$fa1[2]
final.spider$fa[4] <- 100-final.spider$fa1[4]
final.spider$fa[5] <- 100-final.spider$fa1[6]
final.spider$fa[6] <- 100-final.spider$fa1[5]
final.spider$fa[7] <- 100-final.spider$fa1[7]
final.spider$fa[8] <- 100-final.spider$fa1[9]
final.spider$fa[9] <- 100-final.spider$fa1[8]
final.spider$fa[10] <- 100-final.spider$fa1[10]
final.spider$fa[11] <- 100-final.spider$fa1[12]
final.spider$fa[12] <- 100-final.spider$fa1[11]
final.spider$fa[13] <- 100-final.spider$fa1[13]
final.spider$fa[14] <- 100-final.spider$fa1[15]
final.spider$fa[15] <- 100-final.spider$fa1[14]

final.spider <- final.spider %>%
  select(ql, sf, cf, ef, sl, nv, ap, co, di, fi, dy, pa, fa, rf, pf)

span <- rbind(rep(100, 15), rep(0, 15))
span <- as.data.frame(span)
span <- span %>%
  rename(ql = V1,
         sf = V2,
         cf = V3,
         ef = V4,
         sl = V5,
         nv = V6,
         ap = V7,
         co = V8,
         di = V9,
         fi = V10,
         dy = V11,
         pa = V12,
         fa = V13,
         rf = V14,
         pf = V15)
final.spider100 <- rbind(span, final.spider)

final.spider100 <- as.data.frame(final.spider100)



mynewradarchart <- function (df, axistype = 0, seg = 3, pty = 32, 
                             pcol = c(alpha("#E41A1C", 1), NA, NA, 
                                      alpha("#377EB8", 1), NA, NA, 
                                      alpha("#4DAF4A", 1), NA, NA, 
                                      alpha("#984EA3", 1), NA, NA), 
                             plty = 1, plwd = 2, pdensity = NULL, pangle = 45, pfcol = NA, cglty = 1,
                             cglwd = 1, cglcol = "grey", axislabcol = "grey", 
                             acols = c(alpha("#E41A1C", .6), alpha("#377EB8", .6),
                                       alpha("#4DAF4A", .6), alpha("#984EA3", .6)), 
                             CI = F, title = "", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE,
                             vlabels = NULL, vlcex = NULL, caxislabels = NULL, calcex = NULL, 
                             paxislabels = NULL, palcex = NULL, legendtext = NULL, ...) 
{
  if (!is.data.frame(df)) {
    cat("The data must be given as dataframe.\n")
    return()
  }
  if ((n <- length(df)) < 3) {
    cat("The number of variables must be 3 or more.\n")
    return()
  }
  if (maxmin == FALSE) {
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type = "n", frame.plot = FALSE, 
       axes = FALSE, xlab = "", ylab = "", main = title, asp = 1, ...)
  theta <- seq(90, 450, length = n + 1) * pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) {
    polygon(xx * (i + CGap)/(seg + CGap), yy * (i + CGap)/(seg +  CGap), 
            lty = cglty, lwd = cglwd, border = cglcol)
    if (axistype == 1 | axistype == 3) 
      CAXISLABELS <- paste(i/seg * 100, "(%)")
    if (axistype == 4 | axistype == 5) 
      CAXISLABELS <- sprintf("%3.2f", i/seg)
    if (!is.null(caxislabels) & (i < length(caxislabels))) 
      CAXISLABELS <- caxislabels[i + 1]
    if (axistype == 1 | axistype == 3 | axistype == 4 | axistype == 5) {
      if (is.null(calcex)) 
        text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS, 
             col = axislabcol)
      else text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS, 
                col = axislabcol, cex = calcex)
    }
  }
  if (centerzero) {
    arrows(0, 0, xx * 1, yy * 1, lwd = cglwd, lty = cglty, 
           length = 0, col = cglcol)
  }
  else {
    arrows(xx/(seg + CGap), yy/(seg + CGap), xx * 1, yy * 1, 
           lwd = cglwd, lty = cglty, length = 0, col = cglcol)
  }
  PAXISLABELS <- df[1, 1:n]
  if (!is.null(paxislabels)) 
    PAXISLABELS <- paxislabels
  if (axistype == 2 | axistype == 3 | axistype == 5) {
    if (is.null(palcex)) 
      text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol)
    else text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol, cex = palcex)
  }
  VLABELS <- colnames(df)
  if (!is.null(vlabels)) 
    VLABELS <- vlabels
  if (is.null(vlcex)) 
    text(xx * 1.2, yy * 1.2, VLABELS)
  else text(xx * 1.2, yy * 1.2, VLABELS, cex = vlcex)
  series <- length(df[[1]])
  SX <- series - 2
  if (length(pty) < SX) {
    ptys <- rep(pty, SX)
  }
  else {
    ptys <- pty
  }
  if (length(pcol) < SX) {
    pcols <- rep(pcol, SX)
  }
  else {
    pcols <- pcol
  }
  if (length(plty) < SX) {
    pltys <- rep(plty, SX)
  }
  else {
    pltys <- plty
  }
  if (length(plwd) < SX) {
    plwds <- rep(plwd, SX)
  }
  else {
    plwds <- plwd
  }
  if (length(pdensity) < SX) {
    pdensities <- rep(pdensity, SX)
  }
  else {
    pdensities <- pdensity
  }
  if (length(pangle) < SX) {
    pangles <- rep(pangle, SX)
  }
  else {
    pangles <- pangle
  }
  if (length(pfcol) < SX) {
    pfcols <- rep(pfcol, SX)
  }
  else {
    pfcols <- pfcol
  }
  xxs <- matrix(NA, series-2, n)
  yys <- matrix(NA, series-2, n)
  for (i in 3:series) {
    k <- i-2
    xxs[k, ] <- xx
    yys[k, ] <- yy
    scale <- CGap/(seg + CGap) + (df[i, ] - df[2, ])/(df[1, ] - df[2, ]) * seg/(seg + CGap)
    if (sum(!is.na(df[i, ])) < 3) {
      cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n", i, df[i, ]))
    }
    else {
      for (j in 1:n) {
        if (is.na(df[i, j])) {
          if (na.itp) {
            left <- ifelse(j > 1, j - 1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left > 1, left - 1, n)
            }
            right <- ifelse(j < n, j + 1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right < n, right + 1, 1)
            }
            xxleft <- xx[left] * CGap/(seg + CGap) + 
              xx[left] * (df[i, left] - df[2, left])/(df[1, left] - df[2, left]) * seg/(seg + CGap)
            yyleft <- yy[left] * CGap/(seg + CGap) + 
              yy[left] * (df[i, left] - df[2, left])/(df[1, left] - df[2, left]) * seg/(seg + CGap)
            xxright <- xx[right] * CGap/(seg + CGap) + 
              xx[right] * (df[i, right] - df[2, right])/(df[1, right] - df[2, right]) * seg/(seg + CGap)
            yyright <- yy[right] * CGap/(seg + CGap) + 
              yy[right] * (df[i, right] - df[2, right])/(df[1, right] - df[2, right]) * seg/(seg + CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft
              yytmp <- yyleft
              xxleft <- xxright
              yyleft <- yyright
              xxright <- xxtmp
              yyright <- yytmp
            }
            xxs[k, j] <- xx[j] * (yyleft * xxright - yyright * xxleft)/
              (yy[j] * (xxright - xxleft) - xx[j] * (yyright - yyleft))
            yys[k, j] <- (yy[j]/xx[j]) * xxs[k, j]
          }
          else {
            xxs[k, j] <- 0
            yys[k, j] <- 0
          }
        }
        else {
          xxs[k, j] <- xx[j] * CGap/(seg + CGap) + xx[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, j]) * 
            seg/(seg + CGap)
          yys[k, j] <- yy[j] * CGap/(seg + CGap) + yy[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, j]) * 
            seg/(seg + CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs[k, ], yys[k, ], lty = pltys[i - 2], lwd = plwds[i - 2], 
                border = pcols[i - 2], col = pfcols[i - 2])
      }
      else {
        polygon(xxs[k, ], yys[k, ], lty = pltys[i - 2], lwd = plwds[i - 2], 
                border = pcols[i - 2], density = pdensities[i - 2], 
                angle = pangles[i - 2], col = pfcols[i - 2])
      }
      points(xx * scale, yy * scale, pch = ptys[i - 2], 
             col = pcols[i - 2])
    }
  }
  xxs <- cbind(xxs, xxs[, 1])
  yys <- cbind(yys, yys[, 1])
  if(CI){
    acols2 <- rep(acols, each = 3)
    nr <- c(1:nrow(xxs))
    for (l in sort(nr[seq(3, length(xxs), 3)])){
      polygon(c(xxs[l, ], xxs[l-1, ]), c(yys[l, ], yys[l-1, ]), lty = 1, lwd = 1, 
              border = NA, density = NA, col = acols2[l])
    }
  }
  if(!is.null(legendtext)){
    legend(x = 1.5, y = 1.2, legend = legendtext,
           bty = "n", pch = 20, col = acols,
           text.col = "black", cex = 1, pt.cex = 3, y.intersp = .8)
  }
}


png("rplot.png", width = 2900, height = 2100, res = 300)
mynewradarchart(final.spider100, CI = T,
                cglcol = "grey", cglty = 1, cglwd = 0.8,
                caxislabels = c("0","25","50","75", "100"), axistype = 1, axislabcol = "grey", 
                seg = 4,
                pcol = c(alpha("#00468BFF",.6), NA, NA,
                         alpha("#ED0000FF",.6), NA, NA,
                         alpha("#42B540FF",.6), NA, NA,
                         alpha("#0099B4FF",.6), NA, NA,
                         alpha("#925E9FFF",.6), NA, NA), 
                acols = c(alpha("#00468BFF",.6),
                          alpha("#ED0000FF",.6),
                          alpha("#42B540FF",.6),
                          alpha("#0099B4FF",.6),
                          alpha("#925E9FFF",.6)),
                plwd = 4,
                pty = 32, plty=1, 
                vlabels = c("global
health status",
                            "social
functioning",
                            "cognitive
functioning",
                            "emotional
functioning",
                            "insomnia",
                            "nausea/
vomiting",
                            "appetite
loss",
                            "constipation",
                            "diarrhea",
                            "financial
problems",
                            "dyspnea",
                            "pain",
                            "fatigue",
                            "role
functioning",
                            "physical
functioning"), vlcex = 1, legendtext = c("class 1", "class 2", "class 3", "class 4", "class 5"))
dev.off()







# Figure 3.
var1 <- rep(1:5, each=5)
var2 <- rep(1:5, 5)
var3 <- c(.869, .106, .000, .024, .014,
          .111, .768, .219, .197, .051,
          .000, .087, .692, .055, .154,
          .019, .036, .030, .671, .025,
          .001, .004, .059, .054, .758)

mydata <- cbind(var1, var2)
mydata <- cbind(mydata, var3)
mydata <- as.data.frame(mydata)


plot1 <- ggplot(mydata, aes(var1, var2, fill= var3)) + 
  geom_tile() +
  geom_text(aes(label = c(".87", ".11", ".00", ".02", ".01",
                          ".11", ".77", ".22", ".20", ".05",
                          ".00", ".09", ".69", ".06", ".15",
                          ".02", ".04", ".03", ".67", ".03",
                          ".00", ".00", ".06", ".05", ".76"))) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_x_continuous(position = "top") +
  scale_y_reverse() +
  labs(x = "follow-up", y = "baseline") +
  coord_equal() +
  theme_ipsum() +
  theme(legend.title = element_blank(), legend.key.height = unit(3.7, "line"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), axis.text = element_text(colour = "black"),
        axis.title.x = element_text(size = 12, hjust = .95),
        axis.title.y = element_text(size = 12, vjust = 2, hjust = .95))

setwd("adjust wd")
png("rplot.png", width = 1800, height = 1800, res = 300)
plot1
dev.off()
