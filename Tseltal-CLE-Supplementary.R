################################################################################
####          RANDOM BASELINE SUPPLEMENTARY MODELS          ####
####                    ODS                    ####
# gaussian glm of log transformed DV (w/ zeroes) #----
ods.rand.glmr <- glmmTMB(log(ods_mph+1) ~
                           tchiyr.std + mated.bi.std +
                           hsz.std + nsk.std +
                           I(stthr.std^2) +
                           tchiyr.std:mated.bi.std +
                           tchiyr.std:nsk.std +
                           tchiyr.std:hsz.std +
                           tchiyr.std:I(stthr.std^2) +
                           (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                         data=quantity.rand)
#summary(ods.rand.glmr)
#nsk.std                    1.44765    0.10109  14.320  < 2e-16 ***
#I(stthr.std^2)             1.84545    1.06075   1.740  0.08190 .  
#tchiyr.std:hsz.std         0.36835    0.16802   2.192  0.02836 *  
#tchiyr.std:I(stthr.std^2)  3.01524    1.15157   2.618  0.00884 **
ods.rand.glmr.resid <- residuals(ods.rand.glmr)
ods.rand.glmr.resid.plot <- qqplot.data(ods.rand.glmr.resid)
ods.rand.glmr.resid.resids <- data.frame(resid = ods.rand.glmr.resid)
ods.rand.glmr.resid.dist.plot <- ggplot(ods.rand.glmr.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)

# gaussian glm of log transformed DV (w/o zeroes) #----
ods.rand.glmr.nz <- glmmTMB(log(ods_mph+1) ~
                              tchiyr.std + mated.bi.std +
                              hsz.std + nsk.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bi.std +
                              tchiyr.std:nsk.std +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
                              (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                            data=subset(quantity.rand, ods_mph > 0))
#summary(ods.rand.glmr.nz)
#nsk.std                    1.03287    0.12910   8.000 1.24e-15 ***
#tchiyr.std:I(stthr.std^2)  3.21939    1.13552   2.835  0.00458 **
ods.rand.glmr.nz.resid <- residuals(ods.rand.glmr.nz)
ods.rand.glmr.nz.resid.plot <- qqplot.data(ods.rand.glmr.nz.resid)
ods.rand.glmr.nz.resid.resids <- data.frame(resid = ods.rand.glmr.nz.resid)
ods.rand.glmr.nz.resid.dist.plot <- ggplot(ods.rand.glmr.nz.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)

#AICtab(ods.rand.zinb, ods.rand.glmr, ods.rand.glmr.nz)

# binomial glm of zero vs. non-zero cases #----
ods.rand.bglm.nz <- glmmTMB(ods_mph.nz ~
                              tchiyr.std +
                              hsz.std + nsk.std +
                              I(stthr.std^2) +
#                              tchiyr.std:nsk.std + # no cnvgnce with this term
#                              tchiyr.std:hsz.std + # no cnvgnce with this term
#                              tchiyr.std:I(stthr.std^2) + # no cnvgnce with this term
#                              I(stthr.std^2):hsz.std + # no cnvgnce with this term
                              I(stthr.std^2):nsk.std +
                              (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                            data=quantity.rand,
                            family = "binomial")
#summary(ods.rand.bglm.nz)
# no significant effects


####                    TDS                    ####
# gaussian glm of log transformed DV (w/ zeroes) #----
tds.rand.glmr <- glmmTMB(log(tds_mph+1) ~
                           tchiyr.std + mated.bi.std +
                           hsz.std + nsk.std +
                           I(stthr.std^2) +
                           tchiyr.std:mated.bi.std +
                           tchiyr.std:nsk.std +
                           tchiyr.std:hsz.std +
                           tchiyr.std:I(stthr.std^2) +
                           (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                         data=quantity.rand)
#summary(tds.rand.glmr)
#tchiyr.std                 0.38960    0.19482   2.000  0.04552 *  
#nsk.std                    0.27780    0.12386   2.243  0.02491 *  
#tchiyr.std:nsk.std         0.33295    0.15543   2.142  0.03218 *  
#tchiyr.std:I(stthr.std^2) -4.39596    1.41272  -3.112  0.00186 ** 
tds.rand.glmr.resid <- residuals(tds.rand.glmr)
tds.rand.glmr.resid.plot <- qqplot.data(tds.rand.glmr.resid)
tds.rand.glmr.resid.resids <- data.frame(resid = tds.rand.glmr.resid)
tds.rand.glmr.resid.dist.plot <- ggplot(tds.rand.glmr.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)
# okay residuals, but a bit skewed

# gaussian glm of log transformed DV (w/o zeroes) #----
tds.rand.glmr.nz <- glmmTMB(log(tds_mph+1) ~
                              tchiyr.std + mated.bi.std +
                              hsz.std + nsk.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bi.std +
                              tchiyr.std:nsk.std +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
                             (1+I(stthr.std^2)|aclew_child_id),
                            data=subset(quantity.rand, tds_mph > 0))
#summary(tds.rand.glmr.nz)
#tchiyr.std                 0.42272    0.23051   1.834  0.06669 .  
#mated.bi.std6+             0.47091    0.24874   1.893  0.05833 .  
#nsk.std                   -0.26152    0.13595  -1.924  0.05440 .  
#tchiyr.std:nsk.std         0.41451    0.16041   2.584  0.00976 ** 
#tchiyr.std:I(stthr.std^2) -2.85718    1.46193  -1.954  0.05066 .  
tds.rand.glmr.nz.resid <- residuals(tds.rand.glmr.nz)
tds.rand.glmr.nz.resid.plot <- qqplot.data(tds.rand.glmr.nz.resid)
tds.rand.glmr.nz.resid.resids <- data.frame(resid = tds.rand.glmr.nz.resid)
tds.rand.glmr.nz.resid.dist.plot <- ggplot(tds.rand.glmr.nz.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)

#AICtab(tds.rand.zinb, tds.rand.glmr, tds.rand.glmr.nz)

# binomial glm of zero vs. non-zero cases #----
tds.rand.bglm.nz <- glmmTMB(tds_mph.nz ~
                              tchiyr.std +
                              hsz.std + nsk.std +
                              I(stthr.std^2) +
                              tchiyr.std:nsk.std +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
                              I(stthr.std^2):hsz.std +
                              I(stthr.std^2):nsk.std +
                              (1+I(stthr.std^2)|aclew_child_id),
                            data=quantity.rand,
                            family = "binomial")
#summary(tds.rand.bglm.nz)
#hsz.std                    -2.1306     0.9278  -2.296   0.0217 *
#nsk.std                     2.1014     0.8366   2.512   0.0120 *


####                    TDS-SA                    ####
# gaussian glm of log transformed DV (w/ zeroes) #----
tds.rand.sa.glmr <- glmmTMB(log(tds_mph.sa+1) ~
                              tchiyr.std + mated.bi.std + SpkrAge + 
                              hsz.std + nsk.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bi.std +
                              tchiyr.std:nsk.std +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
                              SpkrAge:tchiyr.std +
                              SpkrAge:nsk.std +
                              (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                            data=quantity.rand.sa)
#summary(tds.rand.sa.glmr)
#tchiyr.std                 0.27459    0.12558   2.187  0.02877 *  
#SpkrAgeChild              -0.58599    0.10904  -5.374 7.69e-08 ***
#nsk.std                    0.15064    0.08105   1.859  0.06309 .  
#tchiyr.std:nsk.std         0.16156    0.07618   2.121  0.03395 *  
#tchiyr.std:I(stthr.std^2) -2.22645    0.79476  -2.801  0.00509 ** 
tds.rand.sa.glmr.resid <- residuals(tds.rand.sa.glmr)
tds.rand.sa.glmr.resid.plot <- qqplot.data(tds.rand.sa.glmr.resid)
tds.rand.sa.glmr.resid.resids <- data.frame(resid = tds.rand.sa.glmr.resid)
tds.rand.sa.glmr.resid.dist.plot <- ggplot(tds.rand.sa.glmr.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)
# very skewed residuals

# gaussian glm of log transformed DV (w/o zeroes) #----
tds.rand.sa.glmr.nz <- glmmTMB(log(tds_mph.sa+1) ~
                                 tchiyr.std + mated.bi.std + SpkrAge + 
                                 hsz.std + nsk.std +
                                 I(stthr.std^2) +
                                 tchiyr.std:mated.bi.std +
                                 tchiyr.std:nsk.std +
                                 tchiyr.std:hsz.std +
                                 tchiyr.std:I(stthr.std^2) +
                                 SpkrAge:tchiyr.std +
                                 SpkrAge:nsk.std +
                                 (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                               data=subset(quantity.rand.sa, tds_mph.sa > 0))
#summary(tds.rand.sa.glmr.nz)
#SpkrAgeChild              -0.47411    0.24805  -1.911   0.0560 .  
#nsk.std                   -0.25804    0.11906  -2.167   0.0302 *
tds.rand.sa.glmr.nz.resid <- residuals(tds.rand.sa.glmr.nz)
tds.rand.sa.glmr.nz.resid.plot <- qqplot.data(tds.rand.sa.glmr.nz.resid)
tds.rand.sa.glmr.nz.resid.resids <- data.frame(resid = tds.rand.sa.glmr.nz.resid)
tds.rand.sa.glmr.nz.resid.dist.plot <- ggplot(tds.rand.sa.glmr.nz.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)

# binomial glm of zero vs. non-zero cases #----
tds.rand.sa.bglm.nz <- glmmTMB(tds_mph.sa.nz ~
                                 tchiyr.std + SpkrAge +
                                 hsz.std + nsk.std +
                                 I(stthr.std^2) +
                                 tchiyr.std:nsk.std +
                                 tchiyr.std:hsz.std +
                                 tchiyr.std:I(stthr.std^2) +
                                 SpkrAge:tchiyr.std +
                                 SpkrAge:nsk.std +
                                 (1+I(stthr.std^2)|aclew_child_id),
                               data=quantity.rand.sa,
                               family = "binomial")
#summary(tds.rand.sa.bglm.nz)
#hsz.std                    -2.1306     0.9278  -2.296   0.0217 *
#nsk.std                     2.1014     0.8366   2.512   0.0120 *


####                    PRP-TDS                    ####
# gaussian glm of log transformed DV (w/ zeroes) #----
prptds.rand.glmr <- glmmTMB(log(prop_tds+1) ~
                              tchiyr.std + mated.bi.std + 
                              hsz.std + nsk.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bi.std +
                              tchiyr.std:nsk.std +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
                              (1+I(stthr.std^2)|aclew_child_id),
                            data=quantity.rand.prpdata)
#summary(prptds.rand.glmr)
#nsk.std                   -0.166023   0.034005  -4.882 1.05e-06 ***
#tchiyr.std:I(stthr.std^2) -0.798645   0.320035  -2.495   0.0126 *
prptds.rand.glmr.resid <- residuals(prptds.rand.glmr)
prptds.rand.glmr.resid.plot <- qqplot.data(prptds.rand.glmr.resid)
prptds.rand.glmr.resid.resids <- data.frame(resid = prptds.rand.glmr.resid)
prptds.rand.glmr.resid.dist.plot <- ggplot(prptds.rand.glmr.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)


####                    N-SKRS                    ####
# poisson regression #----
nspkrs.rand.pois <- glmmTMB(n_spkrs_clip ~
                              tchiyr.std + mated.bin +
                              hsz.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bin +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
                              (1|aclew_child_id),
                            data=quantity.rand,
                            family="poisson")
#summary(nspkrs.rand.pois)
#mated.bin6+                0.41524    0.22444   1.850   0.0643 .  
#I(stthr.std^2)             1.74301    0.77140   2.260   0.0238 *  

# gaussian glm of log transformed DV (w/ zeroes) #----
nspkrs.rand.glmr <- glmmTMB(log(n_spkrs_clip+1) ~
                              tchiyr.std + mated.bi.std +
                              hsz.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bi.std +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
                              (1+I(stthr.std^2)|aclew_child_id),
                            data=quantity.rand)
#summary(nspkrs.rand.glmr)
#mated.bi.std6+             0.71653    0.25329   2.829  0.00467 **
nspkrs.rand.glmr.resid <- residuals(nspkrs.rand.glmr)
nspkrs.rand.glmr.resid.plot <- qqplot.data(nspkrs.rand.glmr.resid)
nspkrs.rand.glmr.resid.resids <- data.frame(resid = nspkrs.rand.glmr.resid)
nspkrs.rand.glmr.resid.dist.plot <- ggplot(nspkrs.rand.glmr.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)
# A bit skewed
################################################################################










################################################################################
####          TT SUPPLEMENTARY MODELS          ####
####                    ODS                    ####
# gaussian glm of log transformed DV (w/ zeroes) #----
ods.nonrand.tt.glmr <- glmmTMB(log(ods_mph+1) ~
                           tchiyr.std + mated.bi.std +
                           hsz.std + nsk.std +
                           I(stthr.std^2) +
                           tchiyr.std:mated.bi.std +
                           tchiyr.std:nsk.std +
                           tchiyr.std:hsz.std +
                           tchiyr.std:I(stthr.std^2) +
#                           offset(segment_dur) + # weird outcomes
                           (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                         data=quantity.nonrand.tt)
#summary(ods.nonrand.tt.glmr)
#hsz.std                   -0.31719    0.17229  -1.841   0.0656 .  
#nsk.std                    1.03206    0.12868   8.020 1.05e-15 ***
#I(stthr.std^2)            -3.24741    1.67292  -1.941   0.0522 .
ods.nonrand.tt.glmr.resid <- residuals(ods.nonrand.tt.glmr)
ods.nonrand.tt.glmr.resid.plot <- qqplot.data(ods.nonrand.tt.glmr.resid)
ods.nonrand.tt.glmr.resid.resids <- data.frame(resid = ods.nonrand.tt.glmr.resid)
ods.nonrand.tt.glmr.resid.dist.plot <- ggplot(ods.nonrand.tt.glmr.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)

# gaussian glm of log transformed DV (w/o zeroes) #----
ods.nonrand.tt.glmr.nz <- glmmTMB(log(ods_mph+1) ~
                              tchiyr.std + mated.bi.std +
                              hsz.std + nsk.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bi.std +
                              tchiyr.std:nsk.std +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
#                              offset(segment_dur) + # weird outcomes
                              (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                            data=subset(quantity.nonrand.tt, ods_mph > 0))
#summary(ods.nonrand.tt.glmr.nz)
#tchiyr.std                -0.40680    0.19990  -2.035   0.0418 *  
#hsz.std                   -0.36809    0.16518  -2.228   0.0259 *  
#nsk.std                    0.88867    0.12548   7.082 1.42e-12 ***
#I(stthr.std^2)            -3.66676    1.95356  -1.877   0.0605 .  
ods.nonrand.tt.glmr.nz.resid <- residuals(ods.nonrand.tt.glmr.nz)
ods.nonrand.tt.glmr.nz.resid.plot <- qqplot.data(ods.nonrand.tt.glmr.nz.resid)
ods.nonrand.tt.glmr.nz.resid.resids <- data.frame(resid = ods.nonrand.tt.glmr.nz.resid)
ods.nonrand.tt.glmr.nz.resid.dist.plot <- ggplot(ods.nonrand.tt.glmr.nz.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)

#AICtab(ods.nonrand.tt.zinb, ods.nonrand.tt.glmr, ods.nonrand.tt.glmr.nz)

# binomial glm of zero vs. non-zero cases #----
ods.nonrand.tt.bglm.nz <- glmmTMB(ods_mph.nz ~
                              tchiyr.std +
                              hsz.std + nsk.std +
                              I(stthr.std^2) +
#                              tchiyr.std:nsk.std + # no cnvgnce with this term
#                              tchiyr.std:hsz.std + # no cnvgnce with this term
                              tchiyr.std:I(stthr.std^2) +
#                              I(stthr.std^2):hsz.std + # no cnvgnce with this term
#                              I(stthr.std^2):nsk.std + # no cnvgnce with this term
#                              offset(segment_dur) + # no cnvgnce with this term
                              (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                            data=quantity.nonrand.tt,
                            family = "binomial")
#summary(ods.nonrand.tt.bglm.nz)
# no significant effects


####                    TDS                    ####
# gaussian glm of log transformed DV (w/ zeroes) #----
tds.nonrand.tt.glmr <- glmmTMB(log(tds_mph+1) ~
                           tchiyr.std + mated.bi.std +
                           hsz.std + nsk.std +
                           I(stthr.std^2) +
                           tchiyr.std:mated.bi.std +
                           tchiyr.std:nsk.std +
                           tchiyr.std:hsz.std +
                           tchiyr.std:I(stthr.std^2) +
#                           offset(segment_dur) +
                           (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                         data=quantity.nonrand.tt)
#summary(tds.nonrand.tt.glmr)
#mated.bi.std6+             0.50615    0.20436   2.477   0.0133 *  
#nsk.std                   -0.23926    0.11965  -2.000   0.0455 *
tds.nonrand.tt.glmr.resid <- residuals(tds.nonrand.tt.glmr)
tds.nonrand.tt.glmr.resid.plot <- qqplot.data(tds.nonrand.tt.glmr.resid)
tds.nonrand.tt.glmr.resid.resids <- data.frame(resid = tds.nonrand.tt.glmr.resid)
tds.nonrand.tt.glmr.resid.dist.plot <- ggplot(tds.nonrand.tt.glmr.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)

# gaussian glm of log transformed DV (w/o zeroes) #----
tds.nonrand.tt.glmr.nz <- glmmTMB(log(tds_mph+1) ~
                              tchiyr.std + mated.bi.std +
                              hsz.std + nsk.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bi.std +
                              tchiyr.std:nsk.std +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
#                              offset(segment_dur) +
                             (1|aclew_child_id),
                            data=subset(quantity.nonrand.tt, tds_mph > 0))
#summary(tds.nonrand.tt.glmr.nz)
#mated.bi.std6+             0.40435    0.19524   2.071   0.0384 *  
#nsk.std                   -0.20225    0.11321  -1.787   0.0740 .  
#tchiyr.std:hsz.std        -0.34299    0.18846  -1.820   0.0688 .  
tds.nonrand.tt.glmr.nz.resid <- residuals(tds.nonrand.tt.glmr.nz)
tds.nonrand.tt.glmr.nz.resid.plot <- qqplot.data(tds.nonrand.tt.glmr.nz.resid)
tds.nonrand.tt.glmr.nz.resid.resids <- data.frame(resid = tds.nonrand.tt.glmr.nz.resid)
tds.nonrand.tt.glmr.nzresid.dist.plot <- ggplot(tds.nonrand.tt.glmr.nz.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)
# somewhat skewed residuals

#AICtab(tds.nonrand.tt.zinb, tds.nonrand.tt.glmr, tds.nonrand.tt.glmr.nz)

# binomial glm of zero vs. non-zero cases #----
# not sensible: only one non-zero case


####                    TDS-SA                    ####
# gaussian glm of log transformed DV (w/ zeroes) #---- tds.nonrand.tt.sa.zinb
tds.nonrand.tt.sa.glmr <- glmmTMB(log(tds_mph.sa+1) ~
                              tchiyr.std + mated.bi.std + SpkrAge + 
                              hsz.std + nsk.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bi.std +
                              tchiyr.std:nsk.std +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
                              SpkrAge:tchiyr.std +
                              SpkrAge:nsk.std +
#                              offset(segment_dur) +
                              (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                            data=quantity.nonrand.tt.sa)
#summary(tds.nonrand.tt.sa.glmr)
#SpkrAgeChild              -1.38508    0.18282  -7.576 3.56e-14 ***
#tchiyr.std:nsk.std         0.21601    0.11232   1.923   0.0545 .  
#tchiyr.std:SpkrAgeChild    0.38339    0.18280   2.097   0.0360 *  
tds.nonrand.tt.sa.glmr.resid <- residuals(tds.nonrand.tt.sa.glmr)
tds.nonrand.tt.sa.glmr.resid.plot <- qqplot.data(tds.nonrand.tt.sa.glmr.resid)
tds.nonrand.tt.sa.glmr.resid.resids <- data.frame(resid = tds.nonrand.tt.sa.glmr.resid)
tds.nonrand.tt.sa.glmr.resid.dist.plot <- ggplot(tds.nonrand.tt.sa.glmr.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)

# gaussian glm of log transformed DV (w/o zeroes) #----
tds.nonrand.tt.sa.glmr.nz <- glmmTMB(log(tds_mph.sa+1) ~
                                 tchiyr.std + mated.bi.std + SpkrAge + 
                                 hsz.std + nsk.std +
                                 I(stthr.std^2) +
                                 tchiyr.std:mated.bi.std +
                                 tchiyr.std:nsk.std +
                                 tchiyr.std:hsz.std +
                                 tchiyr.std:I(stthr.std^2) +
                                 SpkrAge:tchiyr.std +
                                 SpkrAge:nsk.std +
#                              offset(segment_dur) +
                                 (1|aclew_child_id), # max: no cnvgnce w/ stthr slopes
                               data=subset(quantity.nonrand.tt.sa, tds_mph.sa > 0))
#summary(tds.nonrand.tt.sa.glmr.nz)
#tchiyr.std                -0.31820    0.16694  -1.906   0.0566 .  
#SpkrAgeChild              -0.94293    0.21079  -4.473 7.71e-06 ***
#nsk.std                   -0.21429    0.12947  -1.655   0.0979 .  
#tchiyr.std:hsz.std        -0.36036    0.15337  -2.350   0.0188 *  
#tchiyr.std:SpkrAgeChild    0.92202    0.20188   4.567 4.95e-06 ***
#SpkrAgeChild:nsk.std       0.48064    0.25005   1.922   0.0546 .  
tds.nonrand.tt.sa.glmr.nz.resid <- residuals(tds.nonrand.tt.sa.glmr.nz)
tds.nonrand.tt.sa.glmr.nz.resid.plot <- qqplot.data(tds.nonrand.tt.sa.glmr.nz.resid)
tds.nonrand.tt.sa.glmr.nz.resid.resids <- data.frame(resid = tds.nonrand.tt.sa.glmr.nz.resid)
tds.nonrand.tt.sa.glmr.nz.resid.dist.plot <- ggplot(tds.nonrand.tt.sa.glmr.nz.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)

# binomial glm of zero vs. non-zero cases #----
tds.nonrand.tt.sa.bglm.nz <- glmmTMB(tds_mph.sa.nz ~
                                 tchiyr.std + SpkrAge +
                                 hsz.std + nsk.std +
                                 I(stthr.std^2) +
                                 tchiyr.std:nsk.std +
                                 tchiyr.std:hsz.std +
                                 tchiyr.std:I(stthr.std^2) +
                                 SpkrAge:tchiyr.std +
                                 SpkrAge:nsk.std +
#                                 offset(segment_dur) +
                                 (1|aclew_child_id),
                               data=quantity.nonrand.tt.sa,
                               family = "binomial")
#summary(tds.nonrand.tt.sa.bglm.nz)
#tchiyr.std                  1.2604     0.6778   1.860 0.062945 .  
#SpkrAgeChild               -2.4922     0.5916  -4.213 2.52e-05 ***
#hsz.std                     1.0093     0.5289   1.908 0.056354 .  
#nsk.std                     2.3990     0.9539   2.515 0.011904 *  
#tchiyr.std:nsk.std          1.1208     0.6343   1.767 0.077198 .  


####                    PRP-TDS                    ####
# gaussian glm of log transformed DV (w/ zeroes) #----
prptds.nonrand.tt.glmr <- glmmTMB(log(prop_tds+1) ~
                              tchiyr.std + mated.bi.std + 
                              hsz.std + nsk.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bi.std +
                              tchiyr.std:nsk.std +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
                              offset(segment_dur) +
                              (1|aclew_child_id),
                            data=quantity.nonrand.tt.prpdata)
#summary(prptds.nonrand.tt.glmr)
#nsk.std                   -0.70183    0.23866  -2.941  0.00327 ** 
#tchiyr.std:nsk.std        -0.54074    0.28578  -1.892  0.05847 . 
prptds.nonrand.tt.glmr.resid <- residuals(prptds.nonrand.tt.glmr)
prptds.nonrand.tt.glmr.resid.plot <- qqplot.data(prptds.nonrand.tt.glmr.resid)
prptds.nonrand.tt.glmr.resid.resids <- data.frame(resid = prptds.nonrand.tt.glmr.resid)
prptds.nonrand.tt.glmr.resid.dist.plot <- ggplot(prptds.nonrand.tt.glmr.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)
# very bad distribution: bimodal


####                    N-SKRS                    ####
# poisson regression #----
nspkrs.nonrand.tt.pois <- glmmTMB(n_spkrs_clip ~
                              tchiyr.std + mated.bin +
                              hsz.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bin +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
                              (1|aclew_child_id),
                            data=quantity.nonrand.tt,
                            family="poisson")
#summary(nspkrs.nonrand.tt.pois)
#mated.bin6+                0.48389    0.17399   2.781 0.005417 ** 
#tchiyr.std:I(stthr.std^2) -2.49781    1.44738  -1.726 0.084394 . 

# gaussian glm of log transformed DV (w/ zeroes) #----
nspkrs.nonrand.tt.glmr <- glmmTMB(log(n_spkrs_clip+1) ~
                              tchiyr.std + mated.bi.std +
                              hsz.std +
                              I(stthr.std^2) +
                              tchiyr.std:mated.bi.std +
                              tchiyr.std:hsz.std +
                              tchiyr.std:I(stthr.std^2) +
                              offset(segment_dur) +
                              (1|aclew_child_id),
                            data=quantity.nonrand.tt)
#summary(nspkrs.nonrand.tt.glmr)
#I(stthr.std^2)             5.36339    3.03256   1.769   0.0770 .
nspkrs.nonrand.tt.glmr.resid <- residuals(nspkrs.nonrand.tt.glmr)
nspkrs.nonrand.tt.glmr.resid.plot <- qqplot.data(nspkrs.nonrand.tt.glmr.resid)
nspkrs.nonrand.tt.glmr.resid.resids <- data.frame(resid = nspkrs.nonrand.tt.glmr.resid)
nspkrs.nonrand.tt.glmr.resid.dist.plot <- ggplot(nspkrs.nonrand.tt.glmr.resid.resids, aes(x=resid)) +
                              geom_histogram(aes(y=..density..),
                                             fill="white", color="black",
                                             binwidth = 0.5) +
                             geom_density(color="firebrick3", lwd=2)
# very bad distribution: bimodal
################################################################################
