## mad_saccade
**Robust saccade threshold estimation**

Code for *Voloh, Watson, Konig, and Womelsdorf. (2020) MAD saccade: statistically robust saccade threshold estimation via the median absolute deviation. Journal of Eye Movement Research*

## Function Description
__Overview__: Artificial scanpath are generated by **Generate_Artificial_Scan_Paths** (Dai et al 2016), which calls **fun_sacc_model** and **sacc_model_type1** to do so. Saccades are then detected via the adaptive velocity threshold algorithm (Nyström & Holmqvist 2010). The performance of the algorithm is evaluated via **EventConfusion_Expanded** (Warby et al 2014).

* Generate_Artificial_Scan_Paths.m
  * generates *num_simulations* scanpaths for each of *noise_levels* with a total of *NSAC* saccades. *min_fix_dur*, *min_sac_amp*, *max_sac_amp* define the characteristics of the saccade. Simullations are saved in the folder *savepath*
* sacc_model_type1.m	
 * [ model, model_amp, model_pvel, model_vel ] = sacc_model_type1( Stype, t, p )
 * *Stype* is the model type (see function for details). "ms" is the model as proposed by Dai et al 2016. *t* is the time, and *p* is a vector of model parameters eg [eta, c, tau, alpha, beta];
* fun_sacc_model.m	
 * [ model, model_amp, model_pvel, model_vel ] = fun_sacc_model( t, p ) 
 * *t* is the time. *p* are the model parameters: [eta, c, tau, alpha, beta]
* adaptive_VT_mad.m	
 * [fixationtimes,saccadetimes,saccadeIdx,info] = adaptive_VT_mad(eyedat,samplingFreq,useMAD,exciseIntersaccade,initVel,lambda,fixedThreshVal)
 *__inputs__: *eyedat* is an Nx2 matrix of x- and y-coordinates. *samplingFreq* is the sampling frequency of eyedat. *useMAD*, if true, uses robust statsutics to estimate the mean and standard deviation for threshold etsimation. *exciseIntersaccade*, if true, excises some data points at the start and end of the inter-saccadic intervals when calculating the threshold. *initVel* is the initial velocity used for threshold calculation. *lambda* is the scaling parameter used to calculate the threshold (=6 in Nyström & Holmqvist 2010). If *fixedThreshVal* is set to a value, a fixed threshold is used instead to determine saccade. In this case, the peakDetectionThreshold=fixedThreshVal, and saccadeVelocityTreshold=fixedThreshVal-10
 *__outputs__: *fixationtimes* is an 2xN matrix of N fixation start and end times. *saccadetimes* is the same for saccades. *saccadeIdx* is a vector of the same length as the data, with 1s if samples belong to a saccade and 0 otehrwise. *info* is a sturtcure containing information about the algorithm performance
* EventConfusion_Expanded.m	
 * [ F1score, recall, precision, TP, FP, FN,missed_saccades,onset_lag,offset_lag,extra_saccades] =...
    EventConfusion_Expanded( reference, predicted )
    * calculates performance metrics comparing a time series ground-truth *reference* and algorithmic *predicted* time-series of 1s and 0s. For example, *predicted* is the *saccadeIdx* from adaptive_VT_mad, and *reference* is the known saccade times from Generate_Artificial_Scan_Paths



## References and prior work

Nyström, M., & Holmqvist, K. (2010). An adaptive algorithm for fixation, saccade, and glissade detection in eyetracking data. Behavior Research Methods, 42(1), 188–204. https://doi.org/10.3758/BRM.42.1.188

Dai, W., Selesnick, I., Rizzo, J.-R., Rucker, J., & Hud-son, T. (2016). A parametric model for saccadic eye movement. 2016 IEEE Signal Processing in Medicine and Biology Symposium (SPMB), 1–6. https://doi.org/10.1109/SPMB.2016.7846860

Warby, S. C., Wendt, S. L., Welinder, P., Munk, E. G. S., Carrillo, O., Sorensen, H. B. D., Jennum, P., Peppard, P. E., Perona, P., & Mignot, E. (2014). Sleep-spindle detection: Crowdsourcing and evaluat-ing performance of experts, non-experts and automat-ed methods. Nature Methods, 11(4), 385–392. https://doi.org/10.1038/nmeth.2855

Watson, M. R., Voloh, B., Thomas, C. J., & Wom-elsdorf, T. (2019). USE: An integrative suite for temporally-precise psychophysical experiments in virtual environments. Journal of Neuroscience Methods, 326. https://doi.org/10.1016/j.jneumeth.2019.108374
