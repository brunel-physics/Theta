from ROOT import gROOT, TCanvas, TF1, TH1F, TFile, TGraph
from array import array

import numpy as n
 
dosysttable  = True
dobiasscan   = False
dologprofile = False
dodiscovery  = True
dolimits     = True
dolimits_bayes   = True
dolimits_cls     = False

options = Options()
options.set('minimizer', 'strategy', 'robust')

optionscls = Options()
optionscls.set('minimizer', 'strategy', 'robust')
  
# -------------- TO CHANGE BY THE USER

signalname = 'tZq'
signal_init_xs =  0.008184
# Need to update with di-lepton x-section
#signal_init_xs = 0.00783  
# -------------- TO CHANGE BY THE USER
# for model building:
def get_model(signalname):

# Read in and build the model automatically from the histograms in the root file. This model will contain all shape uncertainties given according to the templates-which also includes rate changes according to the alternate shapes. For more info about this model and naming conventuion, see documentation of build_model_from_rootfile.

    model = build_model_from_rootfile('MVA_merged_theta.root',  include_mc_uncertainties=True)
    
# If the prediction histogram is zero, but data is non-zero, the negative log-likelihood is infinite which causes problems for some methods. Therefore, we set all histogram bin entries to a small, but positive value:

    model.fill_histogram_zerobins()

# define what the signal processes are. All other processes are assumed to make up the 'background-only' model.

    model.set_signal_processes(signalname)
        
    model.add_lognormal_uncertainty('ZZ_rate',     math.log(1.3), 'ZZ')
# model.add_lognormal_uncertainty('tZq_rate',    math.log(1.1), 'tZq')

# WZsplit -> WZbc and WZl 19/11/2015 
# For di-lepton change abck to WZ   
# and then commented out as this is zero for di-lepton....needs re-checking

#    model.add_lognormal_uncertainty('WZbc_rate',   math.log(2.0), 'WZbc')
#    model.add_lognormal_uncertainty('WZl_rate',     math.log(2.0), 'WZl')
#    model.add_lognormal_uncertainty('WZ_rate',     math.log(2.0), 'WZ')
    
    uncertainty_mu = 5.0
    uncertainty_el = 5.0

# Commented out for di-lepton...needs checking
     
    #for the mumu channel
#    model.add_lognormal_uncertainty('FakeRate_mu', math.log(uncertainty_mu),       'Zjets', 'MVA_mumu')
    #for the mumue channel
#    model.add_lognormal_uncertainty('FakeRate_mu', math.log(uncertainty_mu*2./3.), 'Zjets', 'MVA_emumu')
#    model.add_lognormal_uncertainty('FakeRate_el', math.log(uncertainty_el/3.),    'Zjets', 'MVA_emumu')
    #for the eemuchannel
#    model.add_lognormal_uncertainty('FakeRate_mu', math.log(uncertainty_mu/3.), #   'Zjets', 'MVA_eemu')
#    model.add_lognormal_uncertainty('FakeRate_el', math.log(uncertainty_el*2./3.), 'Zjets', 'MVA_eemu')
    #for the eee channel
#    model.add_lognormal_uncertainty('FakeRate_el', math.log(uncertainty_el), 'Zjets', 'MVA_ee') 
#    model.add_lognormal_uncertainty('FakeRate_mu', math.log(uncertainty_mu),       'Zjets', 'MVA_mumumu')

# CKM - mTW sig
    
    #for the mumumu channel
#    model.add_lognormal_uncertainty('FakeRate_mu', math.log(uncertainty_mu),	   'Zjets', 'mTW_mumumu')
    #for the mumue channel
#    model.add_lognormal_uncertainty('FakeRate_mu', math.log(uncertainty_mu*2./3.), 'Zjets', 'mTW_emumu')
#    model.add_lognormal_uncertainty('FakeRate_el', math.log(uncertainty_el/3.),    'Zjets', 'mTW_emumu')
    #for the eemuchannel
#    model.add_lognormal_uncertainty('FakeRate_mu', math.log(uncertainty_mu/3.),    'Zjets', 'mTW_eemu')
#    model.add_lognormal_uncertainty('FakeRate_el', math.log(uncertainty_el*2./3.), 'Zjets', 'mTW_eemu')
#for the eee channel
#    model.add_lognormal_uncertainty('FakeRate_el', math.log(uncertainty_el), 'Zjets', 'mTW_eee') 
       
#    model.add_lognormal_uncertainty('Zjets_rate',   math.log(1.3), 'Zjets')    
    model.add_lognormal_uncertainty('TTZ_rate',     math.log(1.3), 'TTZ')
    model.add_lognormal_uncertainty('TTW_rate',     math.log(1.3), 'TTW')
    
    #model.add_lognormal_uncertainty('tZq_rate',   math.log(1.1), signalname)
    
# Add some lognormal rate uncertainties. The first parameter is the name of the uncertainty (which will also be the name of the nuisance parameter), the second is the 'effect' as a fraction, the third one is the process name. The fourth parameter is optional and denotes the channl. The default '*' means that the uncertainty applies to all channels in the same way. Note that you can use the same name for a systematic here as for a shape systematic. In this case, the same parameter will be used; shape and rate changes will be 100% correlated.
        
    for p in model.processes:
    	if p == 'Zjets': continue
    	if p == 'TT'   : continue
    	model.add_lognormal_uncertainty('lumi',        math.log(1.026), p)
        #model.add_lognormal_uncertainty('TrigLept',    math.log(1.05), p)
       
    return model

model = get_model(signalname)


print ("------------------------------------------------------------------")
print ("------------------------fix systematics---------------------------")

# set list of systematics to be externalised
fixed_syst_list = ['pdf', 'matching','scale'] 

for fix_uncertainties in (None, ['pdf','matching','scale'] ): 
    
    print "\nuncertainties not fitted: ", fix_uncertainties
# fixed_dist is a Distribution instance which fixes the uncerainties which should not be fitted. It can also be None, if no uncertainties are to be fixed (i.e., all are fitted)
    if fix_uncertainties is None: fixed_dist = None
    else: fixed_dist = get_fixed_dist_at_values(dict([(u, 0.0) for u in fix_uncertainties]))
    
    print fix_uncertainties
    
    print ("------------------------------------------------------------------")
    print ("------------------------------------------------------------------")
        
model_summary(model)

signal_shapes = {signalname: [signalname]}  

print ("----------------------------------------------------------------------")
print ("--------------compute cross-section with fixed parameters-------------")

one_sigma   = 0.682
two_sigma   = 0.9545
three_sigma = 0.9973

print ("measurement of the cross-section")

res    = pl_interval(model, 'data', n=1, cls = [one_sigma], signal_process_groups = signal_shapes, options = options ,  nuisance_constraint = fixed_dist)
res_2s = pl_interval(model, 'data', n=1, cls = [two_sigma], signal_process_groups = signal_shapes, options = options ,  nuisance_constraint = fixed_dist)
res_3s = pl_interval(model, 'data', n=1, cls = [three_sigma], signal_process_groups = signal_shapes, options = options ,nuisance_constraint = fixed_dist)

#res    = pl_interval(model, 'data', n=1, cls = [one_sigma], signal_process_groups = signal_shapes, options = options )
#res_2s = pl_interval(model, 'data', n=1, cls = [two_sigma], signal_process_groups = signal_shapes, options = options )
#res_3s = pl_interval(model, 'data', n=1, cls = [three_sigma], signal_process_groups = signal_shapes, options = options )

print [ "%.5f" % res[signalname][0][0] ,    "%.5f" %res[signalname][one_sigma][0][0] ,      "%.5f" %res[signalname][one_sigma][0][1] ]
print [ "%.5f" % res_2s[signalname][0][0] , "%.5f" %res_2s[signalname][two_sigma][0][0] ,   "%.5f" %res_2s[signalname][two_sigma][0][1] ]
print [ "%.5f" % res_3s[signalname][0][0] , "%.5f" %res_3s[signalname][three_sigma][0][0] , "%.5f" %res_3s[signalname][three_sigma][0][1] ]

syst_down   = 0.00000001
syst_up     = 0.00000001

if (res[signalname][0][0]> 0.0000001):
	syst_down   = (res[signalname][0][0]            - res[signalname][one_sigma][0][0])/res[signalname][0][0]
	syst_up     = (res[signalname][one_sigma][0][1] - res[signalname][0][0])/res[signalname][0][0]

interval = (res[signalname][one_sigma][0][1] - res[signalname][one_sigma][0][0])/2

signal_fit  = signal_init_xs*res[signalname][0][0]
signal_down = signal_init_xs*res[signalname][one_sigma][0][0]
signal_up   = signal_init_xs*res[signalname][one_sigma][0][1] 

signal_down_2S = signal_init_xs*res_2s[signalname][two_sigma][0][0]
signal_up_2S   = signal_init_xs*res_2s[signalname][two_sigma][0][1] 

signal_down_3S = signal_init_xs*res_3s[signalname][three_sigma][0][0]
signal_up_3S   = signal_init_xs*res_3s[signalname][three_sigma][0][1] 

print ["fitted cross section ", "%.5f" %signal_fit]
print ["down variation       ", "%.5f" %signal_down]
print ["up variation         ", "%.5f" %signal_up]

print ["cross section ", "%.5f" %signal_fit, "+%.5f" %(signal_up-signal_fit), "-%.5f" %(signal_fit - signal_down)]

filename = "tables.tex"
texfile = open(filename,'w')

print ("----------------------------------------------------------------------")
print ("-------------------account for externalized systematics---------------")

texfile.write( "the fitted beta signal   %.5f  [%.5f %.5f] \n" % ( res[signalname][0][0], res[signalname][one_sigma][0][0] , res[signalname][one_sigma][0][1] ) )
texfile.write( "the fitted cross-section %.5f  [%.5f %.5f] \n" % ( signal_fit, signal_down, signal_up) )
texfile.write( "Summary of externalize uncertainties \n")
texfile.write( "\\begin{table} \n")
texfile.write( "\\begin{center} \n")
texfile.write( "\\begin{tabular}{ |c|c|c| } \n")
texfile.write( "\\hline \n")
texfile.write( "Uncertainty source & up (\\%) & down (\\%) \\\\ \n")
texfile.write( "\\hline \n")

total_error_up   = (signal_up - signal_fit)**2
total_error_down = (signal_fit - signal_down)**2

fixed_uncert_error_up = 0
fixed_uncert_error_down = 0

for p in fix_uncertainties:

	#print [p, math.sqrt(total_error_up), math.sqrt(total_error_down)]
	model_tmp = model.copy()
	p_mean = model_tmp.distribution.get_distribution(p)['mean']
	p_width = model_tmp.distribution.get_distribution(p)['width']
	#-------------------------
	#for up shift
	#-------------------------
	nuisance_prior_toys = get_fixed_dist_at_values({p: p_mean + p_width})	
	res_up = mle(model_tmp, 'toys-asimov:1.23187', 1, nuisance_prior_toys = nuisance_prior_toys, nuisance_constraint = fixed_dist, options=options)
	beta_signal_fitted_up = res_up[signalname]['beta_signal'][0][0]
        shift_plus = beta_signal_fitted_up - res[signalname][0][0]
	xs_shift_plus= signal_init_xs*shift_plus
	
	
	#-------------------------
	#for down shift
	#-------------------------
	nuisance_prior_toys = get_fixed_dist_at_values({p: p_mean - p_width})
	res_down = mle(model, 'toys-asimov:1.23187', 1, nuisance_prior_toys = nuisance_prior_toys, nuisance_constraint = fixed_dist, options=options)
	beta_signal_fitted_down = res_down[signalname]['beta_signal'][0][0]
        shift_minus    = beta_signal_fitted_down - res[signalname][0][0] 
	xs_shift_minus = signal_init_xs*shift_minus
	
	print [p, beta_signal_fitted_up, beta_signal_fitted_down]
	print [p, shift_plus, shift_minus]
	
        texfile.write( "%s & %.5f  & %.5f \\\\ \n" % (p, shift_plus*100, shift_minus*100) )
	total_error_up   += xs_shift_plus**2
	total_error_down += xs_shift_minus**2
	fixed_uncert_error_up   += xs_shift_plus**2
	fixed_uncert_error_down += xs_shift_minus**2

print [ math.sqrt(total_error_up), math.sqrt(total_error_down)]	
texfile.write( "\\end{tabular}\n")
texfile.write( "\\end{center} \n")
texfile.write( "\\end{table} \n")

total_error_up   =  math.sqrt(total_error_up)
total_error_down =  math.sqrt(total_error_down)

print ["fitted cross section with externalized systematics", "%.5f" %signal_fit, "+ %.5f" %total_error_up, "- %.5f" %total_error_down]

texfile.write( "fitted cross section with externalized systematics %.4f+%.4f-%.4f pb \n" % (signal_fit, total_error_up, total_error_down) )
texfile.write( "up variation   %.5f \n" % total_error_up)

print ("------------------------------------------------------------------")
print ("------------------------------------------------------------------")

### For max. Likelihood Fit results

print ["start fit mle"]
fit = mle(model, input = 'data', n = 1, signal_process_groups = signal_shapes, with_covariance=True, with_error=True, ks = True, chi2 = True, options = options, nuisance_constraint = fixed_dist)
print ["mle fit done"]
# the output is (fitted value, uncertainty)
# The first numbers in the brackets show how far we are from the nominal value (which is 0) after the fit. 
#A value of 1 would mean 1 sigma deviation. So we are below 1 sigma deviation. 
#The second numbers in the brackets illustrates the uncertainty on the fitted value, it should be below 1, 
#and a value close to 1 corresponds to "no sensitivity" on the systematic.

print ("Determine nuisance parameters and their uncertainties")
parameter_values = {}
parameter_uncert = {}
for p in model.get_parameters([]):
    parameter_values[p] = fit[signalname][p][0][0]
    parameter_uncert[p] = fit[signalname][p][0][1]
    
    print [p, "%.4f" %parameter_values[p], "%.4f" %parameter_uncert[p] ]

parameter_values['beta_signal'] =  res[signalname][0][0]

histos = evaluate_prediction(model, parameter_values, include_signal = True)
write_histograms_to_rootfile(histos, 'results/histos-mle_mwt_merged.root')

print ("------------------------------------------------------------------")
print ("------------------------------------------------------------------")

print ("")
parameter_values_sfup = {}
parameter_uncert_sfup = {}
for q in model.get_parameters([]):
	for p in model.get_parameters([]):
    		parameter_values_sfup[p] = fit[signalname][p][0][0]
    		parameter_uncert_sfup[p] = fit[signalname][p][0][1]
    		if p == q : parameter_values_sfup[p] = fit[signalname][p][0][0]+fit[signalname][p][0][1]
	parameter_values_sfup['beta_signal'] =  res[signalname][0][0]
	histos = evaluate_prediction(model, parameter_values_sfup, include_signal = True)
	write_histograms_to_rootfile(histos, 'results/histos-mle_mwt_'+q+'.root')

if dosysttable:
	texfile.write( " \n")

	texfile.write( " \n")
	texfile.write( " \n")

	texfile.write( "\\begin{table} \n")
	texfile.write( "\\begin{center} \n")
	texfile.write( "\\begin{tabular}{ |c|c|c| } \n")
	texfile.write( "\\hline \n")
	texfile.write( "Uncertainty source & pb & \\%  \\\\ \n")
	texfile.write( "\\hline \n")
	
	print ("--------------------------------------------------------------")
	print ("-----------------------do systematic table--------------------")
	syst = 0
	tot_uncert =0

        print ["--------------------------------"]
	print ("Determine the impact of each systematic")
	for p in model.get_parameters([]):
		excluded_syst=False
		for q in fixed_syst_list:
			if(p==q): excluded_syst=True
		if(excluded_syst==True): continue
		model_syst = model.copy()
		model_syst.distribution.set_distribution_parameters(p, width = 0.0, mean = parameter_values[p], range = [parameter_values[p], parameter_values[p]])
		res_syst = pl_interval(model_syst, 'data', n=1, cls = [one_sigma], signal_process_groups = signal_shapes, options = options  )
	
       		
		interval_syst = (res_syst[signalname][one_sigma][0][1] - res_syst[signalname][one_sigma][0][0])/2
		syst  = (abs(interval**2 - interval_syst**2))**(0.5)
		tot_uncert = tot_uncert + (syst)**2
		print ["syst effect of (%)  ", p, "%.2f" %(syst*100)]
		print ["syst effect of (pb) ", p, "%.2f" %(syst*signal_fit)]
		print [p, "%.4f" %(syst*signal_fit), "%.4f" %(syst*100)]
		
		texfile.write( "%s & %.5f  & %.5f \\\\ \n" % (p, syst*signal_fit, syst*100) )
		
        	#print ["--------------------------------"]

	#print ["total syst down/up" ,"%.4f" % total_down**(0.5), "%.4f" %total_up**(0.5)]
	tot_uncert = tot_uncert**(0.5)
	print ["total uncert %" , tot_uncert*100]
	#tot_uncert = tot_uncert*signal_fit
	print ["total uncert pb" , tot_uncert*signal_fit]
        
	texfile.write( "\\hline \n")
	texfile.write( "total syst & %.5f  & %.5f \\\\ \n" % (tot_uncert*signal_fit, tot_uncert*100 ) )
	texfile.write( "\\hline \n")

	texfile.write( "\\end{tabular}\n")
	texfile.write( "\\end{center} \n")
	texfile.write( "\\end{table} \n")
	texfile.write( " \n")
	texfile.write( " \n")
	syst_error = ( (tot_uncert*signal_fit)**(2) + fixed_uncert_error_up**(2))**(0.5)
	texfile.write( " %.5f \n" % total_error_up)
	texfile.write( " %.5f \n" % syst_error)
        print ["total_error_up" , total_error_up]
        print  ["syst_error" , syst_error]
	#texfile.write( " %.5f \n" % (total_error_up**(2) - syst_error**(2))**(0.5) )
	stat_error = 0
	texfile.write( "final cross section is %.2f \n" % signal_fit )
	texfile.write( "total uncertainty %.1f %.2f\n" % ( total_error_up, total_error_down ) )
	#texfile.write( "stat %.2f \n" % (total_error_up**(2) - syst_error**(2))**(0.5) )
	texfile.write( "syst %.2f \n" % (syst_error) )


print ("----------------------------------------------------------------------")
print ("-----------------------do bias scan and pull--------------------------")

################################
#### Perform toy MC
################################

theBias0p5 = TH1F('theBias0p5', 'theBias0p5', 100, 0.3, 0.7)
theBias0p6 = TH1F('theBias0p6', 'theBias0p6', 100, 0.4, 0.8)
theBias0p7 = TH1F('theBias0p7', 'theBias0p7', 100, 0.5, 0.9)
theBias0p8 = TH1F('theBias0p8', 'theBias0p8', 100, 0.6, 1.0)
theBias0p9 = TH1F('theBias0p9', 'theBias0p9', 100, 0.7, 1.1)
theBias1p0 = TH1F('theBias1p0', 'theBias1p0', 100, 0.8, 1.2)
theBias1p1 = TH1F('theBias1p1', 'theBias1p1', 100, 0.9, 1.3)
theBias1p2 = TH1F('theBias1p2', 'theBias1p2', 100, 1.0, 1.4)
theBias1p3 = TH1F('theBias1p3', 'theBias1p3', 100, 1.1, 1.5)
theBias1p4 = TH1F('theBias1p4', 'theBias1p4', 100, 1.2, 1.6)
theBias1p5 = TH1F('theBias1p5', 'theBias1p5', 100, 1.3, 1.7)


thePull0p5 = TH1F('thePull0p5', 'thePull0p5', 200, -4, 4)
thePull0p6 = TH1F('thePull0p6', 'thePull0p6', 200, -4, 4)
thePull0p7 = TH1F('thePull0p7', 'thePull0p7', 200, -4, 4)
thePull0p8 = TH1F('thePull0p8', 'thePull0p8', 200, -4, 4)
thePull0p9 = TH1F('thePull0p9', 'thePull0p9', 200, -4, 4)
thePull1p0 = TH1F('thePull1p0', 'thePull1p0', 200, -4, 4)
thePull1p1 = TH1F('thePull1p1', 'thePull1p1', 200, -4, 4)
thePull1p2 = TH1F('thePull1p2', 'thePull1p2', 200, -4, 4)
thePull1p3 = TH1F('thePull1p3', 'thePull1p3', 200, -4, 4)
thePull1p4 = TH1F('thePull1p4', 'thePull1p4', 200, -4, 4)
thePull1p5 = TH1F('thePull1p5', 'thePull1p5', 200, -4, 4)

if dobiasscan:
	print ("--------------------------------------------------------------")
	print ("-----------------------do bias scan and pull------------------")
	theBias0p5 = TH1F('theBias0p5', 'theBias0p5', 100, 0.3, 0.7)
	theBias0p6 = TH1F('theBias0p6', 'theBias0p6', 100, 0.4, 0.8)
	theBias0p7 = TH1F('theBias0p7', 'theBias0p7', 100, 0.5, 0.9)
	theBias0p8 = TH1F('theBias0p8', 'theBias0p8', 100, 0.6, 1.0)
	theBias0p9 = TH1F('theBias0p9', 'theBias0p9', 100, 0.7, 1.1)
	theBias1p0 = TH1F('theBias1p0', 'theBias1p0', 100, 0.8, 1.2)
	theBias1p1 = TH1F('theBias1p1', 'theBias1p1', 100, 0.9, 1.3)
	theBias1p2 = TH1F('theBias1p2', 'theBias1p2', 100, 1.0, 1.4)
	theBias1p3 = TH1F('theBias1p3', 'theBias1p3', 100, 1.1, 1.5)
	theBias1p4 = TH1F('theBias1p4', 'theBias1p4', 100, 1.2, 1.6)
	theBias1p5 = TH1F('theBias1p5', 'theBias1p5', 100, 1.3, 1.7)


	thePull0p5 = TH1F('thePull0p5', 'thePull0p5', 200, -4, 4)
	thePull0p6 = TH1F('thePull0p6', 'thePull0p6', 200, -4, 4)
	thePull0p7 = TH1F('thePull0p7', 'thePull0p7', 200, -4, 4)
	thePull0p8 = TH1F('thePull0p8', 'thePull0p8', 200, -4, 4)
	thePull0p9 = TH1F('thePull0p9', 'thePull0p9', 200, -4, 4)
	thePull1p0 = TH1F('thePull1p0', 'thePull1p0', 200, -4, 4)
	thePull1p1 = TH1F('thePull1p1', 'thePull1p1', 200, -4, 4)
	thePull1p2 = TH1F('thePull1p2', 'thePull1p2', 200, -4, 4)
	thePull1p3 = TH1F('thePull1p3', 'thePull1p3', 200, -4, 4)
	thePull1p4 = TH1F('thePull1p4', 'thePull1p4', 200, -4, 4)
	thePull1p5 = TH1F('thePull1p5', 'thePull1p5', 200, -4, 4)
	for i in range(11):
		print ("perform toy MC for bias scan")
		fixed_dist = get_fixed_dist(model.distribution)
		#mle(model, "toys:1.0", 1000, nuisance_prior_toys = fixed_dist)
		if i==0:
			res_toy = pl_interval(model, 'toys:0.5', n=10100, cls = [one_sigma], signal_process_groups = signal_shapes, nuisance_prior_toys = fixed_dist, options = options )
			#res_toy = pl_interval(model, 'toys:0.5', n=10100, cls = [one_sigma], signal_process_groups = signal_shapes , options = options )
			for j in range(10000):
				theBias0p5.Fill(res_toy[signalname][0][j])
				pull = (res_toy[signalname][0][j] - 0.5)/(0.5*(res_toy[signalname][one_sigma][j][1] -res_toy[signalname][one_sigma][j][0]))
				thePull0p5.Fill(pull)
			#print ["%.4f" % res_toy[signalname][0][0], res_toy[signalname][one_sigma][0][1], res_toy[signalname][one_sigma][0][1]]
	        if i==1:
		        res_toy = pl_interval(model, 'toys:0.6', n=10100, cls = [one_sigma], signal_process_groups = signal_shapes  , options = options)
		        for j in range(10000):
		       		theBias0p6.Fill(res_toy[signalname][0][j])
				pull = (res_toy[signalname][0][j] - 0.6)/(0.5*(res_toy[signalname][one_sigma][j][1] -res_toy[signalname][one_sigma][j][0]))
				thePull0p6.Fill(pull)
	        if i==2:
		        res_toy = pl_interval(model, 'toys:0.7', n=10100, cls = [one_sigma], signal_process_groups = signal_shapes , options = options )
		        for j in range(10000):
			        theBias0p7.Fill(res_toy[signalname][0][j])
				pull = (res_toy[signalname][0][j] - 0.7)/(0.5*(res_toy[signalname][one_sigma][j][1] -res_toy[signalname][one_sigma][j][0]))
				thePull0p7.Fill(pull)
	        if i==3:
		        res_toy = pl_interval(model, 'toys:0.8', n=10100, cls = [one_sigma], signal_process_groups = signal_shapes , options = options)
		        for j in range(10000):
			        theBias0p8.Fill(res_toy[signalname][0][j])
				pull = (res_toy[signalname][0][j] - 0.8)/(0.5*(res_toy[signalname][one_sigma][j][1] -res_toy[signalname][one_sigma][j][0]))
				thePull0p8.Fill(pull)
	        if i==4:
		        res_toy = pl_interval(model, 'toys:0.9', n=10100, cls = [one_sigma], signal_process_groups = signal_shapes , options = options)
		        for j in range(10000):
			        theBias0p9.Fill(res_toy[signalname][0][j])
				pull = (res_toy[signalname][0][j] - 0.9)/(0.5*(res_toy[signalname][one_sigma][j][1] -res_toy[signalname][one_sigma][j][0]))
				thePull0p9.Fill(pull)
	        if i==5:
		        res_toy = pl_interval(model, 'toys:1.0', n=10100, cls = [one_sigma], signal_process_groups = signal_shapes , options = options )
		        for j in range(10000):
			        theBias1p0.Fill(res_toy[signalname][0][j])
				pull = (res_toy[signalname][0][j] - 1.0)/(0.5*(res_toy[signalname][one_sigma][j][1] -res_toy[signalname][one_sigma][j][0]))
				thePull1p0.Fill(pull)
	        if i==6:
		        res_toy = pl_interval(model, 'toys:1.1', n=10100, cls = [one_sigma], signal_process_groups = signal_shapes  , options = options)
		        for j in range(10000):
			        theBias1p1.Fill(res_toy[signalname][0][j])
				pull = (res_toy[signalname][0][j] - 1.1)/(0.5*(res_toy[signalname][one_sigma][j][1] -res_toy[signalname][one_sigma][j][0]))
				thePull1p1.Fill(pull)
	        if i==7:
		        res_toy = pl_interval(model, 'toys:1.2', n=10100, cls = [one_sigma], signal_process_groups = signal_shapes  , options = options)
		        for j in range(10000):
			        theBias1p2.Fill(res_toy[signalname][0][j])
				pull = (res_toy[signalname][0][j] - 1.2)/(0.5*(res_toy[signalname][one_sigma][j][1] -res_toy[signalname][one_sigma][j][0]))
				thePull1p2.Fill(pull)
	        if i==8:
		        res_toy = pl_interval(model, 'toys:1.3', n=10101, cls = [one_sigma], signal_process_groups = signal_shapes  , options = options)
		        for j in range(10000):
			        theBias1p3.Fill(res_toy[signalname][0][j])
				pull = (res_toy[signalname][0][j] - 1.3)/(0.5*(res_toy[signalname][one_sigma][j][1] -res_toy[signalname][one_sigma][j][0]))
				thePull1p3.Fill(pull)
	        if i==9:
		        res_toy = pl_interval(model, 'toys:1.4', n=10100, cls = [one_sigma], signal_process_groups = signal_shapes , options = options)
		        for j in range(10000):
			        theBias1p4.Fill(res_toy[signalname][0][j])
				pull = (res_toy[signalname][0][j] - 1.4)/(0.5*(res_toy[signalname][one_sigma][j][1] -res_toy[signalname][one_sigma][j][0]))
				thePull1p4.Fill(pull)         
	        if i==10:
		        res_toy = pl_interval(model, 'toys:1.5', n=10101, cls = [one_sigma], signal_process_groups = signal_shapes  , options = options)
		        for j in range(10000):
			        theBias1p5.Fill(res_toy[signalname][0][j])
				pull = (res_toy[signalname][0][j] - 1.5)/(0.5*(res_toy[signalname][one_sigma][j][1] -res_toy[signalname][one_sigma][j][0]))
				thePull1p5.Fill(pull)
	       
	        print [ "%.4f" %res_toy[signalname][0][0] , "%.4f" %res_toy[signalname][one_sigma][two_sigma][0] , "%.4f" %res_toy[signalname][one_sigma][two_sigma][1] ] 



outputfile = TFile("results/output_pseudo.root", "recreate");
outputfile.cd()

if dobiasscan:
	theBias0p5.Write()
	theBias0p6.Write()
	theBias0p7.Write()
	theBias0p8.Write()
	theBias0p9.Write()
	theBias1p0.Write()
	theBias1p1.Write()
	theBias1p2.Write()
	theBias1p3.Write()
	theBias1p4.Write()
	theBias1p5.Write()

	thePull0p5.Write()
	thePull0p6.Write()
	thePull0p7.Write()
	thePull0p8.Write()
	thePull0p9.Write()
	thePull1p0.Write()
	thePull1p1.Write()
	thePull1p2.Write()
	thePull1p3.Write()
	thePull1p4.Write()
	thePull1p5.Write()
	
if dologprofile:
	resnll = nll_scan(model, 'data', n=1,  npoints=100, range=[0., 5.], options = options)
	finalplot = resnll[signalname][0]
	xprofil =  array("d", resnll[signalname][0].x)
	yprofil =  array("d", resnll[signalname][0].y)
	graph_profile = TGraph(len(xprofil), xprofil, yprofil)
	histo_graph = TH1F("histo_graph", "histo_graph", 100, xprofil[0], xprofil[len(xprofil)-1]);
	histo_graph.SetMaximum(graph_profile.GetMaximum());
	histo_graph.Write()
	graph_profile.Write()

if dodiscovery:
	#discovery(model, options = options, use_data=True)
	#print ("-----------------use data true -----------------")
	discovery(model, options = options, input_expected='toys:1.0', Z_error_max=0.1, use_data=True, nuisance_constraint=fixed_dist, verbose=True)
	print ("-----------------use data false-----------------")
	#discovery(model, options = options, input_expected='toys:1.0', Z_error_max=0.1, use_data=False, nuisance_constraint=fixed_dist)

if dolimits:
	if dolimits_bayes:
		plot_exp_bayes, plot_obs_bayes = bayesian_limits(model, 'all', n_toy = 10000, n_data = 20,nuisance_constraint=fixed_dist)
		plot_exp_bayes.write_txt('results/bayesian_limits_expected.txt')
		plot_obs_bayes.write_txt('results/bayesian_limits_observed.txt')
	if dolimits_cls:
		plot_exp_cls, plot_obs_cls  = cls_limits(model,  nuisance_prior = fixed_dist, options = options)
		plot_exp_cls.write_txt('results/cls_limits_expected.txt')
		plot_obs_cls.write_txt('results/cls_observed.txt')




report = model_summary(model)
