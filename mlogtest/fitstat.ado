*! 4.1.10 | 2017-01-05 | freese long | fixed bug in tobit (for v <14 and 15+)
* 4.1.9 | 2017-01-05 | freese long | fixed bug in tobit (for v 15 only)
* 4.1.8 | 2015-07-14 | freese long | use double types; improve documentation
* 4.1.7 | 2015-07-14 | freese long | esample bug
* 4.1.6 | 2015-03-24 | freese long | regress bug
* 4.1.5 | 2015-01-26 | freese long | dif bug
* 4.1.4 | 2014-09-22 cert | freese long | 11.2 and gologit2
* 4.1.3 | 2014-07-05 | freese long | fixed fweight bug
* 4.1.2 | 2014-05-12 | freese long | fixed sample size issue
* 4.1.1 | 2014-02-28 | freese long | spost13 release (reviewed)

program define fitstat, rclass

	* if before version 11, pass over to old fitstat
	if _caller() < 11 {
        local opts "`0'"
        local opts : subinstr local opts "spost9" "", all
		fitstat9 `opts'
		exit
	}

	local esample "e(sample) == 1" // flag for esample
	
	local cmd  "`e(cmd)'"

	* use old fitstat for some commands not supported by new version
	if "`cmd'"=="clogit" | "`cmd'"=="rologit" | ///
		"`cmd'"=="cnreg" | "`cmd'"=="gologit" {
        local opts "`0'"
        local opts : subinstr local opts "spost9" "", all
		fitstat9 `0'
		exit
	}

	version 11.2

//  parse syntax

    syntax [if] [in] ///
        [, Saving(string) Using(string) /// save/retrieve by name
		SAVE Diff /// save/retrieve without name
		FORCE ///
		Bic ic /// only ic measures (synonyms)
        Debug(integer 0) /// what does this do?
        spost9 /// run fitstat9
        ]

	* prohibit if svy used
	if "`e(prefix)'" == "svy" {
		di as err "quantities cannot be computed when svy is used"
		error 999
	}

	* prohibit if not on list of allowed commands
    if ( "`cmd'"!="regress"  & "`cmd'"!="logit"    & "`cmd'"!="probit" ///
       & "`cmd'"!="logistic" & "`cmd'"!="cloglog"  & "`cmd'"!="ologit" ///
       & "`cmd'"!="oprobit"  & "`cmd'"!="mlogit" ///
       & "`cmd'"!="tobit"    ///
       & "`cmd'"!="intreg"   & "`cmd'"!="poisson"  & "`cmd'"!="nbreg"  ///
       & "`cmd'"!="zip"      & "`cmd'"!="zinb"     ///
       & "`cmd'"!="ztp"      & "`cmd'"!="ztnb"     & "`cmd'"!="slogit" ///
       & "`cmd'"!="mprobit"  & "`cmd'"!="tnbreg" & "`cmd'"!="tpoisson" ///
       & "`cmd'"!="gologit2" ///
       ) {
        display as error ///
        "fitstat does not work with the last model estimated"
        exit
    }

    if "`spost9'"=="spost9" {
        local opts "`0'"
        local opts : subinstr local opts "spost9" "", all
        fitstat9 `opts  '
        exit
    }

	* allow ic as synonym for bic
	if ("`ic'"!="") local bic "bic"

	* define macros for saving and using options
    local SaveName "`saving'"
    local UsingName "`using'"
    if ("`SaveName'"=="" & "`save'" != "") local SaveName 0
    if ("`diff'"!="") local UsingName 0
		tempname UsingResults

//  COMPUTATIONS

    * classify model
	local dv = "`e(depvar)'" // name of dependent variable
	local isPAIweight "no" // pweight, aweight, or iweight?
	if "`e(wtype)'"!="" { 
		local wtype = "`e(wtype)'"
        local wtis "[`e(wtype)'`e(wexp)']"
		if 	"`wtype'" == "pweight" | "`wtype'" == "aweight" | ///
			"`wtype'" == "iweight" {
			local isPAIweight "yes"
		}
    }

	* binary outcome model
		local isBinary "no"
		if ("`cmd'"=="probit")      local isBinary "yes"
		if ("`cmd'"=="logit")       local isBinary "yes"
		if ("`cmd'"=="logistic")    local isBinary "yes"
		if ("`cmd'"=="cloglog")     local isBinary "yes"


    * number of observations
	tempname N
	cap sca `N' = e(N)

    * number of parameters
	tempname Nparm
	cap sca `Nparm' = e(rank)

    * log likelihood and deviance
	tempname ll
	cap sca `ll' = e(ll)

	if (`ll' != .) {
		tempname dev
		sca `dev' = -2*`ll'
		cap local dev_df = `N' - `Nparm'
	}

    * log likelihood of null model
	tempname ll_0
	cap sca `ll_0' = e(ll_0)

	* special for zip/zinb -- what stata does in not what we want so instead
	* must estimate a model with two intercepts
	if ("`cmd'" == "zip" | "`cmd'" == "zinb") {
		tempname insample realresults
		quietly {
			* recover information so can re-run model with two intercepts
			local inftyp "`e(inflate)'"
			if "`inftyp'"=="logit" local inftyp "" // null if logit
			* construct line of zip/zinb to execute
			local doit ///
				"`cmd' `dv' if e(sample)==1 `wtis' , inf(_cons) `inftyp'"
			* re-execute zip/zinb with two intercepts
			 estimates store `realresults'
				`doit'
				scalar `ll_0' = e(ll) // NOTE: not e(ll_0)
			estimates restore `realresults'
			estimates drop `realresults'
		}
	}


//  R2

	tempname r2
	cap sca `r2' = e(r2)

//  Adjusted R2

	tempname r2_adj
	cap sca `r2_adj' = e(r2_a)

//  LR/Wald chi-square

	* degrees of freedom
	cap local x2_df = e(df_m) // local not scalar so as to print properly
	* note: different definition of df used for zip/zinb
	if "`cmd'" == "zip"  	local x2_df = e(rank) - 2
	if "`cmd'" == "zinb" 	local x2_df = e(rank) - 3


	* test statistic
	tempname lrx2 waldx2
	if "`e(chi2type)'" == "LR" {
		cap sca `lrx2' = e(chi2)
	}
	if "`e(chi2type)'" == "Wald" {
		cap sca `waldx2' = e(chi2)
		if (`waldx2' != .) {
		}
	}
	* if zip/zinb then a LRX2 calculated using our definition of ll_0
	if "`cmd'" == "zip" | "`cmd'" == "zinb" {
		cap sca `waldx2' = .
		sca `lrx2' = 2*(`ll' - `ll_0')
	}

	* p-value
	tempname x2p
	if "`e(chi2type)'" == "LR" | "`e(chi2type)'" == "Wald" {
	    cap sca `x2p' = e(p)
	}

//  McFadden's R2 / McFadden's adjusted R2

	tempname r2_mf r2_mfadj
	cap sca `r2_mf' = 1-(`ll'/`ll_0')
	cap sca `r2_mfadj' = 1-((`ll'-`Nparm')/`ll_0')

//  Maximum likelihood R2	/ Cragg-Uhler (Nagelkerke)

	tempname r2_ml r2_cu m_r2_ml m_r2_cu
	if "`isPAIweight'" == "no" {
		* ML R2
		cap scalar `r2_ml' = 1 - exp(2*(`ll_0'-`ll')/`N')
		* Cragg-Uhler R2
		cap scalar `r2_cu' = (`r2_ml')/(1-exp(2*`ll_0'/`N'))
	}

//  BIC/AIC

	quietly cap estat ic // returns matrix r(S)

	if _rc == 0 {

		* BIC degrees of freedom (k)
		local bic_k = el(r(S), 1, 4) // local not scalar so it prints properly

		* BIC (Stata)
		tempname bic
		scalar `bic' = el(r(S), 1, 6)

		* BIC-deviance and BIC-prime
		tempname bic_dev bic_p
		scalar `bic_dev' = (-2*`ll')-((`N'-`bic_k')*ln(`N'))
		scalar `bic_p' = (-2*(`ll'-`ll_0'))+(`x2_df'*ln(`N'))

		* AIC
		tempname aic aicbyn
		scalar `aic' = el(r(S), 1, 5)
		if "`isPAIweight'" == "no" {
			scalar `aicbyn' = `aic' / `N'
		}

	}

//  Tjur and Efron (binary models)

	if "`isBinary'" == "yes" & ///
		("`isPAIweight'" == "no") {

		* probabilities of a 0 and 1
		tempvar pr0 pr1 // probability of a 0/1
		qui predict double `pr1' if `esample', p
		qui gen double `pr0' = 1 - `pr1'

		* dependent variable as 0 and 1 (in case is 0 and !0)
		tempvar dvbinary
		qui gen `dvbinary' = 1 - (`e(depvar)' == 0) if `esample'

		* Tjur's coefficient of determination
		tempname tjur
		tempvar tjur_1 tjur_0
		qui sum `pr1' `wtis' if `dvbinary' == 1 & `esample'
		sca `tjur_1' = r(mean)
		qui sum `pr1' `wtis' if `dvbinary' == 0 & `esample'
		sca `tjur_0' = r(mean)
		sca `tjur' = abs(`tjur_1' - `tjur_0')

		* Efron's R2
		tempname r2_ef
		tempvar  squareDifPr1 squareDifMean
		tempname dvmean sumSquareDifPr1 sumSquareDifMean
		qui sum `dvbinary' `wtis' if `esample'
		sca `dvmean' = r(mean)
		qui gen double `squareDifPr1' = (`dvbinary' - `pr1')^2 if `esample'
		qui gen double `squareDifMean' = (`dvbinary' - `dvmean')^2 if `esample'
		qui sum `squareDifPr1' `wtis' if `esample'
		sca `sumSquareDifPr1'= r(sum)
		qui sum `squareDifMean' `wtis' if `esample'
		sca `sumSquareDifMean' = r(sum)
		sca `r2_ef' = 1 - (`sumSquareDifPr1'/`sumSquareDifMean')
	}

//  count R2 and adjusted count R2 for ordinal/nominal

	if ("`isPAIweight'" == "no") {

		tempvar prediction correct currentprob maxprob
		tempname baseline catcount

		* binary models (uses code from above)
		if "`isBinary'" == "yes" {
			qui gen double `prediction' = (`pr1' > .5) if `esample'
			qui gen double `correct' = (`prediction' == `dvbinary') if `esample'
			sca `baseline' = `dvmean'
			if `dvmean' < .5 {
				sca `baseline' = 1-`dvmean'
			}
		}

		if 	("`cmd'" == "ologit" | "`cmd'" == "oprobit" | ///
			 "`cmd'" == "mlogit" | "`cmd'" == "mprobit" | ///
			 "`cmd'" == "gologit2" ) ///
			{

			qui levelsof `dv' if `esample'
			qui gen `maxprob' = 0 if `esample' // will set to maximum during loop
			sca `baseline' = 0
			qui gen `prediction' = . if `esample'
			foreach i in `r(levels)' {
				capture drop `currentprob'
				qui predict `currentprob' if `esample', outcome(`i')
				qui replace `prediction' = `i' if `esample' & `currentprob' > `maxprob'
				qui replace `maxprob' = `currentprob' if `esample' & `currentprob' > `maxprob'
				qui count if `dv' == `i' & `esample'
				sca `catcount' = r(N)
				if `catcount'/`N' > `baseline' {
					sca `baseline' = `catcount'/`N'
				}
			}

			qui gen `correct' = (`prediction'==`dv')
		}

		tempname r2_ct
		cap qui sum `correct' `wtis' if `esample'
		if _rc == 0 {
			sca `r2_ct' = r(mean)
		}

		tempname r2_ctadj
		cap sca `r2_ctadj' = (`r2_ct' - `baseline') / (1 - `baseline')

	}

//  variance of y* and e

	* variance of xB
	tempname xb Vxb temp
	if "`cmd'" == "logit" | "`cmd'" == "logistic" | "`cmd'" == "ologit" | ///
		"`cmd'" == "probit" | "`cmd'" == "oprobit" {
		quietly {
			cap predict `xb' if `esample', xb
			if _rc == 0 {
				estimates store `temp'
				mean `xb' `wtis' if `esample'
				estat sd
				cap sca `Vxb' = el(r(variance), 1, 1)
				estimates restore `temp'
				estimates drop `temp'
			}
		}
	}


	* variance of e
	tempname Verr
	if "`cmd'" == "logit" | "`cmd" == "logistic" | "`cmd'" == "ologit" {
		sca `Verr' = (_pi^2)/3
	}
	if "`cmd'" == "probit" | "`cmd'" == "oprobit" {
		sca `Verr' = 1
	}
	if "`e(cmd)'"=="intreg" {
        sca `Verr' = e(sigma)*e(sigma)
    }
    if "`e(cmd)'"=="tobit" {
		tempname b
		mat `b' = e(b)
		
		** NOTE: Stata implemented a change in tobit between versions 14 and 15
		** in which Stata returns the variance in the estimation output in 15+
		** as opposed to sigma in earlier version.  -fitstat- needs to detect
		** which version was used to fit the tobit model to decide whether to
		** square the quantity in this column.  AFAIK, this can only be done by 
		** looking at the labels in the column headers.
		
		local tobitnames : coleq `b'
		local _colsofb = colsof(`b')
		local teststring : word `_colsofb' of `tobitnames'
		if "`teststring'" == "sigma" { // pre-15 version of tobit used
			sca `Verr' = `b'[1, colsof(`b')]^2
		}
		else { // 15+ version of tobit used
			sca `Verr' = `b'[1, colsof(`b')]
		}		
		
//        sca `Verr' = `b'[1,"sigma:_cons"]*`b'[1,"sigma:_cons"]
    }

	* variance of y-star
	tempname Vystar
	cap sca `Vystar' = `Vxb' + `Verr' // Verr only exists for supported models

//  McKelvey & Zavoina R2

	tempname r2_mz
	cap sca `r2_mz' = (`Vxb'/`Vystar')

//  RETURN RESULTS

    * arrange statistics for returns & output
	local ancillary "N Nparm" // unprinted but still returned
	local likelihoods "ll ll_0"
	local deviancechi2 "dev lrx2 waldx2 x2p"
	local bicdefault "aic aicbyn bic"
	local bicadditional "bic_dev bic_p"
	local r2s "r2 r2_adj r2_mf r2_mfadj r2_mz r2_ml r2_cu r2_ef tjur"
	local countr2s "r2_ct r2_ctadj"
	local variances "Verr Vystar"
	local dfs "dev_df x2_df bic_k"

    * return computed quantities -- also make matrix of all returns
	tempname ResultsToSave

	* stored as scalars
	foreach type in `ancillary' `likelihoods' `deviancechi2' `bicdefault' ///
		`bicadditional' `r2s' `countr2s' `variances' `dfs' {

		cap return scalar `type' = ``type''	// only returns if it exists
		tempname tempmat // erases tempmat each time through loop
		cap mat `tempmat' = ``type'' // only saves if it exists
		if _rc != 0 {
			mat `tempmat' = .
		}
		mat rownames `tempmat' = "`type'"
		mat `ResultsToSave' = nullmat(`ResultsToSave') \ `tempmat'
	}
	mat colnames `ResultsToSave' = "`cmd'"

    * combine matrix with using matrix if dif option specified
	tempname CombinedResults Differences
	if "`UsingName'" != "" { // will have name if dif option specified

		* create matrix with differences between saved and using results
		mat `Differences' = `ResultsToSave' - _fitstat_`UsingName'
		mat colnames `Differences' = "Difference"

		* combine saved, using, and difference matrices
		mat `CombinedResults' = ///
			`ResultsToSave' , _fitstat_`UsingName', `Differences'
	}

//  compute lrtest for differences

	tempname lrtestp
	if "`UsingName'" != "" {
		cap lrtest . _fsest_`UsingName', `force'
		if _rc == 0 {
			tempname _p _chi2 _df
			sca `lrtestp' = r(p)
		}
	}

//  save matrix

	if "`SaveName'" != "" {
		mat _fitstat_`SaveName' = `ResultsToSave'
		estimates store _fsest_`SaveName'
	}

//  OUTPUT

    * ASSIGN NAMES FOR OUTPUT
	tempname saved dif
*	local eq_N "Number of"
*	local nm_N "observations"
	local eq_N "N"
	local nm_N "N"
	local eq_Nparm "Number of"
	local nm_Nparm "parameters"
	local eq_ll "Log-likelihood"
	local nm_ll "Model"
	local eq_ll_0 "Log-likelihood"
	local nm_ll_0 "Intercept-only"
	local eq_dev "Chi-square"
	local nm_dev "Deviance (df=`dev_df')"
	local eq_lrx2 "Chi-square"
	local nm_lrx2 "LR"
	if "`UsingName'" == "" local nm_lrx2 "LR (df=`x2_df')"
	local eq_waldx2 "Chi-square"
	local nm_waldx2 "Wald"
	if "`UsingName'" == "" local nm_waldx2 "Wald (df=`x2_df')"
	local eq_x2p "Chi-square"
	local nm_x2p "p-value"
	local eq_aic "IC"
	if "`ic'" != "" local eq_aic "AIC"
	local nm_aic "AIC"
	local eq_aicbyn "IC"
	if "`ic'" != "" local eq_aicbyn "AIC"
	local nm_aicbyn "AIC divided by N"
	if "`ic'" != "" local nm_aicbyn "(divided by N)"
	local eq_bic "IC"
	if "`ic'" != "" local eq_bic "BIC"
	local nm_bic "BIC (df=`bic_k')"
	local eq_bic_dev "IC"
	if "`ic'" != "" local eq_bic_dev "BIC"
	local nm_bic_dev "BIC (based on deviance)"
	local eq_bic_p "IC"
	if "`ic'" != "" local eq_bic_p "BIC"
	local nm_bic_p "BIC' (based on LRX2)"
	local eq_r2 "R2"
	local nm_r2 "R2"
	local eq_r2_adj "R2"
	local nm_r2_adj "Adjusted R2"
	local eq_r2_mf "R2"
	local nm_r2_mf "McFadden"
	local eq_r2_mfadj "R2"
	local nm_r2_mfadj "McFadden (adjusted)"
	local eq_r2_mz "R2"
	local nm_r2_mz "McKelvey & Zavoina"
	local eq_r2_ml "R2"
	local nm_r2_ml "Cox-Snell/ML"
	local eq_r2_cu "R2"
	local nm_r2_cu "Cragg-Uhler/Nagelkerke"
	local eq_r2_ef "R2"
	local nm_r2_ef "Efron"
	local eq_tjur "R2"
	local nm_tjur "Tjur's D"
	local eq_r2_ct "R2"
	local nm_r2_ct "Count"
	local eq_r2_ctadj "R2"
	local nm_r2_ctadj "Count (adjusted)"
	local eq_Vxb "Variance of"
	local nm_Vxb "xb"
	local eq_Verr "Variance of"
	local nm_Verr "e"
	local eq_Vystar "Variance of"
	local nm_Vystar "y-star"

	* detect differences in sample size
	if "`UsingName'" != "" {
		mat `tempmat' = `CombinedResults'["N", 1...]
		local difN = el(`tempmat', 1, 3)
		if `difN' != 0 { // if sample sizes are different
			if "`force'" == "" { // return error if force not specified
				di as err ///
					"different Ns between saved and current model (must use -force- option)"
				error 999
			}
			else { // if force is specified, display N
				local displayN "N"
			}
		}
	}

	* handle presenting degrees of freedom
	if "`UsingName'" != "" {
		tempname tempmat
		* deviance df
		mat `tempmat' = `CombinedResults'["dev_df", 1...]
		local currentdf = el(`tempmat', 1, 1)
		local saveddf = el(`tempmat', 1, 2)
		local difdf = el(`tempmat', 1, 3)
		local nm_dev "D (df=`currentdf'/`saveddf'/`difdf')"
		* bic
		mat `tempmat' = `CombinedResults'["bic_k", 1...]
		local currentdf = el(`tempmat', 1, 1)
		local saveddf = el(`tempmat', 1, 2)
		local difdf = el(`tempmat', 1, 3)
		local nm_bic "BIC (df=`currentdf'/`saveddf'/`difdf')"
		* lr chi2
		mat `tempmat' = `CombinedResults'["x2_df", 1...]
		local currentdf = el(`tempmat', 1, 1)
		local saveddf = el(`tempmat', 1, 2)
		local difdf = el(`tempmat', 1, 3)
		local nm_lrx2 "LR (df=`currentdf'/`saveddf'/`difdf')"
		* wald chi2
		mat `tempmat' = `CombinedResults'["x2_df", 1...]
		local currentdf = el(`tempmat', 1, 1)
		local saveddf = el(`tempmat', 1, 2)
		local difdf = el(`tempmat', 1, 3)
		local nm_waldx2 "Wald (df=`currentdf'/`saveddf'/`difdf')"
	}

    * determine which quantities will be printed
	if "`ic'" == "" {
		local PrintList ///
		"`likelihoods' `deviancechi2' `r2s' `countr2s' `bicdefault' `variances' `displayN'"
	}
	else {
		local PrintList ///
			"`bicdefault' `bicadditional' `displayN'"
	}

    * print quantities
	tempname ResultsToPrint

	* if using is not specified
	foreach type in `PrintList' {

		* set row of output
		tempname tempmat // erases each time through loop
		if "`UsingName'" == "" {
			cap mat `tempmat' = ``type''
			if _rc != 0 {
				mat `tempmat' = .
			}
		}
		else {
			mat `tempmat' = `CombinedResults'["`type'", 1...]
		}

		if el(`tempmat', 1, 1) != . { // has to exist in current results
			mat rownames `tempmat' = "`nm_`type''"
			mat roweq `tempmat' = "`eq_`type''"
			mat `ResultsToPrint' = nullmat(`ResultsToPrint') \ `tempmat'
		}
	}

	if "`UsingName'" == "" {
		mat colnames `ResultsToPrint' = "`cmd'"
	}
	else {
		mat colnames `ResultsToPrint' = "Current" "Saved" "Difference"
	}

    * LR TEST
	if "`UsingName'" != "" & "`ic'" == "" {
		* first change the value in question to missing
		cap mat `ResultsToPrint'[rownumb(`ResultsToPrint', "`nm_x2p'"), 3] = .
		cap mat `ResultsToPrint'[rownumb(`ResultsToPrint', "`nm_x2p'"), 3] ///
				= `lrtestp'
		if _rc == 0 {
				mat `tempmat' = `CombinedResults'["x2_df", 1...]
				local difdf = el(`tempmat', 1, 3)
				if `difdf' > 0 {
					local lrtestmessage ///
						"Note: Likelihood-ratio test assumes saved model nested in current model."
				}
				else {
					local lrtestmessage ///
						"Note: Likelihood-ratio test assumes current model nested in saved model."
				}
		}
	}

    * print matrix of results
	matlist `ResultsToPrint', twidth(24) format(%11.3f)

	if "`e(vce)'" == "robust" {
		di _n "Note: Some measures based on pseudolikelihoods."
	}

	if "`lrtestmessage'" != "" {
		di _n "`lrtestmessage'"
	}

    * BIC COMPARISONS
    if "`UsingName'" != "" {

		tempname tempmat savedbic bicdif savedll
		mat `tempmat' = `CombinedResults'["bic", 1...]
		sca `savedbic' = el(`tempmat', 1, 2) // get saved bic
		sca `bicdif' = abs(el(`tempmat', 1, 3)) // get difference in bic
		mat `tempmat' = `CombinedResults'["ll_0", 1...]
		sca `savedll' = el(`tempmat', 1, 2) // get saved ll_0

		if `bicdif' != . { // skip even if force specified
			if (abs(`ll'-`savedll') > .0001) /// only do if ll_0s are the same
			| "`force'" != ""            ///  or force option is specified
			{
				local support "very strong"
				if `bicdif' <= 2 		local support "weak"
				else if `bicdif' <= 6 	local support "positive"
				else if `bicdif' <= 10	local support "strong"
				local preferred "current"
				if `savedbic' < `bic'	local preferred "saved"
				if `bicdif' < .000000001 {
					di _n "No difference in BIC between models."
				}
				else {
					di _n "Difference of " %8.3f in y `bicdif' ///
					   in g " in BIC provides " in y "`support'" ///
					   in g " support for " in y "`preferred'" in g " model."
				}
			}
		}
    }

end
exit

/*

Notes:

*  	y* for tobit and intreg could be added for this command but
	have not been implemented yet

*	Commands outside those references in book are not tested for this command.
	Users are responsible for interpretation.

/*




*/
