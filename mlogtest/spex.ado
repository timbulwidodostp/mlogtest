*! version 4.0.9 2014-08-19 | long freese | allow any valid data name

capture program drop spex
program define spex

	* online data location
	local data_url "http://www.indiana.edu/~jslsoc/stata/spex_data/"

	/*
	SPECIFICATIONS OF DIFFERENT MODELS
	details of specific models/examples should be handled here
	[cmd]_data -> 	name of dataset to be used
	[cmd]_model -> 	model to be estimated
	[cmd]_pre# -> 	any commands to be executed between loading data
					and estimating model
	[cmd]_post# -> 	any commands to be executed between loading data
					and estimating model

	post-estimation commands can be included -- see below
		for post-estimation commands, options will be executed with first cmd
		after logit

	*/

	* binary
	foreach cmd in logit probit cloglog {
		local `cmd'_data "binlfp4"
		local `cmd'_model "`cmd' lfp k5 k618 i.agecat i.wc i.hc lwg inc"
	}

	* ordinal
	foreach cmd in ologit oprobit {
		local `cmd'_data "gssclass4"
		local `cmd'_model ///
			"`cmd' class i.female i.white i.year i.educ c.age##c.age income"
	}

	* nominal
	foreach cmd in mlogit mprobit {
		local `cmd'_data "partyid4"
		local `cmd'_model "`cmd' party age income i.black i.female i.educ"
	}

	* count
	foreach cmd in poisson nbreg {
		local `cmd'_data "couart4"
		local `cmd'_model "`cmd' art i.female i.married kid5 phd mentor"
	}

	* regress
	foreach cmd in regress reg {
		local `cmd'_data "regjob3"
		local `cmd'_model "`cmd' job i.fem phd ment i.fel art cit"
	}

	* truncated count
	foreach cmd in ztnb ztp { // earlier version
		local `cmd'_data "`poisson_data'"
		local `cmd'_model "`cmd' art i.female i.married kid5 phd mentor if art > 1"
	}
	foreach cmd in tpoisson tnbreg { // later version
		local tcou_data "`ztp_data'"
		local tcou_model "`cmd' `ztp_model', ll(0)"
	}

	* zero-inflated
	foreach cmd in zip zinb {
		local `cmd'_data "`poisson_data'"
		local zi_counteq "i.female i.married kid5 phd mentor" // allow to be different
		local zi_infeq "i.female i.married kid5 phd mentor"
		local `cmd'_model "`cmd' art `zi_counteq', inflate(`zi_infeq')"
	}

	* slogit
	local slogit_data "`nom_data'"
	local slogit_model "slogit party age income i.black i.female i.educ"

	* asclogit
	local asclogit_data "travel4"
	local asclogit_model "asclogit choice time, alt(mode) case(id) nolog"

	* asmprobit
	local asmprobit_data "travel4"
	local asmprobit_model ///
		"asmprobit choice time, case(id) alternatives(mode) casevars(hinc psize) nolog"

	* svy example
	local svy_data "svyhrs4"
	local svy_pre1 `"svyset secu [pweight=kwgtr], strata(stratum)"'
    local svy_model "svy: logit arthritis male ib3.ed3cat c.age##c.age"

	* rologit example
    local rologit_data "wlsrank4"
    local rologit_model "rologit rank high low ib4.jobchar##(c.score female), group(id) reverse noomitted"

	* tobit
	local tobit_data "tobjob2"
    local tobit_model "tobit jobcen i.fem phd ment i.fel art cit, ll(1)"

	* intreg
    local intreg_data "nels_censored2"
    local intreg_model "intreg minscor maxscor bymomed bydaded black hispanic"

	// POST-ESTIMATION COMMANDS

	* set data to be the logit data
	foreach cmd in fitstat listcoef margins mchange mtable mgen mlincom {
		*cmd_post="yes" tells stata to associate spex options with
		*first cmd after est cmd, not with est cmd
		local `cmd'_post "yes"
		local `cmd'_data "`logit_data'"
		local `cmd'_model "quietly `logit_model'"
	}

	* fitstat
	local fitstat_post1 "fitstat"

	* listcoef
	local listcoef_post1 "listcoef"

	* margins
	local margins_post1 "margins"

	* mchange
	local mchange_post1 "mchange"

	* mtable
	local mtable_post1 "mtable, at(wc=(0 1))"

	* mlincom
	local mlincom_post1 "estimates store mymodel"
    local mlincom_post2 "mtable, at(wc=(0 1)) post"
	local mlincom_post3 "mlincom 1 - 2"
	local mlincom_post4 "estimates restore mymodel"

	* mgen
	local mgen_post1 "mgen, at(inc=(0(10)100)) atmeans"
    local mgen_post2 graph twoway ///
        (rarea _ll1 _ul1 _inc, col(gs12)) ///
        (connected _pr1 _inc, lpat(solid) msym(i)), ///
        ylab(0(.25)1) ytitle(Pr(LFP)) ///
        xlab(0(10)100) xtitle(Income)

	/*
	EXECUTION
	no code below should be specific to a particular model
	*/

	* parameters that can be tweaked

	local look_first "local"

	* parse input

*	syntax name [, Web local user copy *]
*    tokenize "`namelist'"
    syntax anything [, Web local user copy *]
    tokenize "`anything'"
	local arg "`1'"

	local done "no"

	* options about locations
	if "`web'" != "" {
		local look_first "web"
	}
	if "`user'" != "" | "`local'" != "" {
		local look_first "local"
	}

	* handle case where argument is a model

	if "``arg'_data'" != "" {

		* get data
		_SpexUse ``arg'_data', url(`data_url') lookfirst(`look_first') `copy'

		* execute any pre-commands
		local precommands_done "no"
		local pre_count = 1
		while "`precommands_done'" == "no" {
			if "``arg'_pre`pre_count''" == "" {
				local precommands_done "yes"
			}
			else {
				local cmdtext "``arg'_pre`pre_count''"
				di as txt ". `cmdtext'"
				`cmdtext'
                window push `cmdtext'
				local pre_count = `pre_count' + 1
			}
		}

		* estimate model
		local cmdtext "``arg'_model'"
		if "`options'" != "" & "``arg'_post'" != "yes" { // if user has specified options
			if strpos("`cmdtxt'", ",") == 0 { // no options in default model
				local cmdtext "`cmdtext', `options'"
			}
			else { // are options in default model, so no comma
				local cmdtext "`cmdtext' `options'"
			}
		}
		di as txt ". `cmdtext'"
		`cmdtext'
		di
        window push `cmdtext'
		local done "yes"

		* execute any post-commands
		local postcommands_done "no"
		local post_count = 1
		while "`postcommands_done'" == "no" {
			if "``arg'_post`post_count''" == "" {
				local postcommands_done "yes"
			}
			else {
				local cmdtext "``arg'_post`post_count''"
				if "`options'" != "" & "``arg'_post'" == "yes" & `post_count' == 1 {
					if strpos("`cmdtxt'", ",") == 0 { // no options in default model
						local cmdtext "`cmdtext', `options'"
					}
					else { // are options in default model, so no comma
						local cmdtext "`cmdtext' `options'"
					}
				}

				di as txt ". `cmdtext'"
				`cmdtext'
				di
                window push `cmdtext'
				local post_count = `post_count' + 1
			}
		}

	}

	* if not handled as command, try loading as dataset

	if "`done'" == "no" {
		di as txt "(`arg' is not supported command; attempting to load as dataset)"
		_SpexUse `arg', url("`data_url'") lookfirst(`look_first') `copy'
		local done "yes"
	}

end

capture program drop _SpexUse
program define _SpexUse

*	syntax namelist, url(string) lookfirst(string) [ copy ]
*    tokenize "`namelist'"

    syntax anything, url(string) lookfirst(string) [ copy ]
    tokenize "`anything'"

	local data "`1'"
	local look_first "`lookfirst'"

// look first locally

	if "`look_first'" == "local" {

		capture sysuse `data'.dta, clear
		if _rc != 0 {
			capture use "`url'`data'", clear
			if _rc != 0 {
				di as err `"program cannot find file `data' locally or online"'
				error 999
			}
			else {
				di _n as txt `". use "`where'`data'.dta", clear"'
				di
                window push use `where'`data'.dta, clear
			}
		}
		else {
			di _n as txt `". sysuse "`data'.dta", clear"'
			di
            window push sysuse `where'`data'.dta, clear
		}

	}

//  look first online

	if "`look_first'" == "web" {

		capture use "`url'`data'", clear
		if _rc != 0 {
			 capture sysuse `data'.dta, clear
			 if _rc != 0 {
				  di as err `"program cannot find file `data' locally or online"'
				  error 999
			 }
			 else {
				 di _n as txt `". sysuse "`data'.dta", clear"'
                 window push sysuse `data'.dta, clear
			 }
		}
		else {
			di _n as txt `". use "`url'`data'.dta", clear"'
            window push use `url'`data'.dta, clear
		}

	}

//  copy to working directory

    if "`copy'"=="copy" {
        capture copy "`url'`data'.dta" `data'.dta
        di "note: file `data'.dta copied to working directory"
    }

end
exit
* version 4.0.8 2014-08-07 | long freese | copy fix
* version 4.0.7 2014-07-28 | long freese | no regjob3
* version 4.0.6 2014-07-01 | long freese | copy
* version 4.0.5 2014-06-28 | long freese | mgen
* version 4.0.4 2014-06-09 | long freese | window push
* version 4.0.3 2014-06-08 | long freese | additional examples
* version 4.0.1 2014-05-13 | long freese | full rewrite for spost13 release
