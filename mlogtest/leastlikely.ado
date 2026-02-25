*! 3.1.0 2014-02-14 | freese long | spost13 release

//  least likely predictions

program define leastlikely
    version 11.0

    syntax [varlist(default=none fv)] [if] [in] [, n(integer 5) Generate(string) *]

    if "`e(cmd)'"=="clogit" | "`e(cmd)'"=="nlogit" | "`e(cmd)'"=="xtlogit" | /*
    */ "`e(cmd)'"=="blogit" | "`e(cmd)'"=="bprobit" | "`e(cmd)'"=="glogit" | /*
    */ "`e(cmd)'"=="gprobit" {
        di _n as err "leastlikely not intended for use after `e(cmd)'"
        exit 198
    }

    tempname prob touse
    local depvar = "`e(depvar)'"
    local depvarlabel : value label `depvar'

    quietly {
        gen `touse' = e(sample) `if' `in'
        gen `prob' = .

		** get values of outcome

			levelsof `depvar' if e(sample)==1, local(values)
			local numcats : word count `values'
			forvalues i = 1(1)`numcats' {
				local value`i' : word `i' of `values'
			}

		** determine if binary variable

			local isbin "no"
			if `numcats'==2 & `value1'==0 & `value2'==1 {
				local isbin "yes"
			}

		di "isbin is `isbin'"

		** predict outcome that occurred (different syntax for binary variables)

			if "`isbin'" == "yes" {
				tempname temp
				predict `temp', p
				replace `prob' = `temp' if `depvar'==1 & `touse'==1
				replace `prob' = (1-`temp') if `depvar'==0 & `touse'==1
			}

			if "`isbin'" == "no" {
				forvalues i = 1(1)`numcats' {
					tempname temp
					predict `temp', outcome(`value`i'')
					replace `prob' = `temp' if `depvar'==`value`i'' & `touse'==1
				}
			}

    }

    if "`generate'"=="" {
        local generate "Prob"
		local erase "yes"
    }
	quietly rename `prob' `generate'

    forvalues i = 1(1)`numcats' {
        local vallabel ""
        if "`depvarlabel'"!="" {
            local vallabel : label `depvarlabel' `value`i''
            if "`vallabel'"!="" {
				local vallabel = "(`vallabel')"
			}
        }
        di _n as txt "Outcome: " as res `value`i'' " `vallabel'"

        tempname temp
        quietly egen `temp' = rank(`generate') if `depvar'==`value`i'' & `touse'==1, track
        list `generate' `varlist' if `temp'<=`n' , `options'
    }


    if "`erase'"=="yes" {
		drop `generate'
	}

end
exit
