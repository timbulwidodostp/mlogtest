*! version 3.5.0 2016-10-25 | long freese | problem in Stata 14.2 smhsiao
//  mlogit tests
//      : small-hsiao test is based on Nick Winter's -smhsiao-

capture program drop mlogtest
program define mlogtest, rclass

    version 11.2
    tempvar sample tmp touse samp tempdepvar denom lnL
    tempname omit pval v rowres matsh Vals b0a b0b b0ab b1b
    tempname matiia matsuest matlr matwald matcomb matlrc nextrow testn newbase

    syntax [varlist(default=none fv)][, detail iia Hausman SMhsiao SUest ///
        Lr Wald COMBine LRCOMBine set(string) all NOISily ]
    local base "base" // base is now default always
    local Dst = 3 // decimals for stats
    local Dpv = 3 // decimals for pv

    if "`e(cmd)'" != "mlogit" {
        display as error "mlogtest only works after mlogit"
        exit
    }
    qui estimates store _Orig
    local n "(N=`e(N)')"

//  decode options

    if ("`*'"=="") local wald wald
    if "`iia'`hausman'`smhsiao'`suest'`lr'`combine'`lrcombine'`all'"=="" ///
        local wald wald

    local dolr = "no"
    if ("`lr'"!="") local dolr = "yes"
    local doiia = "no"
    if ("`hausman'"!="" | "`iia'"!="") local doiia = "yes"
    local dosmhsiao = "no"
    if ("`smhsiao'"!="" | "`iia'"!="") local dosmhsiao = "yes"
    local dosuest = "no"
    if ("`suest'"!="" | "`iia'"!="") local dosuest = "yes"
    local dowald = "no"
    if ("`wald'"!="") local dowald = "yes"
    local docombine = "no"
    if ("`combine'"!="") local docombine = "yes"
    local dolrcombine = "no"
    if ("`lrcombine'"!="") local dolrcombine = "yes"
    if "`all'"!="" {
        local dolr = "yes"
        local doiia = "yes"
        local dosmhsiao = "yes"
        local dosuest = "yes"
        local dowald = "yes"
        local docombine = "yes"
        local dolrcombine = "yes"
    }

//  get weight info from last mlogit

    local wtis ""
    if ("`e(wtype)'"!="") local wtis "[`e(wtype)'`e(wexp)']"

//  check that estimation sample matches n from regression

    qui gen `sample' = e(sample)
    if "`e(wtype)'"=="" | "`e(wtype)'"=="aweight" | "`e(wtype)'"=="pweight" {
        qui count if `sample' == 1
        scalar `testn' = r(N)
    }
    if "`e(wtype)'"=="fweight" | "`e(wtype)'"=="iweight" {
        local wtexp = substr("`e(wexp)'", 3, .)
        gen `tmp' = (`wtexp')*`sample'
        qui su `tmp', meanonly
        scalar `testn' = round(r(sum),1)
    }
    if e(N) != `testn' {
        display _new as error "data has been altered since mlogit was estimated"
        exit
    }

//  information about the mlogit

    local depvar = "`e(depvar)'"
    local printdepvar = abbrev("`printdepvar'", 8)
    _rm_modelinfo2
    local rhscmdline "`r(rhs_cmdline)'" // this will have i.catvar
    local rhsbetabase "`r(rhs_betanmsbase)'" // beta's including base vars
    local rhsbetanobase "`r(rhs_betanmsnobase)'" // without base vars

    local issvy = 0
    if "`r(cmdsvy)'"=="1" local issvy = 1
    local rhsnam "`r(rhsnms)'" // this will have 2.catvar 3.catvar...
    local rhsN = `r(rhsn)'
    local catnms8 "`r(lhscatnms8)'"
    local catvals "`r(lhscatvals)'"
    local catnms "`r(lhscatnms)'"
    qui `noisily' display "catvals `catvals'"
    qui `noisily' display `"catnms: `catnms' "'
    local catsN = `r(lhscatn)'
    qui `noisily' display "catsN: " `catsN'
    local basecat = e(baseout)
    local ipos : list posof "`basecat'" in catvals
    qui `noisily' display "basecat: `basecat'"
    local refnm : word `ipos' of `catnms'
    qui `noisily' display "refnm: `refnm'"
    local check : word count `catvals'
    if `catsN' != real("`check'") {
        display as error _new ///
        "problem determining number of outcome categories"
        exit
    }
    if `issvy' {
        local svylr = 0
        if "`dolr'"=="yes" {
            local svylr = 1
            local dolr "no"
        }
        if "`dolrcombine'"=="yes" {
            local svylr = 1
            local dolrcombine "no"
        }
        local dolf = 0
        if `svylr' {
            display "LR tests cannot be run with svy estimation"
            local dolf = 1
        }
        if "`doiia'"=="yes" | "`dosuest'"=="yes" | "`dosmhsiao'"=="yes" {
            local doiia "no"
            local dosuest "no"
            local dosmhsiao "no"
            display "IIA test are not available for svy estimation"
            local dolf = 1
        }
        if `dolf' display
    }

    qui `noisily' display "catsN: " `catsN'
    qui `noisily' display "basecat:  `basecat'"
    qui `noisily' display "refnm:    `refnm'"

//  parse sets

    if "`set'" != "" {
        * add spaces so symbols are words
        local set : subinstr local set "=" " = ", all
        local set : subinstr local set "\" " \ ", all
        local setN = wordcount("`set'")
        tokenize "`set'"
        local setis = 1
        forvalues iword = 1/`setN' {
            local word "``iword''"
            local nextword = `iword' + 1
            local nextword "``nextword''"
            if ("`nextword'"=="=") local set`setis'nm "`word'"
            else { // current is not set name
                if ("`word'"=="\") local ++setis
                else if ("`word'"=="=") local dog ""
                else local set`setis' "`set`setis''`word' "
            }
        }
        local numset = `setis'
        forvalues iset = 1/`numset' {
            if ("`set`iset'nm'"=="") local set`iset'nm "set_`iset'"
            foreach nm in `set`iset'' {
                local ipos : list posof "`nm'" in rhsnam
                if `ipos'==0 {
                    display as error "variable `nm' in set() is not in model"
                    exit
                }
            }
            return local set_`iset' "`set`iset''"
        }
    }
    else local numset = 0

//  LR test of independent variables

    if "`dolr'"=="yes" & `rhsN' == 0 {
        display _new in r "LR test cannot be computed for an intercept-only model"
    }
    else if "`dolr'"=="yes" {
        display _new as res "LR tests for independent variables `n'"
        display _new as text ///
        "  Ho: All coefficients associated with given variable(s) are 0"
        if "`varlist'"!="" {
            local nvlist : word count `varlist'
            local ntests = `nvlist' // use this as the # of rhsvars
            local testnms "`varlist'"
        }
        else {
            local testnms "`rhsnam'"
            local ntests = `rhsN'
        }
        tokenize "`rhsnam'"
        local count = 1
        local rspecis "&|"
        while `count' <= `ntests'+`numset' {
            local rspecis "`rspecis'&"
                * local testvar : word `count' of `testnms'
            if `count' <= `ntests' { // for single parameter
                local var`count' : word `count' of `testnms'

                * remove i form i#.name 2014-05-30

                local testvar "`var`count''"
                local nonum = regexr("`testvar'","[0-9]+\.",".")
                if "`testvar'"!="`nonum'" {
                * 2014-06-19 quick remove i in i1, i2 etc.
                    forvalues i = 1/99 {
                        local testvar : subinstr local testvar "i`i'" "`i'"
                    }
                    local var`count' `testvar'
                }
                * error for i.educ 2014-05-30
                local ipos : list posof "`var`count''" in rhsbetabase
                if `ipos'==0 {
                    display as error "`var`count'' is not a coefficient in the  model"
                    exit
                }
                local lrrhs : list rhsbetabase - var`count'
            }
            if `count' > `ntests' { // for sets of parameters
                local thisset = `count'-`ntests'
                local var`count' "set_`thisset'"
                local lrrhs : list rhsbetabase - set`thisset'
            }
            qui mlogit `depvar' `lrrhs' `wtis' if `sample' == 1, b(`basecat')
            qui lrtest . _Orig
            mat `nextrow' = r(chi2), r(df), r(p)
            mat rownames `nextrow' = "`var`count''"
            local setnum = `count' - `ntests'
            if (`setnum'>0) mat rownames `nextrow' = "`set`setnum'nm'"
            mat `matlr' = nullmat(`matlr') \ `nextrow'
            local ++count
        }
        mat colna `matlr' = chi2 df "P>chi2"
        qui estimates restore _Orig
        local cspecis "& %16s | %9.`Dst'f & %4.0f & %7.`Dpv'f &"
        matlist `matlr', cspec(`cspecis') rspec(`rspecis')
        return matrix lrtest `matlr'
        if `numset'>0 display
        forvalues iset = 1/`numset' {
            display as text "   `set`iset'nm' contains: `set`iset''"
        }

    } // if "`dolr'"=="yes"

//  WALD test of independent variables

    if "`varlist'"!="" {
        local nvlist : word count `varlist'
        local testnms "`varlist'"
        local ntests = `nvlist' // use this as the # of rhsvars
    }
    else {
        local ntests = `rhsN'
        local testnms "`rhsnam'"
    }
    if "`dowald'"=="yes" & `ntests' == 0 {
        display _n as error ///
        "Wald test cannot be computed on intercept-only model"
        exit
    }
    else if "`dowald'"=="yes" {
        tokenize "`testnms'"
        display _new as res "Wald tests for independent variables `n'"
        display _new as text ///
        "  Ho: All coefficients associated with given variable(s) are 0"
        local count = 1
        local rspecis "&|"
        while `count' <= `ntests'+`numset' {
            local rspecis "`rspecis'&"
            if `count' <= `ntests' { // individual variable
                local var`count' : word `count' of `testnms'
                qui test `var`count''
            }
            if `count' > `ntests' { // sets
                local thisset = `count'-`ntests'
                local var`count' "`set`thisset'nm'"
                qui test `set`thisset''
            }
            if `issvy'==0 {
                mat `nextrow' = r(chi2), r(df), r(p)
                mat rownames `nextrow' = "`var`count''"
                mat `matwald' = nullmat(`matwald') \ `nextrow'
                mat colnames `matwald' = chi2 df "P>chi2"
                local cspecis "& %16s | %9.`Dst'f & %4.0f & %7.`Dpv'f &"
            }
            else {
                local cspecis ///
                "& %16s | %9.`Dst'f & %4.0f & %5.0f & %7.`Dpv'f &"
                mat `nextrow' = r(F), r(df), r(df), r(p)
                mat rownames `nextrow' = "`var`count''"
                mat `matwald' = nullmat(`matwald') \ `nextrow'
                mat colnames `matwald' = F df df_r "P>F"
            }
            local ++count
        } // while `count' <= `ntests'
        matlist `matwald', cspec(`cspecis') rspec(`rspecis')
        return matrix wald `matwald'
        if `numset'>0 display
        forvalues iset = 1/`numset' {
            display as text "   `set`iset'nm' contains: `set`iset''"
        }
    } // if "`dowald'"=="yes"

//  IIA test

    if ("`doiia'"=="yes" | "`dosmhsiao'"=="yes" | "`dosuest'"=="yes") ///
            & `catsN' == 2 {
        di _n as error "IIA tests requires at least 3 dependent categories"
        exit
    }
    if "`doiia'"=="yes" {
        display _new as res "Hausman tests of IIA assumption `n'"
        display _new as text ///
"  Ho: Odds(Outcome-J vs Outcome-K) are independent of other alternatives"
        tokenize "`catvals'"
        local count = 1
        local anyneg "no"
        local rspecis "&|"
        while real("`count'") <= `catsN' {
            local rspecis "`rspecis'&"
            local lab`count' : word `count' of `catnms'
            local slab`count' : word `count' of `catnms8'
            scalar `omit' = real("``count''")
            if "`detail'" == "detail" {
                if "`lab`count''"!="``count''" {
                    display as text _new ///
                        "Hausman test with omitted alternative " ///
                        `omit' " (`lab`count'')"
                }
                else {
                    display as text _new ///
                        "Hausman test with omitted alternative " `omit'
                }
            }
            if `omit'==real("`basecat'") & "`base'"!="" {
                * IIA for basecat requires estimating mlogit with
                * new basecategory. Make new basecat the largest category
                * of the dependent variable that is not the original basecat
                local maxcnt = 0
                local count2 = 1
                while `count2' <= `catsN' {
                    if real("``count2''") != `omit' {
                        qui count if `depvar'==``count2'' & `sample'==1
                        if r(N)>`maxcnt' {
                            local newbase = ``count2''
                            local maxcnt = r(N)
                            local tmplab : word `count2' of `catnms'
                        }
                    }
                    local ++count2
                } // while count2 <= `catsN'
                qui mlogit `depvar' `rhscmdline' `wtis' ///
                        if `sample' == 1, b(`newbase')
                qui estimates store _OrigNewBase
                qui mlogit `depvar' `rhscmdline' `wtis' ///
                        if `sample'==1 & `depvar'!=`omit', b(`newbase')
                if "`detail'"=="detail" {
                    if "`lab`count''"!="``count''" {
                        display as text ///
"(Using category `newbase' (`tmplab') as comparison group)"
                    }
                    else {
                        display as text ///
                        "(Using category `newbase' as comparison group)"
                    }
                    hausman . _OrigNewBase, alleq constant
                }
                else qui hausman . _OrigNewBase, alleq constant
                scalar `pval' = r(p)
                if r(chi2)<0 {
                    local anyneg "yes"
                    scalar `pval' = .
                }
                mat `nextrow' = r(chi2), r(df), `pval'
                mat rownames `nextrow' = "`slab`count''"
                mat `matiia' = nullmat(`matiia') \ `nextrow'
                qui estimates restore _Orig
            } // if `omit'==real("`basecat'")
            else if `omit'!=real("`basecat'") {
                qui mlogit `depvar' `rhscmdline' `wtis' ///
                    if `sample' == 1 & `depvar' != `omit', b(`basecat')
                if "`detail'"=="detail" hausman . _Orig, alleq constant
                else qui hausman . _Orig, alleq constant
                scalar `pval' = r(p)
                if r(chi2)<0 {
                    local anyneg "yes"
                    scalar `pval' = .
                }
                mat `nextrow' = r(chi2), r(df), `pval'
                mat rownames `nextrow' = "`slab`count''"
                mat `matiia' = nullmat(`matiia') \ `nextrow'
            } // else if `omit'!=real("`basecat'")
            local ++count
        } // while real("`count'") <= `catsN'
        mat colnames `matiia' = chi2 df "P>chi2"
        local cspecis "& %16s | %9.`Dst'f & %4.0f & %7.`Dpv'f &"
        matlist `matiia', cspec(`cspecis') rspec(`rspecis')
        if ("`detail'"=="detail") display as res _new "Summary of results"
        display _new as text ///
        "  Note: A significant test is evidence against Ho."
        if "`anyneg'" == "yes" {
            display as text ///
"  Note: If chi2<0, the estimated model does not meet asymptotic assumptions."
        }
        qui estimates restore _Orig
        return matrix hausman `matiia'
    } // if "`doiia'"=="no"

//  suest-based HAUSMAN IIA Test

    * see if cluster used; if so incorporate in suest
    local clustvar = e(clustvar)
    if "`clustvar'" != "." & "`clustvar'" != "" {
        display _n as err ///
"suest-based Hausman IIA test must be estimated differently when cluster()"
        display as err "is specified; see [R] suest"
    }
    else if "`dosuest'"=="yes" {
        display as res _new "suest-based Hausman tests of IIA assumption `n'"
        display _new as text ///
        "  Ho: Odds(Outcome-J vs Outcome-K) are independent of other alternatives"
        tokenize "`catvals'"
        local count = 1
        local rspecis "&|"
        while real("`count'") <= `catsN' {
            local rspecis "`rspecis'&"
            local lab`count' : word `count' of `catnms'
            local slab`count' : word `count' of `catnms8'
            scalar `omit' = real("``count''")
            * base is option to run also with base category included
            if `omit'==real("`basecat'") & "`base'" != "" {
                * IIA for basecategory requires estimating mlogit with
                * new basecategory. Make new basecat the largest category
                * of the dependent variable that is not the original basecat
                local maxcnt = 0
                local count2 = 1
                while `count2' <= `catsN' {
                    if real("``count2''") != `omit' {
                        qui count if `depvar'==``count2'' & `sample'==1
                        if r(N)>`maxcnt' {
                            local newbase = ``count2''
                            local maxcnt = r(N)
                            local tmplab : word `count2' of `catnms'
                        }
                    }
                    local ++count2
                } // while count2 <= `catsN'
                * estimates of old logit are held in _Orig
                qui {
                    mlogit `depvar' `rhscmdline' `wtis' ///
                        if `sample' == 1, b(`newbase')
                    estimates store _OrigNewBase
                    mlogit `depvar' `rhscmdline' `wtis' ///
                        if `sample'==1 & `depvar'!=`omit', b(`newbase')
                    estimates store _Aux
                    local eqnames = e(eqnames)
                    estimates store _Aux
                    suest _Aux _OrigNewBase
                    estimates store _suest
                }
                foreach name in `eqnames' {
                    qui test [_OrigNewBase_`name' = _Aux_`name'], cons accum
                }
                if "`detail'"=="detail" {
                    suest
                    test
                }
                mat `nextrow' = r(chi2), r(df), r(p)
                qui estimates restore _Orig
                mat rownames `nextrow' = "`slab`count''"
                mat `matsuest' = nullmat(`matsuest') \ `nextrow'
            } // if `omit'==real("`basecat'")
            else if `omit'!=real("`basecat'") {
                qui {
                    mlogit `depvar' `rhscmdline' `wtis' ///
                        if `sample' == 1 & `depvar' != `omit', b(`basecat')
                    local eqnames = e(eqnames)
                    estimates store _Aux
                    suest _Aux _Orig
                    estimates store _suest
                }
                foreach name in `eqnames' {
                   * common needed for factor variables
                   qui test [_Orig_`name' = _Aux_`name'], cons accum common
                }
                if "`detail'"=="detail" {
                    suest
                    test
                }
                mat `nextrow' = r(chi2), r(df), r(p)
                qui estimates restore _Orig
                mat rownames `nextrow' = "`slab`count''"
                mat `matsuest' = nullmat(`matsuest') \ `nextrow'
            } // else if `omit'!=real("`basecat'")
            local ++count
        } // while real("`count'") <= `catsN'
        mat colnames `matsuest' = chi2 df "P>chi2"
        local cspecis "& %16s | %9.`Dst'f & %4.0f & %7.`Dpv'f &"
        matlist `matsuest', cspec(`cspecis') rspec(`rspecis')
        display _new as text ///
        "  Note: A significant test is evidence against Ho."
        qui estimates restore _Orig
        return matrix suest `matsuest'
    }  // if "`doiia'"=="no"

//  Small-Hsiao test of iia

    * adapted from Nick Winter's -smhsiao-
    if "`dosmhsiao'"=="yes" {
        display as res _new "Small-Hsiao tests of IIA assumption `n'"
        display _new as text ///
"  Ho: Odds(Outcome-J vs Outcome-K) are independent of other alternatives"

        * randomly divide sample
        qui gen `touse' = e(sample)
        qui gen `samp' = round(uniform(),1)+1 if `touse'
        local y `e(depvar)'
        qui gen `tempdepvar' = `depvar'
        local varlist "`tempdepvar' `rhscmdline'"
        * verify subsamples include all outcomes
        qui ta `tempdepvar' if `touse' & `samp'==1
        local cat1 `r(r)'
        qui ta `tempdepvar' if `touse' & `samp'==2
        local cat2 `r(r)'
        if `cat1'!=`catsN' | `cat2'!=`catsN' {
            display as error ///
"random draw yielded empty cells for outcomes; try different set seeed #"
            exit
        }

        local count = 1
        local countto = `catsN' - 1
        * if base used with, count to #cats - 1, else to # of cats
        if ("`base'"!="") local countto = `catsN'
        local rspecis "&|"
        while `count' <= `countto' { // eliminating categories one at a time
            local rownm : word `count' of `catnms'
            local rspecis "`rspecis'&"
            local basecatopt "b(`basecat')" // use original basecat
            * if count is base category, change base
            if (`count'==`basecat') { // current count is base
                if `basecat'!=1  { // if base is not 1, use 1 as base
                    local bval : word 1 of `catvals'
                    local basecatopt "b(`bval')"
                }
                else if `basecat'==1 { // otherwise use 1 as base
                    local bval : word 2 of `catvals'
                    local basecatopt "b(`bval')"
                }
            }
            * get values of outcome since they might have gaps
            qui tab `tempdepvar' if `touse', matrow(`Vals')
            * value to eliminate from current mlogit
            local elimval : word `count' of `catvals'
            local catsNM1 = `catsN' - 1
            local Yvals = ""
            local EYvals = ""
            local EYeqs = ""
            local i = 1
            * get list of values excluding elimval
            while `i' <= `catsN' {
                * value #i for y outcomes
                local Yval`i' = `Vals'[`i',1]
                local Yvals "`Yvals' `Yval`i''" // list of y values
                * if not eliminated value add to EYval list
                if `Yval`i'' != `elimval' {
                    local EYvals "`EYvals' `Yval`i''"
                    local EYeqs "`EYeqs' `i'"
                    local Ylab`i' `Yval`i''
                }
                local ++i
            }
            tempvar lnL denom
            tempname b0a b0b b0ab b1b
            * sample 2 w/o eliminated value
            local quiis qui
            `quiis' mlogit `varlist' ///
                if (`touse' & `samp'==2 & `tempdepvar'!=`elimval') ///
                `wtis',  `basecatopt' nolog
            mat `v' = e(V)
            _ms_omit_info `v' // count o. columns in v
            local dof = rowsof(`v') - r(k_omit) // v305
            local lnL_1 = e(ll)

            `quiis' mlogit `varlist' if `touse' & `samp'==1 `wtis', `basecatopt'
            mat `b0a' = e(b)
            `quiis' mlogit `varlist' if `touse' & `samp'==2 `wtis', `basecatopt'
            mat `b0b' = e(b)
            * Zhang & Hoffman equation 9
            mat `b0ab' = (0.70710678)*(`b0a') + (0.29289322)*(`b0b')

            * get XBs & assemble denominator
            qui gen double `denom' = 0 if `touse'
            local i 1
            * cycle through values (w/o eliminated one)
            while `i' <= (`catsN'-1) {
                local cury : word `i' of `EYvals'
                local cureq : word `i' of `EYeqs'
                tempvar xb`cureq'
                if `cury' != `basecat' & "`cury'" != "`tmpbcat'" {
                    matrix score double `xb`cureq'' = `b0ab' ///
                        if `touse', eq(`Ylab`cureq'')
                }
                else {
                    qui gen double `xb`cureq''=0  // because (exp(0)=1)
                }
                qui replace `denom' = `denom' + exp(`xb`cureq'') if `touse'
                local ++i
            }
            * Log likelihood using amalgamated coeff
            qui gen double `lnL' = . if `touse'
            local i 1
            while `i'<=`catsNM1' {
                local cureq : word `i' of `EYeqs'
                local cury : word `i' of `EYvals'
                qui replace `lnL' = ln(exp(`xb`cureq'')/(`denom')) ///
                    if `tempdepvar'==`cury' & `touse'
                local i=`i'+1
            }
            * LnL for obs in 2nd sample w/o eliminated obs
            sum `lnL' if `touse' & `samp'==2 & `tempdepvar'!=`elimval', ///
                meanonly
            local lnL_0 = r(sum)
            local SH = -2 * (`lnL_0' - `lnL_1')
            local p = chiprob(`dof',`SH')
            mat `rowres' = `lnL_0' , `lnL_1' , `SH' , `dof' , `p'
            mat rownames `rowres' = "`rownm'"
            mat `matsh' = nullmat(`matsh') \ `rowres'
            local ++count
        }
        matrix colna `matsh' = "lnL(full)" "lnL(omit)" chi2   df  "P>chi2"
        local cspecis ///
        "& %16s | %9.`Dst'f & %9.`Dst'f &  %9.`Dst'f & %4.0f & %7.`Dpv'f &"
        matlist `matsh', cspec(`cspecis') rspec(`rspecis')
        display _new as text ///
        "  Note: A significant test is evidence against Ho."
        qui estimates restore _Orig
        return matrix smhsiao `matsh'
    } // if "`dosmhsiao'"=="no"

//  Wald tests for combining categories

    if ("`docombine'"=="yes" |"`dolrcombine'"=="yes") & `catsN' == 2 {
        display as error "tests of combining categories requires at least 3 categories"
        local nocomb "yes"
    }
    else if "`docombine'"=="yes" {
        display as res _new "Wald tests for combining alternatives `n'"
        display
        display as text ///
        "  Ho: All coefficients except intercepts associated with a given pair"
        display as text ///
        "      of alternatives are 0 (i.e., alternatives can be combined)"
        tokenize "`catvals'"
        qui estimates restore _Orig
        local count1 = 1
        local rspecis "&|"
        while `count1' <= (`catsN'-1) {
            local count2 = `count1' + 1
            while `count2' <= `catsN' {
                local rspecis "`rspecis'&"
                if ("``count1''"=="`basecat'") qui test [``count2'']
                else if "``count2''"=="`basecat'" qui test [``count1'']
                else qui test [``count2''=``count1'']
                local s1 : word `count1' of `catnms8'
                local s2 : word `count2' of `catnms8'
                local rownm "`s1' & `s2'"
                if `issvy'==0 {
                    mat `nextrow' = r(chi2), r(df), r(p)
                    mat rownames `nextrow' = "`rownm'"
                    mat `matcomb' = nullmat(`matcomb') \ `nextrow'
                    mat colnames `matcomb' = chi2 df "P>chi2"
                    local cspecis "& %16s | %9.`Dst'f & %4.0f & %7.`Dpv'f &"
                }
                else {
                    local cspecis ///
                    "& %16s | %9.`Dst'f & %4.0f & %5.0f & %7.`Dpv'f &"
                    mat `nextrow' = r(F), r(df), r(df), r(p)
                    mat rownames `nextrow' = "`rownm'"
                    mat `matcomb' = nullmat(`matcomb') \ `nextrow'
                    mat colnames `matcomb' = F df df_r "P>F"
                }
                local ++count2
            }
            local ++count1
        }
        matlist `matcomb', cspec(`cspecis') rspec(`rspecis')
        qui estimates restore _Orig
        return matrix combine `matcomb'
    } // if "`docombine'"=="yes"

//  LR tests for combining categories

    if "`dolrcombine'"=="yes" & "`nocomb'"=="" {

        display as res _new "LR tests for combining alternatives `n'"
        display _new as text ///
        "  Ho: All coefficients except intercepts associated" ///
        " with a given pair"
        display as text ///
        "      of alternatives are 0 (i.e., alternatives can be collapsed)"
        tokenize "`catvals'"
        local Nconstraint = 1
        local rspecis "&|"
        while `Nconstraint' <= (`catsN'-1) {
            local Nbase = `Nconstraint'+1
            while `Nbase' <= `catsN' {
                local rspecis "`rspecis'&"
                local catnm : word `Nconstraint' of `catnms'
                local catnm : subinstr local catnm " " "_", all // 312
                constraint define 999 [`catnm']
                qui mlogit `depvar' `rhscmdline' `wtis' ///
                    if `sample' == 1, base(``Nbase'') constr(999)
                qui lrtest . _Orig
                matrix `nextrow' = r(chi2), r(df), r(p)
                local s1 : word `Nconstraint' of `catnms8'
                local s2 : word `Nbase' of `catnms8'
                local rownm "`s1' & `s2'"
                matrix rowna `nextrow' = "`rownm'"
                mat `matlrc' = nullmat(`matlrc') \ `nextrow'
                local ++Nbase
            }
            local ++Nconstraint
        } // while `Nconstraint' <= `catsN'
        mat colnames `matlrc' = chi2 df "P>chi2"
        local cspecis "& %25s | %9.`Dst'f & %4.0f & %7.`Dpv'f &"
        matlist `matlrc', cspec(`cspecis') rspec(`rspecis')
        qui estimates restore _Orig
        return matrix lrcomb `matlrc'
    } // if "`dolrcombine'"=="yes"

    qui estimates restore _Orig
    capture estimates drop _Orig
    capture estimates drop _OrigNewBase
    capture estimates drop _Aux
    capture estimates drop _suest
end
exit
* version 3.3.3 2014-06-19 | long freese | fix i# cleanup
* version 3.3.2 2014-05-30 | long freese | select fv with mlogtest, lr
* version 3.3.1 2014-03-25 | long freese | all fv
* version 3.3.0 2014-02-18 | long freese | spost13 release
* version 3.3.4 2015-10-20 | long freese | SH with 0 base
