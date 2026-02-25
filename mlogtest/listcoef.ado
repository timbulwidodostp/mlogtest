*! version 5.0.1 2020-07-19 | long freese | zinb/zip when inf() is subset of x's
*! version 5.0.0 2017-08-23 | long freese | stata 15 parameter naming fix

//  list coefficients from last model

capture program drop listcoef
program define listcoef, rclass

    version 12 // requires 12 for _matrix_table
    tempname b v bvec b1 b12 b2 bnocon bnocon2
    tempname bS bSx bSy con2 conp2 conz2 conb conz conp
    tempname expb expbS expb2 expbS2 pctb pctb2 pctbS pctbS2
    tempname dft factr lnalpha p p2 sd sd2 sdb sdb2 sdx sdy sdyobs
    tempname sdystar se vare vx z zval zval2 slcat noomit pval
    tempname Mvec Mout M2vec M2out mnlMout

    local Mvecnms ""
    local M2vecnms ""
    * old formats for output
    local bF "%9.5f"
    local bzF "%7.3f"
    local bpF "%7.3f"
    local bSF "%7.4f"
    local bexpF "%8.4f"
    local bpctF "%8.1f"
    local bsdF "%8.4f"
    * newer formats with fewer decimals
    local bF "%9.4f"
    local bzF "%7.3f"
    local bpF "%7.3f"
    local bSF "%7.3f"
    local bexpF "%8.3f"
    local bpctF "%8.1f"
    local bsdF "%8.3f"

    syntax [varlist(default=none fv)] ///
        [,  PValue(real 0) Factor Percent Help ///
            CONStant constantoff NOConstant ///
            Odds Std Reverse gt lt ADJacent NOLabel EXpand ///
            Delta(real 1) vsquish NOFVLABEL ///
            POSitive NEGative ///
           ]

    local issvy = 0
    if "`e(prefix)'"=="svy" local issvy = 1

    local cmd "`e(cmd)'"
    if "`cmd'" == "" {
        di _n as error "listcoef must be run after a model is estimated"
        exit
    }
    if ("`cmd'"=="logistic") local cmd "logit"

    if "`positive'"!="" & "`negative'"!="" {
        di as erorr "options -positive- and -negative- cannot be used together"
        exit
    }
    if ("`positive'"!="" | "`negative'"!="") ///
        & ("`gt'"=="gt" | "`lt'"=="lt" | "`adjacent'"=="adjacent") {
        di as error ///
"-positive- and -negative- cannot be used with -gt-, -lt- or -adjacent-"
        exit
    }

    * check if constant in model
    model_has_cons #1
    if (`s(hascons)'!=1) local constantoff "constantoff"

    if _caller() > 11 {
        * with varlist, no constant unless constant option
        if ("`varlist'"!="" & "`constant'"!="constant") ///
                local constantoff "constantoff"
        if ("`constant'"=="") local constant constant
        if ("`constantoff'"=="constantoff") local constant ""
        if ("`noconstant'"=="noconstant") local constant ""
        if ("`cmd'"=="mlogit") local constant ""
        if ("`cmd'"=="slogit") local constant ""
    }
    if _caller() < 12 {
        di as error "note: annotation of factor variables requires Stata 12 or later"
        listcoef11
        exit
    }

//  expand variable list

    matrix `b' = e(b)
    local ebcolsN = colsof(`b')
    local ebnms : colnames `b'
    local keepnms ""
    foreach vnm in `varlist' {
    * remove i. ibn. 1. etc.
        local ipos = strpos("`vnm'",".")
        local vnm = substr("`vnm'",`ipos'+1,.)
        forvalues i = 1/`ebcolsN' {
            local ebnm : word `i' of `ebnms'
            _ms_parse_parts `ebnm'
            *2013-12-17
            if r(omit)!=1 { // drop if o.varname
                local name "`r(name)'"
                if ("`name'"=="`vnm'") local keepnms "`keepnms'`ebnm' "
                forvalues ipart = 1/9 {
                    local name "`r(name`ipart')'"
                    if ("`name'"=="`vnm'") local keepnms "`keepnms'`ebnm' "
                }
            }
        }
    }
    local varlist : list uniq keepnms

    local fvops = "`s(fvops)'"=="true" | _caller() >= 11
    if `fvops' {
        local eqns = 1 + cond("`iszero'"=="yes",1,0)
        forvalues i = 1/`eqns' {
            if (`eqns'==1) _ms_extract_varlist `varlist'
            else _ms_extract_varlist `varlist', noomit eq(#`i') nofatal
            local list `list' `r(varlist)'
        }
        local varlist `list'
    }

    local noLTpvalue = "yes"

*   coeftyp: which types of coefficients can be computed
*       bSx       beta x-std          bSx
*       bSy       beta y-std          bSy
*       bS        beta xy-std         bS
*       expb        exp(beta)           expb
*       expbS     exp(beta x-std)     expstd
*       pctb        %                   pct
*       pctbS     % xstd              pctstd
*       byopt       ?
*   modltyp: model class
*       tdist       use t not z for p-values
*       count       count model
*       zero        zero-inflated model
*       ystar       latent dependent variable
*       ystd:
*       special     own loop for calculating odds/std coefs
*       nosdx       do not report sdx
*       nocon       do not allow constant option
*   defhead: default header type: std odds count

    if "`cmd'" == "regress" {
        local coeftyp "bSx bSy bS"
        local modltyp "tdist ystd"
        local defhead "std"
    }
    if "`cmd'" == "logit" | "`cmd'" == "logistic" | "`cmd'" == "ologit" {
        local coeftyp "expb expbS byopt pctb pctbS bSx bSy bS"
        local modltyp "ystar ystd reverse"
        local defhead "odds"
    }
    if "`cmd'" == "probit" | "`cmd'" == "oprobit" {
        local coeftyp "bSx bSy bS"
        local modltyp "ystar ystd"
        local defhead "std"
    }
    if "`cmd'" == "cloglog" {
        local coeftyp "bSx"
        local modltyp ""
        local defhead "std"
    }
    if "`cmd'" == "mlogit" {
        local coeftyp "expb expbS byopt pctb pctbS"
        local modltyp "special nocon"
        local defhead "odds"
    }
    if "`cmd'" == "mprobit" {
        local coeftyp "bSx"
        local modltyp "special nocon"
        local defhead "std"
    }
    if "`cmd'" == "slogit" {
        local coeftyp "expb expbS byopt pctb pctbS"
        local modltyp "nocon"
        local defhead "odds"
        if e(k_dim) != 1 {
            di as err "listcoef only works for slogit ..., dim(1)"
            exit
        }
        matrix `slcat' = e(outcomes)
        local slncats = e(k_out)
        local c1 = `slcat'[1,1]
        foreach i of numlist 2/`slncats'  {
            local c1p1 = `c1' + 1
            local c2 = `slcat'[`i',1]
            if `c2'!=`c1p1' {
                di as err ///
            "listcoef with slogit does not allows skips in outcome categories"
                exit
            }
            local c1 = `c2'
        }
    }
    if "`cmd'" == "clogit" | "`cmd'" == "rologit" {
        local coeftyp "expb byopt pctb"
        local modltyp "nosdx nocon reverse"
        local defhead "odds"
    }
    if "`cmd'" == "tobit" | "`cmd'" == "cnreg" | "`cmd'" == "intreg" {
        local coeftyp "bSx bSy bS"
        local modltyp "tdist ystar ystd"
        local defhead "std"
    }
    if "`cmd'" == "poisson" | "`cmd'" == "nbreg" | "`cmd'" == "ztp" ///
        | "`cmd'" == "ztnb" | "`cmd'" == "tpoisson" | "`cmd'" == "tnbreg" {
        local coeftyp "expb expbS byopt pctb pctbS"
        local modltyp "count"
        local defhead "count"
    }
    if "`cmd'" == "zip" | "`cmd'" == "zinb" {
        local coeftyp "expb expbS byopt pctb pctbS"
        local modltyp "count zero"
        local defhead "count"
    }
    if "`coeftyp'" == "" {
        display as error "listcoef does not work with `e(cmd)'"
        exit
    }

//  unpack coeftyp

    local defopt = "yes"
    local ncoeftyps : word count `coeftyp'
    forvalues count = 1/`ncoeftyps' {
        local type : word `count' of `coeftyp'
        if ("`type'"=="bSx") local dobSx "`defopt'"
        if ("`type'"=="bSy") local dobSy "`defopt'"
        if ("`type'"=="bS") local dobS "`defopt'"
        if ("`type'"=="pctb") local dopctb "`defopt'"
        if ("`type'"=="pctbS") local dopctbS "`defopt'"
        if ("`type'"=="expb") local doexpb "`defopt'"
        if ("`type'"=="expbS") local doexpbS "`defopt'"
        if ("`type'"=="byopt") local defopt "option"
    }

//  parse options and check for errors

    local opterr "no"
    if "`std'" == "std" {
        if ("`dobSx'"=="option") local dobSx "yes"
        else if ("`dobSx'"=="") local opterr "std"
        if ("`dobSy'"=="option") local dobSy "yes"
        if ("`dobS'"=="option") local dobS "yes"
        if ("`cmd'"=="ologit") local defhead "std"
        if ("`cmd'"=="logit") local defhead "std" // 2014-07-29
        local doexpb ""
        local doexpbS ""
        local dopctb ""
        local dopctbS ""
    }

    if "`percent'"=="percent" & "`factor'"=="factor" {
        di as error "options percent and factor cannot both be used together"
        exit
    }
    if "`std'"=="std" & "`factor'"=="factor" {
        di as error "options std and factor cannot both be used together"
        exit
    }
    if "`std'" == "std" & "`percent'" == "percent" {
        di as error "options std and percent cannot both be used together"
        exit
    }
    if "`percent'" == "percent" {
        if "`dopctb'"=="option" {
            local dopctb "yes"
            if ("`doexpb'"=="yes") local doexpb ""
        }
        else if ("`dopctb'"=="") local opterr "percent"
        if "`dopctbS'"=="option" {
            local dopctbS "yes"
            if ("`doexpbS'"=="yes") local doexpbS ""
        }
    }
    if "`factor'" == "factor" {
        if "`doexpb'"=="option" {
            local doexpb "yes"
            if ("`dopctb'"=="yes") local dopctb ""
        }
        else if ("`doexpb'"=="") local opterr "odds"
        if "`doexpbS'"=="option" {
            local doebSx "yes"
            if ("`dopctbS'"=="yes") local dopctbS ""
        }
    }
    if "`opterr'" != "no" {
        di as error "option `opterr' not allowed after `cmd'"
        exit
    }

//  unpack mldtyp: define what to do for each type of model

    if (index("`modltyp'","tdist")==0) local t_or_z "z"
    else local t_or_z "t"
    if `issvy' local t_or_z "t" // svy is always t

    if (index("`modltyp'","count")==0) local iscount "no"
    else local iscount "yes"
    if (index("`modltyp'","zero")==0) local iszero "no"
    else local iszero "yes"
    if (index("`modltyp'","ystar")==0) local isystar "no"
    else local isystar "yes"
    if (index("`modltyp'","ystd")==0) local isystd "no"
    else local isystd "yes"
    if (index("`modltyp'","special")==0) local isspec "no"
    else local isspec "yes"
    if (index("`modltyp'","nosdx")==0) local isnosdx "no"
    else local isnosdx "yes"
    if (index("`modltyp'","nocon")==0) local iscon "yes"
    else local iscon "no"
    if (index("`modltyp'","reverse")==0) local isrev "no"
    else local isrev "yes"

    if "`constant'"!="" & "`iscon'"=="no" {
        di as error "constant option not allowed for this model"
        exit
    }
    if "`reverse'"!="" & "`isrev'"=="no" {
        di as error "reverse option not allowed for this model"
        exit
    }
    if ("`iszero'"=="yes") & ("`e(inflate)'"!="logit") {
       display as error ///
       "listcoef requires logit inflation for `cmd'"
        exit
    }

//  get model information

    local nobs = e(N)
    if "`t_or_z'"=="t" scalar `dft' = `nobs' - e(df_m) - 1
    if `issvy' scalar `dft' = e(N_psu) - e(N_strata)

* 5.0.0 2017-08-23 stata 15 parameter naming fix

    * 5.0.0 2017-08-23 stata 15 parameter naming fix
    if (_caller()>=15) version 15: _rm_modelinfo2
    else _rm_modelinfo2

    local rhsnms "`r(rhsnms)'"
    local nrhs "`r(rhsn)'"
    local rhsnmsbase "`r(rhs_betanmsbase)'"
    local rhsnms2 "`r(rhsnms2)'"
    local nrhs2 "`r(rhsn2)'"
    local lhsnam2 "`r(rhs2_cmdnms)'"
    local lhsnm "`r(lhsnm)'"
    local ncats `r(lhscatn)'

    if "`e(cmd)'"=="tobit"   | "`e(cmd)'"=="intreg"  | "`e(cmd)'"=="cnreg" ///
     | "`e(cmd)'"=="regress" | "`e(cmd)'"=="poisson" | "`e(cmd)'"=="nbreg" ///
     | "`e(cmd)'"=="ztp"     | "`e(cmd)'"=="ztnb"    | "`e(cmd)'"=="zip" ///
     | "`e(cmd)'"=="zinb"    | "`cmd'" == "tpoisson" | "`cmd'" == "tnbreg" ///
    {
        local ncats = 2
    }
    local lhscatnms "`r(lhscatnms)'"
    local lhscatvals "`r(lhscatvals)'"
    local refcat "`r(lhsbaseval)'"
    if "`cmd'"=="mprobit" & `ncats'==2 {
        display as error "use probit instead of mprobit with binary outputs"
        exit
    }
    if "`cmd'"=="logit" | "`e(cmd)'"=="logistic" | "`e(cmd)'"=="clogit" {
        local nmFrom : word 2 of `lhscatnms'
        local nmTo   : word 1 of `lhscatnms'
    }
    else if "`cmd'" == "ologit" {
        local nmFrom ">m"
        local nmTo   "<=m"
    }
    else if "`cmd'" == "rologit" {
        local nmFrom "ranked ahead"
        local nmTo   "ranked behind"
    }
    if "`reverse'"!="" {
        local temp   "`nmFrom'"
        local nmFrom "`nmTo'"
        local nmTo   "`temp'"
    }

    if "`cmd'" == "mlogit" | "`cmd'" == "mprobit" {
        local catnms "`lhscatnms'"
        local catnums `lhscatvals'
        if ("`nolabel'"=="nolabel") local catnms `catnums'
    }

    if "`e(wtype)'"=="iweight" {
        di _n as error "listcoef does not work with iweights"
        exit
    }
    local wtis ""
    if ("`e(wtype)'"!="") local wtis "[`e(wtype)'`e(wexp)']"
    if "`e(wtype)'"=="pweight" & !`issvy' {
        local wtis "[aweight`e(wexp)']"
        di "note: pweights not compatible with summarize;" ///
           " weights treated as aweights."
    }

    if "`e(cmd)'"=="mlogit" | "`e(cmd)'"=="slogit" | "`e(cmd)'"=="mprobit" {
        local cmdline `e(cmdline)'
        local isconstraint = regexm("`cmdline'","[(a-zA-Z)]* c[a-z]*\(")
        if `isconstraint'==1 {
            display as error ///
            "listcoef does not allow constrained estimation for this model"
            exit
        }
    }

    if "`e(cmd)'"=="mlogit" | "`e(cmd)'"=="mprobit" {
        _rm_mlogitbv `b' `v'
        matrix `bvec' = vec(`b'')'
    }
    else { // not mlogit or mprobit
        matrix `b' = e(b)
        matrix `v' = e(V)
        local coln : colfullnames `b'
        local cols = colsof(`b')
        _ms_omit_info `b' // names of o. columns
        matrix `noomit' = J(1,`cols',1) - r(omit)
        mata: noomit = st_matrix(st_local("noomit"))
        mata: newb = select(st_matrix(st_local("b")), noomit)
        mata: st_matrix(st_local("b"),newb)
        mata: newv = select(select(st_matrix(st_local("v")), noomit), noomit')
        mata: st_matrix(st_local("v"),newv)
        * reassign column names to b
        foreach var of local coln {
            _ms_parse_parts `var'
            if (!`r(omit)') local coln2 `coln2' `var'
        }
        matrix colnames `b' = `coln2'
        matrix colnames `v' = `coln2'
        matrix rownames `v' = `coln2'
        matrix `sdb' = vecdiag(`v')
        local nb = colsof(`b')
        matrix `bnocon' = `b'[1,1..`nrhs'] // trim _con
        matrix coleq `bnocon' = _
        * is there a constant/
        local bnms : colnames(e(b))
        local conpos : list posof "_cons" in bnms
        if "`iscon'"=="yes" & `conpos'!=0 {
            matrix `conb' = `b'[1,`nrhs'+1..`nrhs'+`ncats'-1]
            matrix `conz' = `conb'
            matrix `conp' = `conb'
            local i = 1
            while `i' < `ncats' {
                matrix `conz'[1,`i'] ///
                    = `conb'[1,`i'] / sqrt(`sdb'[1,`nrhs'+`i'])
                if "`t_or_z'"=="t" ///
                     matrix `conp'[1,`i'] = tprob(`dft',-abs(`conz'[1,`i']))
                else matrix `conp'[1,`i']  = 2*normprob(-abs(`conz'[1,`i']))
                local ++i
            }
        }
    } // not mlogit or mprobit

//  slogit coefficients

    if "`cmd'" == "slogit" {
        tempname slb slV slbeta slphi sltheta slcatnum slthetaV
        local slnvars = e(df_m) // # of rhs variables
        matrix `slb' = e(b)
        matrix `slV' = e(V)
        matrix `slb' = `b' // after omit
        matrix `slV' = `v' // after omit
        local slncats   = e(k_out)
        local slncatsm1 = e(k_out) - 1
        local slnphi   = e(k_out) - 1
        local slntheta = `slnphi'
        local slrefnum = e(i_base) // number of reference category
        local slrefnm `e(out`slrefnum')'
        matrix `slcatnum' = e(outcomes) // values for cats regardless of base
        matrix `slcatnum' = `slcatnum''
        local slrefrow = 0
        foreach i of numlist 1/`slncats'  {
            local cati = `slcatnum'[1,`i']
            if (`slrefnum'==`cati') local slrefrow = `i'
        }
        * if 1 is reference, from category is 2
        local slfromnum = 1
        if `slrefnum'==1 local slfromnum = 2
        local slfromnm `e(out`slfromnum')'
        matrix `slbeta' = `slb'[1,1..`slnvars']
        matrix `slphi' = `b' // after omit
        matrix `slphi' = `slphi'[1,`slnvars'+1..`slnvars'+`slncatsm1'],(0)
        matrix `sltheta' = `b' // after omit
        matrix `sltheta' = ///
          `sltheta'[1,`slnvars'+`slncatsm1'+1..`slnvars'+2*`slncatsm1'],(0)
        matrix `slthetaV' = ///
          `slV'[`slnvars'+`slncatsm1'+1..`slnvars'+2*`slncatsm1',.]
        matrix `slthetaV' = ///
           `slthetaV'[.,`slnvars'+`slncatsm1'+1..`slnvars'+2*`slncatsm1']
        * get theta# and phi#_# names
        local slphinm   : coleq `slphi'
        local slthetanm : coleq `sltheta'
        local nmFrom "`slrefnm'"
        * comparison category for beta's
        local nmTo "`slfromnm'"
    } // slogit

    * coefficients for 2nd equation for zip and zinb
    if "`iszero'"=="yes" {
        scalar `con2' = `b'[1,`nrhs'+2+`nrhs2']
        scalar `conz2' = `con2'/sqrt(`sdb'[1,`nrhs'+2+`nrhs2'])
        if "`t_or_z'"=="t" ///
             scalar `conp2' = tprob(`dft',-abs(`conz2'))
        else scalar `conp2'  = 2*normprob(-abs(`conz2'))
        matrix `bnocon2' = `b'[1,(`nrhs'+2)..(`nrhs'+`nrhs2'+1)]
        matrix coleq `bnocon2' = _
        matrix `sdb2' = `sdb'[1,(`nrhs'+2)..(`nrhs'+`nrhs2'+1)]
        * 5.0.0 2017-08-23 stata 15 parameter naming fix
        * 5.0.0 if (_caller()>=15) version 15: _rm_sum `wtis' if e(sample) == 1
        * 5.0.0 else _rm_sum `wtis' if e(sample) == 1
        * 5.0.1 2020-07-19
        if (_caller()>=15) version 15: _rm_sum `wtis' if e(sample) == 1, two
        else _rm_sum `wtis' if e(sample) == 1, two // for 2nd equation
        matrix `sd2' = r(matsd)
        matrix `sd2' = `sd2'[1,2...] // trim off lhs variable
    }

    * 5.0.0 2017-08-23 stata 15 parameter naming fix
    if (_caller()>=15) version 15: _rm_sum `wtis' if e(sample) == 1
    else _rm_sum `wtis' if e(sample) == 1

    matrix `sd' = r(matsd)
    scalar `sdy' = `sd'[1,1]
    scalar `sdyobs' = `sdy'
    matrix `sd' = `sd'[1,2...]  /* trim lhs variable */

//  parse varlist: names and #s of variables to print

    if "`varlist'" == "" {
        local prt_nms "`rhsnms'"
        forvalues count = 1/`nrhs' {
            local prt_nums "`prt_nums' `count'"
        }
        if "`iszero'"=="yes" {
            local prt_nmsinf "`rhsnms2'"
            forvalues count = 1/`nrhs2' {
                local prt_numinf "`prt_numinf' `count'"
            }
        }
    }

    * varlist specified, print output in varlist order
    else {
        _ms_extract_varlist `rhsnms', eq(#1)
        local rhsnms `r(varlist)'

        if "`iszero'"=="yes" {
            _ms_extract_varlist `rhsnms2', eq(#2)
            local rhsnms2 `r(varlist)'
        }
        local countto : word count `varlist'
        tokenize `varlist'
        forvalues count = 1/`countto' {
            local found = "no"
            forvalues count2 = 1/`nrhs' {
                local rhstmp : word `count2' of `rhsnms'
                if "``count''" == "`rhstmp'" {
                    local prt_nms "`prt_nms' `rhstmp'"
                    local prt_nums "`prt_nums' `count2'"
                    local found = "yes"
                }
            }
            if "`iszero'"=="yes" {
                forvalues count2 = 1/`nrhs2' {
                    local rhstmp : word `count2' of `rhsnms2'
                    if "``count''" == "`rhstmp'" {
                        local prt_nmsinf "`prt_nmsinf' `rhstmp'"
                        local prt_numinf "`prt_numinf' `count2'"
                        local found = "yes"
                    }
                }
            }
        }
    } // if a varlist has been specified

//  parse pvalue option

    local pcutoff = `pvalue'
    if (`pcutoff'>=1 & `pcutoff'<= 100) local pcutoff = `pcutoff' / 100
    if (`pcutoff'==0) local pcutoff = 1.00
    if `pcutoff' < 0 | `pcutoff' > 1 {
        di as error "pvalue() must be a nonzero probability"
        exit
    }

//  model is not special case

    if "`isspec'"=="no" {

        * compute sd(y*)
        if "`isystar'"=="yes" {
            * get cov(rhs) for computing var(y*)
            * nb: need all base names for case of i.catvar
            quietly matrix accum `vx' = `lhsnm' `rhsnmsbase' `wtis' ///
                if e(sample)==1 `in', deviations noconstant
            * removed omitted variables from cross product matrix
            _ms_omit_info `vx'
            local cols = colsof(`vx')
            matrix `noomit' = J(1,`cols',1) - r(omit)
            mata: noomit = st_matrix(st_local("noomit"))
            mata: newvx ///
                = select(select(st_matrix(st_local("vx")), noomit), noomit')
            mata: st_matrix(st_local("vx"),newvx)
            matrix `vx' = `vx'[2...,2...] // trim off lhs variable
            scalar `factr' = 1/(`nobs'-1) // 1 over nobs - 1
            matrix `vx' = `factr' * `vx'
            matrix `sdystar' = `bnocon' * `vx'
            matrix `sdystar' = `sdystar' * `bnocon''
            matrix `vare' = J(1,1,1) // default for probit
            scalar `factr' = 1
            if "`e(cmd)'"=="logit" | "`e(cmd)'"=="ologit" {
                scalar `factr' = (_pi*_pi)/3
            }
            if "`e(cmd)'"=="tobit" | "`e(cmd)'"=="intreg" | ///
                "`e(cmd)'"=="cnreg" {
                scalar `factr' = `b'[1,(`nrhs'+2)]^2
            }
            matrix `vare' = `factr' * `vare'
            matrix `sdystar' = `sdystar' + `vare'
            scalar `sdy' = sqrt(`sdystar'[1,1])
        }

        * compute standardized coefficients, t's, z's and p's
        local nx = colsof(`sd')
        matrix `bSy' = `bnocon'
        matrix `bSx' = `bnocon'
        matrix `bS'  = `bnocon'
        matrix `expb'  = `bnocon'
        matrix `expbS' = `bnocon'
        matrix `pctb'  = `bnocon'
        matrix `pctbS' = `bnocon'
        matrix `zval'  = `bnocon'
        matrix `p' = `bnocon'
        scalar `factr' = 1/`sdy'
        matrix `bSy' = `factr' * `bSy' // y-standardized betas

        forvalues i = 1/`nx' {
            matrix `sdb'[1,`i']  = sqrt(`sdb'[1,`i']) // sd of b's
            scalar `sdx' = `sd'[1,`i']
            if (`delta'!=1) scalar `sdx' = `delta' // if delta, not sd
            matrix `bSx'[1,`i']  = `bnocon'[1,`i']*`sdx' // b*sd_x
            matrix `bS'[1,`i']   = `bSy'[1,`i']*`sdx' // b*sd_x/sd_y)
            scalar `b1' = `b'[1,`i']
            matrix `expb'[1,`i'] = exp(`b1') // factor change
            matrix `expbS'[1,`i'] = exp(`b1'*`sdx') // std factor change
            matrix `pctb'[1,`i'] = (exp(`b1')-1)*100 // % change
            matrix `pctbS'[1,`i'] = (exp(`b1'*`sdx')-1)*100 // % std change
            if "`reverse'"!="" {
                matrix `expb'[1,`i'] = 1/exp(`b1')
                matrix `expbS'[1,`i'] = 1/exp(`b1'*`sdx')
                matrix `pctb'[1,`i'] = ((1/exp(`b1'))-1)*100
                matrix `pctbS'[1,`i'] = ((1/exp(`b1'*`sdx'))-1)*100
            }
            matrix `zval'[1,`i'] = `bnocon'[1,`i']/`sdb'[1,`i'] // t/z of b
            if "`t_or_z'"=="t" ///
                 matrix `p'[1,`i'] = tprob(`dft',-abs(`zval'[1,`i']))
            else matrix `p'[1,`i'] =  2*normprob(-abs(`zval'[1,`i']))
        }

        * coefficients for zip and zinb
        if "`iszero'"=="yes" {
            matrix `zval2' = `bnocon2'
            matrix `p2' = `bnocon2'
            matrix `expb2' = `bnocon2'
            matrix `expbS2' = `bnocon2'
            matrix `pctb2' = `bnocon2'
            matrix `pctbS2' = `bnocon2'
            local nx2 = colsof(`sd2')
            local i = 1
            while `i'<=`nx2' {
                matrix `sdb2'[1,`i'] = sqrt(`sdb2'[1,`i']) // sd of b's
                matrix `zval2'[1,`i'] = `bnocon2'[1,`i']/`sdb2'[1,`i']
                if "`t_or_z'"=="t" ///
                     matrix `p2'[1,`i'] = tprob(`dft',-abs(`zval2'[1,`i']))
                else matrix `p2'[1,`i']  = 2*normprob(-abs(`zval2'[1,`i']))
                matrix `expb2'[1,`i'] = exp(`bnocon2'[1,`i'])
                matrix `expbS2'[1,`i'] = exp(`bnocon2'[1,`i']*`sd2'[1,`i'])
                matrix `pctb2'[1,`i'] = (exp(`bnocon2'[1,`i'])-1)*100
                matrix `pctbS2'[1,`i'] ///
                    = (exp(`bnocon2'[1,`i']*`sd2'[1,`i'])-1)*100
                local ++i
            }
        }

//  print headers

        di _n as res "`cmd' (N=`nobs'): " _continue
        if "`defhead'"=="std" | "`std'"=="std" {
            di as res "Unstandardized and standardized estimates "
        }
        else {
            local header "Factor"
            if ("`percent'"=="percent") local header "Percentage"
            if "`defhead'"=="odds" | "`factor'"=="factor" {
                di as res "`header' change in odds " _continue
            }
            else if "`defhead'"=="count" {
                di as res "`header' change in expected count " _continue
            }
            if `pcutoff' < 1 & `pcutoff' >= 0 {
                di as res "(P<" %3.2f `pcutoff' ")" _continue
            }
            display
        }
        if "`cmd'"=="intreg" {
            di _n as res "    LHS vars: `e(depvar)'" _continue
        }
        if ("`defhead'"=="std" | "`std'"=="std") | ("`defhead'"=="count") {
            di as text _n "  Observed SD: " %7.4f `sdyobs'
            if ("`isystar'"=="yes") di "    Latent SD: " %7.4f `sdy'
        }
        if ("`cmd'"=="regress") di as text "  SD of error: " %7.4f e(rmse)
        if "`cmd'"=="tobit" | "`cmd'"=="cnreg" | "`cmd'"=="intreg" {
            local sde = `b'[1,`nrhs'+2]
            di as text "  SD of Error: " %7.4f `sde'
        }
        if "`defhead'"=="odds" | "`factor'"=="factor" {
            di _n as text "  Odds of: `nmFrom' vs `nmTo'"
        }
        if "`iszero'"=="yes" {
            local header "Factor"
            if ("`percent'"=="percent") local header "Percentage"
            di _n as text "Count equation: `header' change in expected " ///
                "count for those not always 0"
        }

//  header for columns

        local matcnms `"b `t_or_z' "P>|`t_or_z'|""'
        local fmtis "`bF' `bzF' `bpF'"
        if ("`dobSx'"=="yes") local fmtis "`fmtis' `bSF'"
        if ("`dobSy'"=="yes") local fmtis "`fmtis' `bSF'"
        if ("`dobS'"=="yes") local fmtis "`fmtis' `bSF'"
        if ("`doexpb'"=="yes") local fmtis "`fmtis' `bexpF'"
        if ("`doexpbS'"=="yes") local fmtis "`fmtis' `bexpF'"
        if ("`dopctb'"=="yes") local fmtis "`fmtis' `bpctF'"
        if ("`dopctbS'"=="yes") local fmtis "`fmtis' `bpctF'"

        if `delta'!=1 { // delta option
            if ("`dobSx'"=="yes") local matcnms `"`matcnms' bDeltaX"'
            if ("`dobSy'"=="yes") local matcnms `"`matcnms' bStdY"'
            if ("`dobS'"=="yes") local matcnms `"`matcnms' bDeltaStdY"'
            if ("`doexpb'"=="yes") local matcnms `"`matcnms' "e^b""'
            if ("`doexpbS'"=="yes") local matcnms `"`matcnms' "e^bDelta""'
            if ("`dopctb'"=="yes") local matcnms `"`matcnms' %"'
            if ("`dopctbS'"=="yes") local matcnms `"`matcnms' "%StdX""'
            if "`dobSy'"=="yes" { // std coef
                if "`isnosdx'"!="yes" {
                    local matcnms `"`matcnms' Delta"
                    local fmtis "`fmtis' `bSF'"
                }
            }
            else { // e(b)
                if "`isnosdx'"!="yes" {
                    local matcnms "`matcnms' Delta"
                    local fmtis "`fmtis' `bexpF'"
                }
            }
        }
        else{ // no delta option
            if ("`dobSx'"=="yes") local matcnms `"`matcnms' bStdX"'
            if ("`dobSy'"=="yes") local matcnms `"`matcnms' bStdY"'
            if ("`dobS'"=="yes") local matcnms `"`matcnms' bStdXY"'
            if ("`doexpb'"=="yes") local matcnms `"`matcnms' "e^b" "'
            if ("`doexpbS'"=="yes") local matcnms `"`matcnms' "e^bStdX""'
            if ("`dopctb'"=="yes") local matcnms `"`matcnms' "%""'
            if ("`dopctbS'"=="yes") local matcnms `"`matcnms' "%StdX""'
            if ("`isnosdx'"!="yes") local matcnms `"`matcnms' "SDofX""'
        }
        if ("`isnosdx'"=="no") local fmtis "`fmtis' `bsdF'"

//  print coefficients

        tokenize `prt_nums' // beta #'s to be printed
        local count = 1
        while "``count''" != "" {
            local indx: word `count' of `prt_nums'
            local vname: word `count' of `prt_nms'
            local vname2 "`vname'"
            local vname = abbrev("`vname'", 12)
            if `p'[1, `indx'] < `pcutoff' { // skip if not sig
                local noLTpvalue "no"
                matrix `Mvec' ///
                  = `bnocon'[1,`indx'], `zval'[1,`indx'], `p'[1,`indx']
                if "`dobSx'"=="yes"   matrix `Mvec' = `Mvec', `bSx'[1,`indx']
                if "`dobSy'"=="yes"   matrix `Mvec' = `Mvec', `bSy'[1,`indx']
                if "`dobS'"=="yes"    matrix `Mvec' = `Mvec', `bS'[1,`indx']
                if "`doexpb'"=="yes"  matrix `Mvec' = `Mvec', `expb'[1,`indx']
                if "`doexpbS'"=="yes" matrix `Mvec' = `Mvec', `expbS'[1,`indx']
                if "`dopctb'"=="yes"  matrix `Mvec' = `Mvec', `pctb'[1,`indx']
                if "`dopctbS'"=="yes" matrix `Mvec' = `Mvec', `pctbS'[1,`indx']
                if `delta'!=1 {  // delta option
                    if ("`isnosdx'"!="yes") matrix `Mvec' = `Mvec', `delta'
                }
                else {
                    if ("`isnosdx'"!="yes") {
                        matrix `Mvec' = `Mvec', `sd'[1,`indx']
                    }
                }
                local Mvecnms "`Mvecnms' `vname'"
                matrix `Mout' = nullmat(`Mout') \ `Mvec'
                local tabncols = colsof(`Mvec')
            } // if `p'[1, `indx']<`pcutoff'
            local count = `count' + 1
        } // loop through varlist

        if "`constant'"=="constant" & "`iscon'"=="yes" {
            if `ncats'==2 {
                matrix `Mvec' = J(1,`tabncols',.)
                matrix `Mvec'[1,1] = `conb'[1,1]
                matrix `Mvec'[1,2] = `conz'[1,1]
                matrix `Mvec'[1,3] = `conp'[1,1]
                local Mvecnms "`Mvecnms' constant"
                matrix `Mout' = nullmat(`Mout') \ `Mvec'
                if "`matrix'"!="" {
                    return scalar cons = `conb'[1,1]
                    return scalar cons_z = `conz'[1,1]
                    return scalar cons_p = `conp'[1,1]
                }
            }
            else {
                forvalues i = 1/`ncats' {
                    matrix `Mvec' = J(1,`tabncols',.)
                    matrix `Mvec'[1,1] = `conb'[1,`i']
                    matrix `Mvec'[1,2] = `conz'[1,`i']
                    matrix `Mvec'[1,3] = `conp'[1,`i']
                    local Mvecnms "`Mvecnms' constant`i'"
                    matrix `Mout' = nullmat(`Mout') \ `Mvec'
                }
            } // more than one constant
        }
        if "`cmd'"=="zinb" | "`cmd'"=="nbreg" | "`cmd'"=="ztnb" {
          if "`constant'"!="" {
            scalar `lnalpha' = `b'[1,`nb']
            matrix `Mvec' = J(1,`tabncols',.)
            quietly _diparm lnalpha, exp label("alpha") noprob
            matrix `Mvec'[1,1] = `lnalpha'
            local Mvecnms `"`Mvecnms' "alpha:lnalpha" "alpha:alpha""'
            matrix `Mout' = nullmat(`Mout') \ `Mvec'
            matrix `Mvec'[1,1] = r(est)
            matrix `Mout' = nullmat(`Mout') \ `Mvec'
          }
        }

        * phi and theta for slogit models
        if "`e(cmd)'" == "slogit" {
            matrix `b' = e(b)
            matrix `v' = e(V)
            * find location of 1st phi
            local bnms : colnames e(b)
            local icon : list posof "_cons" in bnms
            forvalues num = 1/`slnphi' {
                tempname phi Vphi sdphi zphi pphi
                scalar `phi' = `b'[1, `icon'] // grab phi
                local phinm : word `num' of `slphinm'
                scalar `Vphi' = `v'[`icon', `icon']
                scalar `sdphi' = sqrt(`Vphi')
                scalar `zphi' = `phi'/`sdphi'
                if "`t_or_z'"=="t" ///
                     scalar `pphi' = tprob(`dft',-abs(`zphi'))
                else scalar `pphi'  = 2*normprob(-abs(`zphi'))

                matrix `Mvec' = J(1,`tabncols',.)
                matrix `Mvec'[1,1] = `phi'
                matrix `Mvec'[1,2] = `zphi'
                matrix `Mvec'[1,3] = `pphi'
                local Mvecnms "`Mvecnms' phi:`phinm'"
                matrix `Mout' = nullmat(`Mout') \ `Mvec'
                local ++icon
            }
            forvalues num = 1/`slntheta' {
                tempname thetamat theta thetaVmat thetaV thetasd thetaz thetap
                local thetanm : word `num' of `slthetanm'
                scalar `theta' = `sltheta'[1,`num']
                scalar `thetaV' = `slthetaV'[`num',`num']
                scalar `thetasd' = sqrt(`thetaV')
                scalar `thetaz' = `theta' / `thetasd'
                if "`t_or_z'"=="t" ///
                     scalar `thetap' = tprob(`dft',-abs(`thetaz'))
                else scalar `thetap'  = 2*normprob(-abs(`thetaz'))
                matrix `Mvec' = J(1,`tabncols',.)
                matrix `Mvec'[1,1] = `theta'
                matrix `Mvec'[1,2] = `thetaz'
                matrix `Mvec'[1,3] = `thetap'
                local Mvecnms "`Mvecnms' theta:`thetanm'"
                matrix `Mout' = nullmat(`Mout') \ `Mvec'
            }
        } // end slogit - non expanded output

        display
        matrix colna `Mout' = `matcnms'
        matrix rowna `Mout' = `Mvecnms'
        _matrix_table `Mout', format(`fmtis') `vsquish' `nofvlabel'

        * LR test: code based on nbreg.ado version 3.3.9
        if "`e(cmd)'"=="nbreg" | "`cmd'" == "ztnb" {
          if "`constant'"!="" {

            if ((e(chi2_c)>0.005) & (e(chi2_c)<1e4)) | (ln(e(alpha))<-20) {
                local fmt "%-8.2f"
            }
            else local fmt "%-8.2e"
            scalar `pval' = chiprob(1, e(chi2_c))*0.5
            if (ln(e(alpha))<-20) scalar `pval'= 1
            di as text "  LR test of alpha=0: " ///
                `fmt' e(chi2_c) " Prob>=LRX2 = " %5.3f `pval'
          }
        }
        return matrix table = `Mout'

    } // if "`isspec'"=="no"

//  printing of non-standard models follows

//  mlogit: special model

* mlogit and mprobit: e(b) and e(V) have this structure:
*
*   varorder: v1 v2 ... _con
*   catorder: c1 c2 c3 ...
*   blockc#: c#1_v1 c#_v2 ... c#_con
*
*   bvec: blockc1 blockc2 blockc3  [1X(nvar+1)*(ncat-1)]
*
*   v = blockc1_X_blockc1
*       blockc2_X_blockc1 blockc2_X_blockc2
*       blockc3_X_blockc1 blockc3_X_blockc2 blockc3_X_blockc3

    if "`cmd'" == "mlogit" {
        display _new as res "`cmd' (N=`nobs'): " _continue
        local pcttext "Factor"
        if ("`percent'"=="percent") local pcttext = "Percentage"
        display as res "`pcttext' change in the odds of `lhsnm' " _c
        if `pcutoff' < 1 & `pcutoff' >= 0 {
            di as res "(P<" %3.2f `pcutoff' ")" _continue
        }
        display
        tokenize "`prt_nums'"

        * maximum length of labels
        local icat = 1
        local maxlen_catnm = 0
        while `icat' <= `ncats' {
            local catnum  : word `icat' of `catnums'
            local lcatnum = length("`catnum'")
            local catnm   : word `icat' of `catnms'
            local lcatnm = length("`catnm'")
            if (`lcatnum'>`maxlen_catnm') local maxlen_catnm = `lcatnum'
            if (`lcatnm'>`maxlen_catnm') local maxlen_catnm = `lcatnm'
            local ++icat
        }
        if (`maxlen_catnm'>15) local maxlen_catnm = 15 // to limit to 80 cols
        if (`maxlen_catnm'<11) local maxlen_catnm = 11

        * loop through variables
        local count = 1
        while "``count''" != "" {
          capture matrix drop `Mout'
          local ivar = "``count''"
          scalar `sdx' = `sd'[1,`ivar']
          local sdis = string(`sdx',"%8.3f")
          local vname: word `ivar' of `rhsnms'
          _ms_parse_parts `vname'
          if !`r(omit)' {
            di _n as res "Variable: `vname' (sd=`sdis')"
            local colnms `"b `t_or_z' "P>|`t_or_z'|""'
            local fmtis "`bF' `bzF' `bpF'"
            if "`doexpb'" == "yes" {
                local colnms `"`colnms' "e^b" "e^bStdX""'
                local fmtis "`fmtis' `bexpF' `bexpF'"
            }
            else { // pct
                local colnms `"`colnms' "%" "%StdX""'
                local fmtis "`fmtis' `bpctF' `bpctF'"
            }

            local b1cat 1 // 1...#cats but not necessarily the value of
                          // the category : see b1catnum for that
            local b2cat 0
            local b1block 0
            while `b1cat' <= `ncats' {
              local b1catnum  : word `b1cat' of `catnums'
              local b1catnm   : word `b1cat' of `catnms'
              * truncate long names
              local lnm = length("`b1catnm'")
              if `lnm'>`maxlen_catnm' {
                local b1catnm = substr("`b1catnm'",1,`maxlen_catnm')
              }
              * if not refcat, use next block of coefficients
              if (`b1catnum'!=`refcat') local ++b1block
              * the category # of the variable
              local b2block 0

              while `b2cat' < `ncats' {
                local ++b2cat
                local b2catnum  : word `b2cat' of `catnums'
                local b2catnm : word `b2cat' of `catnms'
                * truncate long names
                local lnm = length("`b2catnm'")
                if (`lnm'>`maxlen_catnm') {
                    local b2catnm = substr("`b2catnm'",1,`maxlen_catnm')
                }
                local c1 = substr(`"`b1catnm'           "',1,12)
                local c2 = substr(`"`b2catnm'           "',1,12)
                local rownm `"`c1' vs `c2'"'

                * if not refcat, use next block of coefficients
                if (`b2catnum'!=`refcat') local ++b2block
                local b1indx = `ivar' + ((`b1block'-1)*(`nrhs'+1))
                local b2indx = `ivar' + ((`b2block'-1)*(`nrhs'+1))
                * if both reference category, don't proccess
                if (`b1catnum'==`refcat') & (`b2catnum'==`refcat') {
                    continue
                }

                * get b for b1cat and b2cat
                if (`b1catnum'==`refcat') scalar `b1' = 0
                else scalar `b1' = `bvec'[1,`b1indx']
                if (`b2catnum'==`refcat') scalar `b2' = 0
                else scalar `b2' = `bvec'[1,`b2indx']
                scalar `b12' = `b1'-`b2'
                scalar `expb'  = exp(`b12')
                scalar `expbS' = exp(`b12'*`sdx')
                scalar `pctb'  = (`expb'-1)*100
                scalar `pctbS' = (`expbS'-1)*100

                * standard error of contrast
                if `b1catnum'!=`refcat' & `b2catnum'==`refcat' { // 2015-10-18
                    scalar `se' = sqrt(`v'[`b1indx',`b1indx'])
                }
                if `b1catnum'==`refcat' & `b2catnum'!=`refcat' {
                    scalar `se' = sqrt(`v'[`b2indx',`b2indx'])
                }
                if `b1catnum'!=`refcat' & `b2catnum'!=`refcat' {
                    scalar `se' = sqrt(`v'[`b1indx',`b1indx'] ///
                        + `v'[`b2indx',`b2indx'] - 2*`v'[`b1indx',`b2indx'])
                }
                scalar `z' = `b12'/`se'
                if "`t_or_z'"=="t" ///
                     scalar `p' = tprob(`dft',-abs(`z'))
                else scalar `p' =  2*normprob(-abs(`z'))

                if `p' <= `pcutoff' {
                  local noLTpvalue = "no"

                  if ("`gt'"!="gt" & "`lt'"!="lt") ///
                        | ("`gt'"=="gt") & (`b1catnum'>`b2catnum') ///
                        | ("`lt'"=="lt") & (`b1catnum'<`b2catnum') {
                    local absdif = abs(`b1catnum'-`b2catnum')
                    if ("`adjacent'"=="adjacent" & `absdif'==1) ///
                        | "`adjacent'"!="adjacent" {
                      * gt, lt, or adj
                      if "`positive'"=="" & "`negative'"=="" {
                        capture drop `Mvec'
                        matrix `Mvec' = `b12', `z', `p'
                        if "`doexpb'"=="yes" matrix `Mvec' = `Mvec', `expb'
                        else matrix `Mvec' = `Mvec', `pctb'
                        if "`doexpbS'"=="yes" mat `Mvec' = `Mvec', `expbS'
                        else matrix `Mvec' = `Mvec', `pctbS'
                        matrix rowna `Mvec' = "`rownm'"
                        matrix `Mout' = nullmat(`Mout') \ `Mvec'
                      } // not pos or neg
                    } // adjacent option
                  } // gt lt option
                  if "`positive'"=="positive" {
                    if `b12'>=0 {
                      capture drop `Mvec'
                      matrix `Mvec' = `b12', `z', `p'
                      if "`doexpb'"=="yes" matrix `Mvec' = `Mvec', `expb'
                      else matrix `Mvec' = `Mvec', `pctb'
                      if "`doexpbS'"=="yes" mat `Mvec' = `Mvec', `expbS'
                      else matrix `Mvec' = `Mvec', `pctbS'
                      matrix rowna `Mvec' = "`rownm'"
                      matrix `Mout' = nullmat(`Mout') \ `Mvec'
                    }
                  } // positive
                  if "`negative'"=="negative" {
                    if `b12'<=0 {
                      capture drop `Mvec'
                      matrix `Mvec' = `b12', `z', `p'
                      if "`doexpb'"=="yes" matrix `Mvec' = `Mvec', `expb'
                      else matrix `Mvec' = `Mvec', `pctb'
                      if "`doexpbS'"=="yes" mat `Mvec' = `Mvec', `expbS'
                      else matrix `Mvec' = `Mvec', `pctbS'
                      matrix rowna `Mvec' = "`rownm'"
                      matrix `Mout' = nullmat(`Mout') \ `Mvec'
                    }
                  } // negative
                } // `p' <= `pcutoff'
              } // while `b2cat' <= `ncats'
              local b2cat 0
              local ++b1cat
            } // b1cat
          }
          local ++count
          * print for current variable
          capture matrix list `Mout'
          if _rc==0 { // in cases no coefs selected with adj or pval optionts
            matrix colna `Mout' = `colnms'
            _matrix_table `Mout', format(`fmtis')  `vsquish' `nofvlabel'
            matrix roweq `Mout' = `vname'
            matrix `mnlMout' = nullmat(`mnlMout') \ `Mout'
          }
        } // while count

        return matrix table = `mnlMout' // all variables

    } // if "`cmd'" == "mlogit"

//  mprobit: special model

    if "`cmd'" == "mprobit" {
        di _n as res "`cmd' (N=`nobs'): " _c
        di as res "Unstandardized and standardized estimates for `lhsnm' " _c
        if `pcutoff' < 1 & `pcutoff' >= 0 {
            di as res "(P<" %3.2f `pcutoff' ")" _c
        }
        display
        tokenize "`prt_nums'"
        local count = 1
        while "``count''" != "" {
            local ivar = "``count''"
            scalar `sdx' = `sd'[1,`ivar']
          local sdis = string(`sdx',"%8.3f")

            local vname: word `ivar' of `rhsnms'
            _ms_parse_parts `vname'
            if !`r(omit)' {
                di _n as res "Variable: `vname' (sd=`sdis')"
                local colnms `"b `t_or_z' "P>|`t_or_z'|" bStdX"'
                local fmtis "`bF' `bzF' `bpF' `bSF'"
                local b1cat 1
                local b2cat 0
                local b1block 0
                while `b1cat' <= `ncats' {
                    local b1catnum  : word `b1cat' of `catnums'
                    local b1catnm   : word `b1cat' of `catnms'
                    * if not refcat, use next block of coefficients
                    if `b1catnum'!=`refcat' {
                        local ++b1block
                    }
                    local b2block 0
                    while `b2cat' < `ncats' {
                        local ++b2cat
                        local b2catnum  : word `b2cat' of `catnums'
                        local b2catnm : word `b2cat' of `catnms'
                        * if not refcat, use next block of coefficients* 3.1.4
                        if `b2catnum'!=`refcat' {
                            local ++b2block
                        }
                        local b1indx = `ivar' + ((`b1block'-1)*(`nrhs'+1))
                        local b2indx = `ivar' + ((`b2block'-1)*(`nrhs'+1))
                        * if both reference category, don't proccess
                        if (`b1catnum'==`refcat') & (`b2catnum'==`refcat') {
                            continue
                        }
                        local c1 = substr(`"`b1catnm'           "',1,12)
                        local c2 = substr(`"`b2catnm'           "',1,12)
                        local rownm `"`c1' vs `c2'"'
                        if (`b1catnum'==`refcat') scalar `b1' = 0
                        else scalar `b1' = `bvec'[1,`b1indx']
                        if (`b2catnum'==`refcat') scalar `b2' = 0
                        else scalar `b2' = `bvec'[1,`b2indx']
                        scalar `b12' = `b1'-`b2'
                        scalar `bSx' = `b12'*`sdx'
                        if `b1catnum'!=`refcat' & `b2catnum'==`refcat' {
                            scalar `se' = sqrt(`v'[`b1indx',`b1indx'])
                        }
                        if `b1catnum'==`refcat' & `b2catnum'!=`refcat' {
                            scalar `se' = sqrt(`v'[`b2indx',`b2indx'])
                        }
                        if `b1catnum'!=`refcat' & `b2catnum'!=`refcat' {
                            scalar `se' = sqrt(`v'[`b1indx',`b1indx'] ///
                                + `v'[`b2indx',`b2indx'] ///
                                - 2*`v'[`b1indx',`b2indx'])
                        }
                        scalar `z' = `b12'/`se'
                        if "`t_or_z'"=="t" ///
                             scalar `p' = tprob(`dft',-abs(`z'))
                        else scalar `p'  = 2*normprob(-abs(`z'))
                        if `p' <= `pcutoff' {
                            local noLTpvalue = "no"
                            if ("`gt'"!="gt" & "`lt'"!="lt") ///
                               | ("`gt'"=="gt") & (`b1catnum'>`b2catnum') ///
                               | ("`lt'"=="lt") & (`b1catnum'<`b2catnum') {
                                local absdif = abs(`b1catnum'-`b2catnum')
                                if ("`adjacent'"=="adjacent" & `absdif'==1) ///
                                    | "`adjacent'"!="adjacent" {
                                    capture drop `Mvec'
                                    matrix `Mvec' = `b12', `z', `p' , `bSx'
                                    matrix rowna `Mvec' = "`rownm'"
                                    matrix `Mout' = nullmat(`Mout') \ `Mvec'
                               } // adjacent option
                            } // gt lt option
                        } // `p' <= `pcutoff'
                    } // while `b2cat' <= `ncats'
                    local b2cat 0
                    local ++b1cat
                } // b1cat
            }
            local ++count
            matrix colna `Mout' = `colnms'
            _matrix_table `Mout', format(`fmtis') `vsquish' `nofvlabel'
            matrix roweq `Mout' = `vname'
            matrix `mnlMout' = nullmat(`mnlMout') \ `Mout'
            capture matrix drop `Mout'
        } // while count
        return matrix table = `mnlMout'
    } // if mprobit

//  special model--slogit // print expanded coefficients

    if "`cmd'" == "slogit" & "`expand'"=="expand" {
        di _n as res "`cmd' (N=`nobs'): " _c
        local pcttext "Factor"
        if "`percent'" == "percent" local pcttext = "Percentage"
        di as res "`pcttext' change in the odds of `lhsnm' " _c
        if `pcutoff' < 1 & `pcutoff' >= 0 {
            di as res "(P<" %3.2f `pcutoff' ")" _c
        }
        di
        local colnms `"b `t_or_z' "P>|`t_or_z'|""'
        local fmtis "`bF' `bzF' `bpF'"
        if "`doexpb'" == "yes" {
            local colnms `"`colnms' "e^b" "e^bStdX""'
            local fmtis "`fmtis' `bexpF' `bexpF'"
        }
        else { // pct
            local colnms `"`colnms' "%" "%StdX""'
            local fmtis "`fmtis' `bpctF' `bpctF'"
        }
        tokenize "`prt_nums'"

        local count = 1
        while "``count''" != "" { // loop through rhs variales
            local ivar = "``count''"
            scalar `sdx' = `sd'[1,`ivar']
            local sdis = string(`sdx',"%8.3f")
            local vname: word `ivar' of `rhsnms'
            * compute (phi_i - phi_j) * beta
            tempname phiBest phiBse estsl est se
            matrix `phiBest' = J(`slncats',`slncats',0)
            matrix `phiBse' = `phiBest'
            local i = 1
            * # cat - 1 since base is not in coef matrix
            foreach i of numlist 1(1)`slncatsm1' { // index for 1st phi
                * get phi # in case different base is used
                local phinm1 : word `i' of `slphinm' // name of ith parameter
                local phirow = substr("`phinm1'",6,.) //
                foreach j of numlist 1(1)`slncatsm1' { // index for 2nd phi
                    local phinm2 : word `j' of `slphinm'
                    local phicol = substr("`phinm2'",6,.)
                    * hold est from current model, leaving values in memory
                    _estimates hold `estsl', copy
                    qui nlcom ([`phinm1']_b[_cons] ///
                        - [`phinm2']_b[_cons])*_b[`vname'], post
                    matrix `est' = e(b)
                    matrix `se' = e(V)
                    _estimates unhold `estsl'
                    matrix `phiBest'[`phirow',`phicol'] = -1*`est'[1,1]
                    matrix `phiBse'[`phirow',`phicol'] = sqrt(`se'[1,1])
                }
                _estimates hold `estsl', copy
                qui nlcom ([`phinm1']_b[_cons]) * _b[`vname'], post
                matrix `est' = e(b)
                matrix `se' = e(V)
                _estimates unhold `estsl'
                matrix `phiBest'[`phirow',`slrefrow'] = -1*`est'[1,1]
                matrix `phiBse'[`phirow',`slrefrow'] = sqrt(`se'[1,1])
                matrix `phiBest'[`slrefrow',`phirow'] = `est'[1,1]
                matrix `phiBse'[`slrefrow',`phirow'] = sqrt(`se'[1,1])
            }
            di _n as res "Variable: `vname' (sd=`sdis')"

            * loop through all contrasts
            local i1 1
            local i2 1
            while `i1' <= `slncats' {
                local b1catnum = `slcatnum'[1,`i1']
                local b1catnm "`e(out`i1')'"
                while `i2' <= `slncats' {
                    scalar `b12' = `phiBest'[`i1',`i2']
                    scalar `se' = `phiBse'[`i1',`i2']
                    scalar `expb' = exp(`b12')
                    scalar `expbS' = exp(`b12'*`sdx')
                    scalar `pctb' = (`expb'-1)*100
                    scalar `pctbS' = (`expbS'-1)*100
                    scalar `z' = `b12'/`se'
                    if "`t_or_z'"=="t" ///
                         scalar `p' = tprob(`dft',-abs(`z'))
                    else scalar `p'  = 2*normprob(-abs(`z'))

                    if `p' <= `pcutoff' {
                        local noLTpvalue = "no"
                        * get outcome value of second outcome
                        local num2 = `slcatnum'[1,`i2']
                        local b2catnm "`e(out`i2')'"
                        if ("`gt'"!="gt" & "`lt'"!="lt") ///
                            | ("`gt'"=="gt") & (`b1catnum'>`num2') ///
                            | ("`lt'"=="lt") & (`b1catnum'<`num2') {
                            local absdif = abs(`b1catnum'-`num2')
                            if ("`adjacent'"=="adjacent" & `absdif'==1) ///
                                | "`adjacent'"!="adjacent" {

                    matrix `Mvec' = `b12', `z', `p'
                    local c1 = substr(`"`b1catnm'           "',1,12)
                    local c2 = substr(`"`b2catnm'           "',1,12)
                    local rownm `"`c1' vs `c2'"'
                    matrix rowna `Mvec' = "`rownm'"
                    if ("`doexpb'"=="yes") {
                        matrix `Mvec' = `Mvec', `expb', `expbS'
                    }
                    else matrix `Mvec' = `Mvec', `pctb', `pctbS'
                    matrix `Mout' = nullmat(`Mout') \ `Mvec'
                            } // adjacent option
                        } // gt lt option
                    } // `p' <= `pcutoff'
                    local ++i2
                    if `i1' == `i2' {
                        local ++i2
                    }
                } // while `i2' <= `slncats'
                local i2 1
                local ++i1
                } // i1

                matrix colna `Mout' = `colnms'
                _matrix_table `Mout', format(`fmtis')  `vsquish' `nofvlabel'
                matrix roweq `Mout' = `vname'
                matrix `mnlMout' = nullmat(`mnlMout') \ `Mout'
                capture matrix drop `Mout'
                local count = `count' + 1
        } // loop through variables

        return matrix table = `mnlMout'

    } // slogit expanded

//  help for first equation (zip and zinb printed below)

    if "`help'"=="help" {
        di "       b = raw coefficient"
        di "       `t_or_z' = `t_or_z'-score for test of b=0"
        di "   P>|`t_or_z'| = p-value for `t_or_z'-test"
        local std "standardized"
        if ("`dobSx'"=="yes") di "   bStdX = x-`std' coefficient"
        if ("`dobSy'"=="yes") di "   bStdY = y-`std' coefficient"
        if ("`dobS'"=="yes") di  "  bStdXY = fully `std' coefficient"
        if "`doexpb'" == "yes"  {
            if "`iscount'"=="yes" {
                di ///
"     e^b = exp(b) = factor change in expected count for unit increase in X"
            }
            else {
                di ///
"     e^b = exp(b) = factor change in odds for unit increase in X"
            }
        }
        if "`doexpbS'" == "yes" {
            if "`iscount'"=="yes" {
                di ///
" e^bStdX = exp(b*SD of X) = change in expected count for SD increase in X"
            }
            else {
                di ///
" e^bStdX = exp(b*SD of X) = change in odds for SD increase in X"
            }
        }
        if "`dopctb'" == "yes" {
            if "`iscount'"=="yes" {
                di ///
"       % = percent change in expected count for unit increase in X"
            }
            else {
                di ///
"       % = percent change in odds for unit increase in X"
            }
        }
        if "`dopctbS'" == "yes" {
            if "`iscount'"=="yes" {
                di ///
"   %StdX = percent change in expected count for SD increase in X"
            }
            else {
            di ///
"   %StdX = percent change in odds for SD increase in X"
            }
        }
        if ("`cmd'"!="mlogit") di in g "   SDofX = standard deviation of X"
    } // help

//  binary equation for zip zinb

    if "`iszero'"=="yes" {
        di _n as text "Binary equation: factor change in odds of always 0"
        local colnms `"b `t_or_z' "P>|`t_or_z'|""'
        local fmtis "`bF' `bzF' `bpF'"
        if "`doexpb'" == "yes" {
            local colnms `"`colnms' "e^b" "e^bStdX""'
            local fmtis "`fmtis' `bexpF' `bexpF'"
        }
        else { // pct
            local colnms `"`colnms' "%" "%StdX""'
            local fmtis "`fmtis' `bpctF' `bpctF'"
        }
        local colnms `"`colnms' SDofX"'
        local fmtis "`fmtis' `bsdF'"
        tokenize `prt_numinf'
        local count = 1
        while "``count''" != "" {
            local indx: word `count' of `prt_numinf'
            local vname: word `count' of `prt_nmsinf'
            local vname2 "`vname'"
            if `p2'[1, `indx'] < `pcutoff' {
                local noLTpvalue "no"
                matrix `M2vec' ///
                = `bnocon2'[1,`indx'], `zval2'[1,`indx'], `p2'[1,`indx']
                matrix rowna `M2vec' = `vname'
                if "`doexpb'" == "yes" {
                    matrix `M2vec' ///
                    = `M2vec', `expb2'[1,`indx'], `expbS2'[1,`indx']
                }
                if "`dopctb'" == "yes" {
                    matrix `M2vec' ///
                    = `M2vec', `pctb2'[1,`indx'], `pctbS2'[1,`indx']
                }
                matrix `M2vec' = `M2vec', `sd2'[1,`indx']
                matrix colna `M2vec' = `colnms'
                matrix `M2out' = nullmat(`M2out') \ `M2vec'
                local M2outcolsN = colsof(`M2out')
            } // if `p'[1, `count'] < `pcutoff'

            local ++count

        } // while count

        if "`constant'"=="constant" {
            matrix `M2vec' = J(1,`M2outcolsN',.)
            matrix `M2vec'[1,1] = `con2'
            matrix `M2vec'[1,2] = `conz2'
            matrix `M2vec'[1,3] = `conp2'
            matrix rowna `M2vec' = "constant"
            matrix `M2out' = nullmat(`M2out') \ `M2vec'
        }
        _matrix_table `M2out', format(`fmtis')  `vsquish' `nofvlabel'
        return matrix table2 = `M2out'

        //  vuong

        if e(vuong)~=. {
            local favor ""
            if e(vuong) > 0 {
                if "`t_or_z'"=="t" ///
                     local p = tprob(`dft',-e(vuong))
                else local p    = normprob(-e(vuong))

                if e(cmd)=="zip" {
                    if (`p'<.10) local favor "favoring ZIP over PRM."
                }
                else {
                    if (`p'<.10) local favor "favoring ZINB over NBRM."
                }
            }
            else {
                if "`t_or_z'"=="t" ///
                     local p = tprob(`dft',e(vuong))
                else local p    = normprob(e(vuong))
                if e(cmd)=="zip" {
                    if (`p'<.10) local favor "favoring PRM over ZIP."
                }
                else {
                    if (`p'<.10) local favor "favoring NBRM over ZINB."
                }
            }
            di "Vuong Test =" %6.2f e(vuong) " (p=" %4.3f `p' ") `favor'"
        }

        if "`help'"=="help" {
            di "       b = raw coefficient"
            di "       `t_or_z' = `t_or_z'-score for test of b=0"
            di "   P>|`t_or_z'  | = p-value for z-test"
            if "`doexpb'" == "yes"  {
                di ///
"     e^b = exp(b) = factor change in odds for unit increase in X"
                di ///
" e^bStdX = exp(b*SD of X) = change in odds for SD increase in X"
            }
            if "`dopctb'" == "yes" {
                di ///
"       % = percent change in odds for unit increase in X"
                di ///
"   %StdX = percent change in odds for SD increase in X"
            }
            di ///
"   SDofX = standard deviation of X"
        }
    } /* is -zip- or -zinb- */

    if "`noLTpvalue'" == "yes" {
        di _n as res "No results in which p < " %3.2f `pcutoff'
        exit
    }

    return local cmd `cmd'
    return scalar pvalue = `pcutoff'

end

program define model_has_cons, sclass // based on code from Jeff Pitblado
         args eqspec
         capture local cons = [`eqspec']_b[_cons]
         local hascons = c(rc) == 0
         if `hascons' {
                 if [`eqspec']_se[_cons] == 0 {
                         local hascons 0
                 }
         }
         sreturn local hascons = `hascons'
end

exit
* version 4.1.1 2014-05-03 | long freese | nocons fix; pos and neg
* version 4.1.0 2014-02-15 | long freese | spost13 release
* version 4.1.2 2014-07-29 | long freese | std label logistic; label if std
* version 4.1.3 2014-10-22 | long freese | improved svy
* version 4.2.0 2014-10-23 | long freese | t z labels; p values for t
* version 4.2.1 2015-10-18 mlogit refcat fix with base 0
* version 4.2.2 2017-08-02 stata 15 ologit fix
