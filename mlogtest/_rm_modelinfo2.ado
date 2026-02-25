*! version 2.0.0 2017-08-23 | long freese | stata 15 parameter naming fix

//  sreturns with information about the last estimation command

program define _rm_modelinfo2, rclass

    version 11.2
    local cmdnm "`e(cmd)'"
    if "`cmdnm'"=="" {
        display as error "estimation results are not in memory"
        exit
    }

//  type of model

    if "`cmdnm'"=="logistic"  local cmdnm "logit" // 2015-02-10

    local cmdunclassified = 1
    if "`cmdnm'"=="logit" | "`cmdnm'"=="probit" | "`cmdnm'"=="cloglog" {
        local cmdbin = 1
        local cmdunclassified = 0
    }
    else local cmdbin = 0

    if "`cmdnm'"=="mlogit" | "`cmdnm'"=="mprobit" | "`cmdnm'"=="slogit" {
        local cmdnrm = 1
        local cmdunclassified = 0
    }
    else local cmdnrm = 0

    if "`cmdnm'"=="poisson"   | "`cmdnm'"=="nbreg" | "`cmdnm'"=="zip" ///
      | "`cmdnm'"=="zinb"     | "`cmdnm'"=="ztp"   | "`cmdnm'"=="ztnb" ///
      | "`cmdnm'"=="tpoisson" | "`cmdnm'"=="tnbreg" {
        local cmdcrm = 1
        local cmdunclassified = 0
    }
    else local cmdcrm = 0

    if ("`e(prefix)'"=="svy") local cmdsvy = 1
    else local cmdsvy = 0

//  basic model information

    local cmdn = e(N)

//  lhs variable

    local lhsnm "`e(depvar)'"
    if "`e(cmd)'"=="intreg" {
        local lhsnm2 : word 2 of `lhsnm'
        local lhsnm  : word 1 of `lhsnm'
    }

    local lhsbaseval = .
    if ("`e(cmd)'"=="mlogit") local lhsbaseval = e(baseout)
    if ("`e(cmd)'"=="mprobit") local lhsbaseval = e(i_base)

    * information on outcome categories
    local lhscatn = e(k_cat)
    if `cmdnrm' local lhscatn = e(k_out)
    if `cmdbin' local lhscatn = 2
    capture qui tab `lhsnm' if e(sample)
    if _rc!=0 local lhsvalues = . // if too many values for tab
    else local lhsvalues = r(r)
    if `lhscatn'==. {
        local lhscatn = 0
        local lhscatnms ""
        local lhscatvals ""
        local lhscatnms8 ""
    }
    else {
        qui sum `lhsnm' if e(sample)
        local lmin = r(min)
        _labels2names `lhsnm' if e(sample), indexfrom(`lmin')
        local lhscatnms `"`s(labels)'"' // names
        _labels2names `lhsnm' if e(sample), indexfrom(`lmin') nolabel
        local lhscatvals `"`s(labels)'"' // values
        local ln = s(n_cat)
        if `lhscatn'!=`ln' {
            di as error "_rm_modelinfo2 error: e(k_cat)!=s(n_cat)"
            exit
        }
        local lhscatnms8 ""
        foreach cnm in `lhscatnms' {
            local cnm = substr("`cnm'",1,8)
            local lhscatnms8 "`lhscatnms8' `"`cnm'"' "
        }
    }

//  rhs information

    * rhs variables in the command line
    local rhscmdline "`e(cmdline)'"
    * remove svy:
    local rhscmdline : list retokenize rhscmdline
    local rhscmdline : subinstr local rhscmdline "svy : " "", all
    local rhscmdline : subinstr local rhscmdline "svy: " "", all
    * remove command and lhs name
    gettoken tmp rhscmdline : rhscmdline
    gettoken tmp rhscmdline : rhscmdline
    local icomma = strpos("`rhscmdline'",",")
    local rhscmdrhs1 "`rhscmdline'"
    if `icomma'>0 {
        local icomma = `icomma' - 1
        local rhscmdrhs1 = substr("`rhscmdline'",1,`icomma')
    }
    return local rhs_cmdline "`rhscmdrhs1'" // legacy -- so it won't break pgms
    return local rhs_cmdnms "`rhscmdrhs1'"
    local rhs_cmdnmsN = wordcount("`rhscmdrhs1'")
    return local rhs_cmdnmsN "`rhs_cmdnmsN'"

    * rhs2 in inf() models
    local cmdrhs2 ""
    local icomma = strpos("`rhscmdline'",",")
    if `icomma'>0 {
        local ++icomma
        local cmdrest = substr("`rhscmdline'",`icomma',.)
        local cmdrest : subinstr local cmdrest "inf(" "inflate(", all
        local cmdrest : subinstr local cmdrest "infl(" "inflate(", all
        local cmdrest : subinstr local cmdrest "infla(" "inflate(", all
        local cmdrest : subinstr local cmdrest "inflat(" "inflate(", all
        local iinf = strpos("`cmdrest'","inflate(")
        if `iinf'>0 {
            local iinf = `iinf' + 8
            local cmdrest = substr("`cmdrest'",`iinf',.)
            local iinf = strpos("`cmdrest'",")")
            local --iinf
            local cmdrhs2 = substr("`cmdrest'",1,`iinf')
        }
    }
    local rhs2_cmdnmsN = wordcount("`cmdrhs2'")
    return local rhs2_cmdnmsN "`rhs2_cmdnmsN'"
    return local rhs2_cmdnms "`cmdrhs2'"

    * e(b) parameter names
    tempname varinfo betainfo betas

    * 2.0.0 2017-08-23 stata 15 parameter naming fix
    if _caller()>=15 {
        _ms_lf_info
        local betanms `r(varlist)' // removes auxillary parameters
    }
    * pre stata 15 with _cons
    else {
        matrix `betas' = e(b)
        local betanms : colnames `betas'
        local nbetanms = wordcount("`betanms'")
        * drop o._cons from ologit type models
        local betanms : subinstr local betanms "o._cons" " ", all
        * get rid of _cons
        local betanms : subinstr local betanms "_cons" " ", all
    }        
    * matrix with info on each beta
    local nbetanms : list sizeof betanms
    matrix `betainfo' = J(13,`nbetanms',.)
    matrix rownames `betainfo' = ///
        1isbin 2iscat 3Ncats 4catbase 5Nomitted 6isfv 7isint 8intself ///
        9intselfonly 10intother 11factor 12omit 13level
    matrix colnames `betainfo' = `betanms'

    * 2.0.0 2017-08-23 stata 15 parameter naming fix
    if _caller()>=15 version 15: _rm_rhsnames rhsnms rhsn rhsnms2 rhsn2
    else _rm_rhsnames rhsnms rhsn rhsnms2 rhsn2

    local tmp : list uniq betanms
    foreach bnm in `tmp' {
        local ipos = strpos("`bnm'","o.")
        if `ipos'==0 local betanmsbase "`betanmsbase'`bnm' "
    }
    local betanmsbase = trim("`betanmsbase'")

    foreach bnm in `betanmsbase' {
        local ipos = strpos("`bnm'","b.")
        if `ipos'==0 local betanmsnobase "`betanmsnobase'`bnm' "
    }
    local betanmsnobase = trim("`betanmsnobase'")

//  fv types

    * 2017-08-23 2.0.0 Stata 15 auxillary parameters renaming fix
    if _caller()>=15 version 15: get_fv_type , rhs_beta(`rhsnms')
    else get_fv_type , rhs_beta(`rhsnms')

    foreach varset in cat bin core typeintprod typefactor ///
            typevariable nfactor {
        local rhs_`varset' "`s(rhs_`varset')'"
        local rhs_`varset' = ltrim("`rhs_`varset''")
    }
    local nrhs_core = wordcount("`rhs_core'")
    matrix `varinfo' = J(6,`nrhs_core',.)
    matrix rownames `varinfo' = ///
        1bin 2cat 3typeint 4typefac 5typevar 6Ncats
    matrix colnames `varinfo' = `rhs_core'

    * for zip and zinb // not affect by stata 15 aux parameter renaming
    if "`rhsnms2'"!="" {

        * 2017-08-23 2.0.0 Stata 15 auxillary parameters renaming fix
        if _caller()>=15 version 15: get_fv_type , rhs_beta(`rhsnms2')
        else get_fv_type , rhs_beta(`rhsnms2')

        foreach varset in cat bin core typeintprod typefactor ///
                typevariable nfactor {
            local rhs2_`varset' "`s(rhs_`varset')'"
            local rhs2_`varset' = ltrim("`rhs2_`varset''")
        }
    }

    * is it binary?
    foreach nm in `s(rhs_bin)' {
        matrix `varinfo'[1,colnumb(`varinfo',"`nm'")] = 1
    }
    * is it categorical?
    foreach nm in `s(rhs_cat)' {
        matrix `varinfo'[2,colnumb(`varinfo',"`nm'")] = 1
    }
    * is it typeint?
    foreach nm in `s(rhs_typeintprod)' {
        matrix `varinfo'[3,colnumb(`varinfo',"`nm'")] = 1
    }
    * is it typefac?
    local i = 0
    foreach nm in `s(rhs_typefactor)' {
        matrix `varinfo'[4,colnumb(`varinfo',"`nm'")] = 1
        local ++i
        local nfac : word `i' of `rhs_nfactor'
        matrix `varinfo'[6,colnumb(`varinfo',"`nm'")] = `nfac'
    }
    * is it typevar?
    foreach nm in `s(rhs_typevariable)' {
        matrix `varinfo'[5,colnumb(`varinfo',"`nm'")] = 1
    }

//  returns

    return matrix varinfo = `varinfo'
    return local cmdnm "`cmdnm'"
    return local cmdbin "`cmdbin'"
    return local cmdnrm "`cmdnrm'"
    return local cmdcrm "`cmdcrm'"
    return local cmdunclassified "`cmdunclassified'" // unclassified
    return local cmdsvy "`cmdsvy'"
    return local cmdn "`cmdn'"

    return local lhsnm "`lhsnm'"
    return local lhsnm2 "`lhsnm2'" // for intreg
    return local lhsbaseval `lhsbaseval'
    return local lhsvalues "`lhsvalues'"
    return local lhscatn "`lhscatn'"
    return local lhscatnms "`lhscatnms'"
    return local lhscatnms8 "`lhscatnms8'"
    return local lhscatvals "`lhscatvals'"

    return local rhsnms "`rhsnms'"
    return local rhsn `rhsn'
    return local rhsnms2 "`rhsnms2'"
    return local rhsn2 `rhsn2'

    return local rhs_core "`rhs_core'"
    return local rhs_bin "`rhs_bin'"
    return local rhs_cat "`rhs_cat'"
    return local rhs_typeintprod "`rhs_typeintprod'"
    return local rhs_typefactor "`rhs_typefactor'"
    return local rhs_nfactor "`rhs_nfactor'"
    return local rhs_typevariable "`rhs_typevariable'"

    return local rhs2_core "`rhs2_core'"
    return local rhs2_bin "`rhs2_bin'"
    return local rhs2_cat "`rhs2_cat'"
    return local rhs2_typeintprod "`rhs2_typeintprod'"
    return local rhs2_typefactor "`rhs2_typefactor'"
    return local rhs2_nfactor "`rhs2_nfactor'"
    return local rhs2_typevariable "`rhs2_typevariable'"

    return local rhs_betanmsnobase "`betanmsnobase'"
    return local rhs_betanmsbase "`betanmsbase'"

end

* syntax: get_fv_type, rhsnms
program define get_fv_type, sclass

    version 11.2
    * currently rhs_beta is not used; instead e(b) is used
    syntax , rhs_beta(string)

    * 2017-08-23 2.0.0 Stata 15 auxillary parameters renaming fix
    if _caller()>=15 {
        _ms_lf_info
        local betanms `r(varlist)' // removes auxillary parameters
    }
    else {
        tempname betas
        matrix `betas' = e(b)
        local betanms : colnames `betas' // aux parms removed later
    }

    local listnms "typevariable typefactor typeintprod core bin cat"

    * loop through parameters and determine variable types
    foreach bnm in `betanms' {
        _ms_parse_parts `bnm'
        * interaction or product
        if inlist(r(type), "interaction", "product") {
            *                # of parts in terms
            forvalues i = 1/`r(k_names)' { // each part in x##y##z...
                local op "`r(op`i')'" // factor operator
                local nm "`r(name`i')'"
                local rhs_core "`rhs_core'`nm' "
                local rhs_typeintprod "`rhs_typeintprod'`nm' "
                * co for omitted
                if ("`op'"=="co") | ("`op'"=="c") /// continuous
                     local rhs_typevariable "`rhs_typevariable'`nm' "
                else local rhs_typefactor "`rhs_typefactor'`nm' "
            }
        }
        * c. variable not part of interaction
        else if r(type) == "variable" {
            local rhs_typevariable "`rhs_typevariable'`r(name)' "
            local rhs_core "`rhs_core'`r(name)' "
        }
        * i. factor not part of interaction
        else if r(type) == "factor" {
            local rhs_typefactor "`rhs_typefactor'`r(name)' "
            local rhs_core "`rhs_core'`r(name)' "
        }
        else {
            di as error "unexpected r(type). fatal error."
            exit
        }
    }

    * binary or categorical factor
    local rhs_typefactor : list uniq rhs_typefactor
    foreach nm in `rhs_typefactor' {
        fvexpand i.`nm'
        local fvlist `r(varlist)'
        local nfactor : list sizeof fvlist
        local rhs_nfactor "`rhs_nfactor'`nfactor' "
        if (`nfactor'==2) local rhs_bin "`rhs_bin'`nm' "
        else local rhs_cat "`rhs_cat'`nm' "
    }
    * cleanup and returns
    foreach lnm in `listnms' {
        local rhs_`lnm' : subinstr local rhs_`lnm' "_cons" " ", all
        local rhs_`lnm' : list uniq rhs_`lnm'
        sreturn local rhs_`lnm' "`rhs_`lnm''"
    }
    sreturn local rhs_nfactor "`rhs_nfactor'"

end
exit

* version 1.0.1 2014-03-03 | long freese | op type "co" and "c" evaluated
* version 1.0.0 2014-02-14 | long freese | spost13 release
* version 1.0.2 2015-02-10 | long freese | logistic treated as logit
* version 1.0.3V3 2017-08-02 | long freese | stata 15 _cut fix
