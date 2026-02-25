*! version 2.0.0 2019-04-02 | long freese | adjust for non-estimable pw comparisons
*! version 1.0.0 2014-02-14 | long freese | spost13 release

//  return names of pw comparisons from RHS variable of factor type

//  Note: If perfect prediction, r(table_vs) has missing for some contrasts,
//  while r(table) has only estimable contrasts

//  all: all pw comparisons
//  est: only estimable comparisons

program define _rm_pwnames, sclass

    syntax , VARnm(string) [ NOISily]
    qui `noisily' pwcompare `varnm'
    tempname PWmat PRmat
    matrix `PWmat' = r(table_vs) // all possible pw comparisons
    local pwN = colsof(`PWmat')
    matrix `PRmat' = r(table) // probabilities
    local prN colsof(`PRmat')
    local PWcolnms : colnames(`PWmat')

    local iest = 0
    forvalues iall = 1/`pwN' { // loop through all possible pairwise comparisons

        if !missing(`PWmat'[1,`iall']) { // only process estimable pwcomparisons
            local ++iest
            _ms_element_info, element(`iall') matrix(`PWmat') compare
            * est pw name from value labels
            sreturn local txtlabel`iest' `"`r(level)'"'       
            * est pw name using category #'s
            local pwvalues : word `iest' of `PWcolnms'       
            local pwvalues : subinstr local pwvalues ".`varnm'" "", all
            local pwvalues : subinstr local pwvalues "bn" "", all
            sreturn local numlabel`iest' "`pwvalues'"        
            local tmp : subinstr local pwvalues "vs" " ", all
            * To and From numbers
            local pwFrom : word 2 of `tmp'
            sreturn local from`iest' = `pwFrom'        
            local pwTo : word 1 of `tmp'
            sreturn local to`iest' = `pwTo'        
        } // if estimable

    } // iall loop
    sreturn local npw = `pwN'
    sreturn local nvarcats = `prN'
end
exit

