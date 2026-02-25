*!  version 1.1.0 \ scott long 2016-06-19
*  version 1.0.0 \ scott long 2007-08-05
*       add names option

//  syntax:     nmlab <list of variables>, column(for labels) number vl
//  task:       list variable names and labels
//  project:    workflow chapter 4
//  author:     scott long \ 2007-08-05


capture program drop nmlab
capture program drop DisplayInCols




program define nmlab
    version 8, missing
    syntax [varlist] [, COLumn(integer 0) NUMber vl names ]
    tokenize `varlist'
    local stop : word count `varlist'

if "`names'"=="names" {
    di
    DisplayInCols `varlist'
    exit
}


    local len = 0
    local i 1
    while `i' <= `stop' {
        local l = length("``i''")
        if `l'>`len' local len = `l'
        local i = `i' + 1
    }
    if `column'==0 local column = `len' + 3

    display
    local i 1
    if "`number'"=="number" {
        local column = `column' + 6
    }
    else {
        local n ""
    }

    * value label location
    if "`vl'"=="vl" {
        local column2 = `column' + 11 // for labels
    }

    while `i' <= `stop' {
        local varlbl :  variable label ``i'' // grab var label
        local vallbl : value label ``i'' // grab value label

        if "`number'"=="number" {
            local n = substr(string(`i',"%4.0f") + ".     ",1,6)
        }
        if "`vl'"!="vl" {
            display in green "`n'``i''" in y _col(`column') "`varlbl'"
        }
        else { // show value label
            display in green "`n'``i''" in white _col(`column') ///
                "`vallbl'" in y _col(`column2') "`varlbl'"
        }
        local i = `i' + 1
    }
end

* txt `indent' `skip' 0 `vlist'
program DisplayInCols /* sty #indent #pad #wid <list>*/

    local indent 0
    local pad    2
    local wid 0  0

	local indent = cond(`indent'==. | `indent'<0, 0, `indent')
	local pad    = cond(`pad'==. | `pad'<1, 2, `pad')
	local wid    = cond(`wid'==. | `wid'<0, 0, `wid')

	local n : list sizeof 0
	if `n'==0 {
		exit
	}

	foreach x of local 0 {
		local wid = max(`wid', length(`"`x'"'))
	}

	local wid = `wid' + `pad'
	local cols = int((`c(linesize)'+1-`indent')/`wid')

	if `cols' < 2 {
		if `indent' {
			local col "_column(`=`indent'+1')"
		}
		foreach x of local 0 {
			di as `sty' `col' `"`x'"'
		}
		exit
	}

	local lines = `n'/`cols'
	local lines = int(cond(`lines'>int(`lines'), `lines'+1, `lines'))

    local ivar = 0
	forvalues i=1(1)`lines' {
		local col = `indent' + 1
		forvalues j=1/`cols' {
            local ++ivar
            local x : word `ivar' of `0'
			di in green _column(`col') "`x'" _c
			local col = `col' + `wid'
		}
		di
	}

end
