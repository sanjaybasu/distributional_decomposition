
* y x


cap program drop figure2
program define figure2
syntax [if] [in], y(varlist) x(varlist max=1) [ xline(numlist) yline(numlist) ci(numlist) legend ///
    yrange(string asis) xrange(string asis) ytitle(string asis) xtitle(string asis) /// 
    title(string asis) line_pattern(string asis) line_width(string asis) save(string asis) ]

    *di "x/yrange: min(incr)max"
* titles later
    
    local text_size medsmall

    if "`line_pattern'" == "" {
        local line_pattern longdash solid shortdash_dot_dot shortdash_dot dash dash_dot shortdash  dot longdash_dot   
    }
    if "`line_width'" == "" {
        local line_width medthin medthin medthin medthin medthin medthin medthin
    }
    

    local yscale yscale(lcolor(black) lpattern(solid) line extend fextend) 
    local xscale xscale(lcolor(black) lpattern(solid) line extend fextend) 

    local title title("`title'", size(`text_size') color(black))
    * Axis label paramters
    local ylab_parm labels labcolor(black) angle(horizontal) labsize(`text_size') nogrid
    *format(%4.3f)
    * grid glwidth(thin)  glcolor(gs14) glpattern(shortdash)
    local ytitle ytitle(`ytitle', angle(horizontal) size(`text_size') margin(small))

    local xlab_parm labcolor(black) angle(horizontal) labsize(`text_size')
    * lcolor(red)
    local xtitle xtitle(`xtitle', size(`text_size') margin(small))

    local graphregion graphregion(margin(small) fcolor(white) lcolor(white) lpattern(solid) ifcolor(none) ilcolor(none) ilwidth(none) ilpattern(solid))
    local plotregion plotregion(margin(zero) fcolor(white) lcolor(white) lpattern(solid) ifcolor(white) ilcolor(white) ilpattern(solid))

    local i = 0
    foreach var in `y' {

        local line_pat `: word `=`i'+1' of `line_pattern'' 
        local line_wid `: word `=`i'+1' of `line_width'' 
        local linelist `linelist' (line `var' `x', lcolor(black) lpattern(`line_pat') lwidth(`line_wid'))

        if "`ci'" != "" {
            local stderr = subinstr("`var'","_b","_se",1)


            qui {
            cap drop `stderr'_h `stderr'_l
            gen double `stderr'_h = `var'+(invttail(`=`ci'-1',0.025))*`stderr'
            gen double `stderr'_l = `var'-(invttail(`=`ci'-1',0.025))*`stderr'
            }
            local cilist `cilist' (rspike `stderr'_h `stderr'_l `x', lpattern(solid) lcolor(gs10) lwidth(thin))
            *local cilist `cilist' (rarea `stderr'_h `stderr'_l `x', fcolor(gs12) fintensity(100) lcolor(gs12))
            local i = `i'+1
        }


        local i = `i'+1

        local label : variable label `var'
        if "`label'" != "" {
            local orderlist `orderlist' `i' "`label'"
        }
        else {
            local orderlist `orderlist' `i' "`var'"
        }

    }

    foreach i in `xline' {
        local xlineset `xlineset' xline(`i', lcolor(red) lpattern(solid) lwidth(thin))
    }

    foreach i in `yline' {
        local ylineset `ylineset' yline(`i', lcolor(black) lpattern(solid) lwidth(thin))
    }


    if "`legend'" != "" {
         * legend(order(1 "dd" )
        local legend legend(on order(`orderlist') nostack cols(1) size(`text_size') keygap(small) colgap(small) symxsize(*.8) bmargin(tiny) color(none) nobox position(7) region(fcolor(none) margin(tiny) lcolor(none) lwidth(none) lpattern(blank))) 
        * position(2) ring(0)
        * 
    }
    else {
        local legend legend(off)
    }

    twoway ///
    `cilist' `linelist' `if' `in', ///
    `yscale' `ylineset' ylabel(`yrange', `ylab_parm') ymtick(none) ///
    `xscale' `xlineset' xlabel(`xrange', `xlab_parm') xmtick(none) ///
    `xtitle' `ytitle' ///
    `title' ///
    `graphregion' `plotregion' `legend' xsize(6) ysize(4) 


    if "`save'" != "" {
        di "`save'"
        graph save `save', replace
    }
    else {
    }


end
