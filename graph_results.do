**************** probably want to be able to state where the legend goes


cap prog drop EraseList
cap prog define EraseList
syntax anything

    foreach i in `anything' {
        erase `i'
    }

end



cap prog drop GraphResults
cap prog define GraphResults

    syntax [anything] [, depvar(string asis) group(string asis) g1(string asis) g2(string asis)  ///
    xline(string asis) result_file(string asis) save(string asis)]

    di _n ""

    * Read result_file
    if ("`result_file'" == "") {
        if r(file) == "" {
            di "Last estimation results not found.  Specify result_file() or run"
            di "GraphResults following decomposition"
            exit
        }
        else {
            use "`=r(file)'", clear
        }
    }
    else {
        use "`result_file'", clear
    }

    di "GraphResults"
    di "Results from `result_file'"_n

    * Read meta-data
    local decomp_level = meta_val[1]
    local g1v = meta_val[4]
    local g2v = meta_val[5]
    local f_g1 = `g1v'`: di _dup(`decomp_level') "`g1v'"'
    local f_g2 = `g2v'`: di _dup(`decomp_level') "`g2v'"'
    local f_gc = `g1v'`: di _dup(`decomp_level') "`g2v'"'
    if regexm("`=meta_val[2]'", "(reps\()([0-9]+)(\))") == 1 local bs_reps = regexs(2)
    local N = `=meta_val[6]'+`=meta_val[7]'


    * Label decomposition results and counterfactual KDE/QFE
    foreach F in `=meta_val[3]' {
        label variable `F'_difference_b "Difference"
        label variable `F'_residual_b "Residual"

        label variable `F'`f_g1' "`g1'"
        label variable `F'`f_g2' "`g2'"
        label variable `F'`f_gc' "What to call this"

        local i = 1
        local list `anything'
        while regexm("`list'", "(\()([A-Za-z0-9_ ]*)(\))(.*)")==1 {
            local E = char(`i'+64)
            local `E' = regexs(2)
            local list = regexs(4)
            local list `list' `=regexs(2)'
            local i = `i' + 1
            *di "Effect `E' (``E'')"
            label var `F'_effect_`E'_b "``E''"
            label var `F'`f_g2'_effect_`E' "``E''"
        }
    }

    * Label axis
    local depvar "`depvar'"
    label variable x "`depvar'"
    label variable p "Percentile"

    * Ranges for graphs
    /*
    qui sum Q`f_g1'
    local min = round(r(min))
    local max = round(r(max))
    qui sum Q`f_g2'
    if r(min) > `min' local min = round(r(min))
    if r(max) < `max' local max = round(r(max))

    local ll = round(r(min)-sqrt(r(max)-r(min)))
    local ul = round(r(max)+sqrt(r(max)-r(min)))
    local ii = round((`ul'-`ll')/10)
    *******yrange(`ll'(`ii')`ul')
    */


    set graphics off
    di "Graphing..."


    ****** Observed KDE/QFE
qui{
    * Density
    figure2, y(f`f_g1' f`f_g2') x(x) title("Observed densities") xline(`xline') legend ytitle(Estimated Probability Density) save(densities_raw)
    figure2, y(f`f_g1' f`f_g2' f`f_gc') x(x) xline(`xline') legend ytitle(Probability Density) save(densities_g1g2gc)
    figure2, y(f`f_g2' f`f_gc') x(x) xline(`xline') legend ytitle(Probability Density) line_pattern(solid shortdash_dot_dot) save(densities_g2gc)

    * QF
    figure2 if (p>=0.05) & (p<=0.95), y(Q`f_g1' Q`f_g2') x(p) title(Estimated Quantile Function) legend ytitle(`depvar') xrange(0(0.1)1) save(qf_raw)
    figure2 if (p>=0.05) & (p<=0.95), y(Q`f_g1' Q`f_g2' Q`f_gc') x(p) legend ytitle(`depvar') xrange(0(0.1)1) save(qf_g1g2gc)
    figure2 if (p>=0.05) & (p<=0.95), y(Q`f_g2' Q`f_gc') x(p) legend ytitle(`depvar') xrange(0(0.1)1) line_pattern(solid shortdash_dot_dot) save(qf_g2gc)

    ****** Counterfactual KDE/QFE

    * Density
    local list
    forvalues i = 1(1)`decomp_level' {
        local list `list' densities_effect`i'.gph
        local E = char(`i'+64)
        figure2, y(f`f_g2' f`f_g2'_effect_`E') x(x) xline(`xline') legend ytitle(Probability Density) line_pattern(solid shortdash_dot_dot) save(densities_effect`i')
    }
    local list densities_g1g2gc.gph `list' densities_g2gc.gph
    qui graph combine `list', col(2) altshrink holes(2) ycommon title("Density functions", size(small)) xsize(8) ysize(11) scheme(s1mono)
    qui graph save densities, replace
    EraseList `list'

    * QF
    local list
    forvalues i = 1(1)`decomp_level' {
        local list `list' qf_effect`i'.gph
        local E = char(`i'+64)
        figure2 if (p>=0.05) & (p<=0.95), y(Q`f_g2' Q`f_g2'_effect_`E') x(p) legend ytitle(`depvar') xrange(0(0.1)1) line_pattern(solid shortdash_dot_dot) save(qf_effect`i')
    }
    local list qf_g1g2gc.gph `list' qf_g2gc.gph
    qui graph combine `list', col(2) altshrink holes(2) ycommon title("Quantile functions", size(small)) xsize(8) ysize(11) scheme(s1mono)
    qui graph save qf, replace
    EraseList `list'


    ****** Decomposition

    * Density
    figure2, y(f_difference_b) x(x) xline(`xline') yline(0) legend ytitle(Difference in densities) ci(`bs_reps') line_pattern(solid) save(diff_densities_raw)
    figure2, y(f_residual_b) x(x) xline(`xline') yline(0) legend ytitle(Difference in densities) ci(`bs_reps') line_pattern(solid) save(diff_densities_resid)

    local list
    forvalues i = 1(1)`decomp_level' {
        local E = char(`i'+64)
        figure2, y(f_effect_`E'_b) x(x) xline(`xline') yline(0) legend ytitle(Difference in densities) ci(`bs_reps') line_pattern(solid) save(diff_densities_`E')
        local list `list' diff_densities_`E'.gph
    }

    local list diff_densities_raw.gph `list' diff_densities_resid.gph
    qui graph combine `list', col(2) hole(2) altshrink ycommon title("Decomposition of density functions", size(medium)) xsize(8) ysize(11) scheme(s1mono)
    qui graph save diff_densities, replace
    EraseList `list'

    * QF
    figure2 if (p>=0.05) & (p<=0.95), y(Q_difference_b) x(p) xrange(0(0.1)1) legend ytitle(Difference in `depvar') xtitle(Percentile) ci(`bs_reps') line_pattern(solid) save(diff_qf_raw)
    figure2 if (p>=0.05) & (p<=0.95), y(Q_residual_b) x(p) xrange(0(0.1)1) legend ytitle(Difference in `depvar') xtitle(Percentile) ci(`bs_reps') line_pattern(solid) save(diff_qf_resid)

    local list
    forvalues i = 1(1)`decomp_level' {
        local E = char(`i'+64)
        figure2 if (p>=0.05) & (p<=0.95), y(Q_effect_`E'_b) x(p) xrange(0(0.1)1) legend ytitle(Difference in `depvar') xtitle(Percentile) ci(`bs_reps') line_pattern(solid) save(diff_qf_`E')
        local list `list' diff_qf_`E'.gph
    }
    local list diff_qf_raw.gph `list' diff_qf_resid.gph
    qui graph combine `list', col(2) hole(2) altshrink title("Decomposition of quantile functions", size(medium)) xsize(8) ysize(11) scheme(s1mono)
    * ycommon
    qui graph save diff_qf, replace
    EraseList `list'

}


    local export_pdf_list densities_raw densities diff_densities qf_raw qf diff_qf

    if "`file'" != "" {
        local file `file'_
    }
    else {
        local file decomp_`=meta_val[8]'_`g1v'_`g2v'_lv`decomp_level'_
    }

    di _n
    foreach i in `export_pdf_list' {
        di "-> `=c(pwd)'/`file'`i'"
        cap erase `file'`i'.pdf
        cap erase `file'`i'.tif

        qui {

            graph use `i'

            if c(os) =="Unix" {
                *graph export `file'`i'.eps, replace
                *sh epstopdf `file'`i'.eps
                *erase `file'`i'.eps
                graph export `file'`i'.tif, replace
            }

            else if c(os) == "Windows" {
                graph export `file'`i'.pdf, replace
            }

            else {
                graph export `file'`i'.tif, replace
            }

            erase "`i'.gph"

        }
    }

    set graphics on
    
end

