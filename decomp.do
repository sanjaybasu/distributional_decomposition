* add [if] [in]



cap program drop DecompEstimation
program define DecompEstimation, eclass

    syntax anything, depvar(varlist max=1) ///
           group(varlist max=1) g1(numlist max=1 integer) g2(numlist max=1 integer) ///
           [, ///
           weight(varlist max=1) dec ///
           ] 

    * `anything':
    * - extract patterns of (a b c) and store as (`x1')
    * - total collection of variables for sample `ungrouped_varlist'
    * -find decomposition order `decomp_level'
    local varlist `anything'
    local grouped_varlist `anything'
    local ungrouped_varlist
    local i = 1
    while regexm("`anything'", "(\()([A-Za-z0-9_ ]*)(\))(.*)")==1 {
        *di regexs(2)    
        local x`i' = regexs(2)      
        local anything = regexs(4)
        local ungrouped_varlist `ungrouped_varlist' `=regexs(2)'
        local i = `i' + 1
    }
    confirm variable `ungrouped_varlist'

    * Display some output to describe decomposition
    local decomp_level = `i' - 1
    di _n "Model: f(`depvar'|z,G=g) = `: di _dup(`decomp_level') "\int"'"_c
    di "f(`depvar',x1,...,x`decomp_level'|`group'=g) dx1,...,dx`decomp_level'"
    di "where z = `grouped_varlist'"    
    forvalues i = 1(1)`decomp_level' {
        di "     x`i' = `x`i''"
    }
    di "     G = `group'" _n
    
    * Dummy variables and inclusion/ exclusion
    tempvar touse touseg1 touseg2 group_g1 group_g2 w
    mark `touse' if ((`group' == `g1') | (`group' == `g2')) & (`depvar'>at_range[1,1] & `depvar'<at_range[rowsof(at_range),1])
    markout `touse' `ungrouped_varlist' `depvar' `weight'
    mark `touseg1' if `touse' & `group' == `g1'
    mark `touseg2' if `touse' & `group' == `g2'
    qui gen `group_g1'=`group'==`g1'
    qui gen `group_g2'=`group'==`g2'

    * Survey weights
    if "`weight'" != "" {
        di "Survey weights (`weight') specified.  Weight phiw:=`weight'*phi "_c
        di "where `weight' normalized to sum to 1."_n
        qui sum `weight' if `touse'
        gen `w' = `weight'/r(sum) if `touse'
    }
    else {
        qui gen `w' = 1/_N if `touse'
    }

    ****************************************************************************
    * Calculate all DFL phi weights                                            *
    ****************************************************************************
    di "Reweighting factors generated from:"_n

    forvalues i = 1(1)`decomp_level' {
        local varlist1
        local varlist2

        forvalues j = `i'(1)`decomp_level' {
            local varlist1 `varlist1' `x`j''
        }
        local i_plus = `i'+1
        forvalues j = `i_plus'(1)`decomp_level' {
            local varlist2 `varlist2' `x`j''
        }
        *di "Varlist 1 `varlist1'"
        *di "Varlist 2 `varlist2'"

        * Calculating weights: phi = (s/t)*(u/v) = (1-t)/t * u/(1-u)
        tempvar t u
        qui{
            logit `group_g1' `varlist1' if `touse'
            predict double `t' if `touse', p
            logit `group_g1' `varlist2' if `touse'
            predict double `u' if `touse', p
        }
        * Create subscripts `g1'_(`g2')`g1'(`g2')
        local subscript `g1'_`: di _dup(`=`i'-1') "`g1'"'`g2'`: di _dup(`=`decomp_level'-`i'') "`g1'"'
        cap drop phi_`subscript'
        qui gen double phi_`subscript' = ((1-`t')/`t')*(`u'/(1-`u'))

        di "phi_`subscript' = P(G=`g2'|`varlist1') * P(G=`g1'|`varlist2')"  /*Numerator*/
        di "`: di _dup(`=length("phi_`subscript'")') " "'   "_c             /*LHS spacing*/
        di "`: di _dup(`=length("P(G=`g2'|`varlist1')")') "-"'   "_c        /*fraction*/
        di "`: di _dup(`=length("P(G=`g1'|`varlist2')")') "-"'   "          /*fraction*/
        di "`: di _dup(`=length("phi_`subscript'")') " "'   "_c             /*LHS spacing*/
        di "P(G=`g1'|`varlist1') * P(G=`g2'|`varlist2')"_n                  /*Denominator*/
    }

    * Generates the list of all possible (derived; composite) reweighting factors and the 
    * equation to find them.  Ex: We have three variables from i,j,k and groups 0,1. 
    * 1. There are 2^n ways we can select group membership for each indep. variable 
    *    ie: x1;g=1,x2;g=0,x3;g=0 and (2^n)-1 without the trivial case where all from group 0
    * 2.  The set of numbers {1...7} and {base2(1...7)}.  In base2, each element represented
    *     by 3 places (011,111) which we can manipulate to form syntax of reweights.    
    * 3.  Check against list of unique weights (found from estimaton)
    *     Ex: In a 3 way decomposition, {1,2,4} = {001,010,100} are thrown
    local throw_count = 1
    forvalues i = 1(1)`=(2^`decomp_level')-1'{
        qui inbase 2 `i'        
        * Add preceding zeros to base2 value according to maximum digit places among range of elements in base2.
        * for base2(7) = 111, we have 3 places.  In base2(1) = 1 -> 001 
        local base2 `: di _dup(`=(`decomp_level' - length("`=r(base)'"))') "0"'`=r(base)'           
        local weight phi_`g1'_`=subinstr(subinstr(subinstr("`base2'","0",")",.),"1","`g2'",.),")","`g1'",.)'  
        * for i=1(1)7, i = 1,2,4 are the unique weights "
        local reweight_list `reweight_list' `weight'        
        if `i' == 2^(`throw_count'-1) {
            local throw_count = `throw_count' + 1           
        }
        else {
            di _n "(`i') `weight' =" _c
            * drop weight
            cap drop `weight'
            qui gen double `weight' = 1
            local multip = 0
            forvalues j = 1(1)`=length("`base2'")' {    
                * accross sequence 1212, we have '2' at position 3,4 -> 1211 and 1112 "
                if substr(subinstr("`weight'","phi_`g1'_","",1),`j',1) == "`g2'" {
                    local star *
                    if `multip' == 0 local star 
                    local multip = `multip' + 1
                    di "`star' phi_`g1'_`: di _dup(`=`j'-1') "`g1'"'`g2'`: di _dup(`=length("`base2'")-`j'') "`g1'"' "_c
                    qui replace `weight' = `weight'*phi_`g1'_`: di _dup(`=`j'-1') "`g1'"'`g2'`: di _dup(`=length("`base2'")-`j'') "`g1'"'
                    *"
                }
            }
        }
    }

    global reweight_list `reweight_list'

    ****************************************************************************
    * Estimate functions + decompose differences                               *
    ****************************************************************************

    *tab `touseg1' `touseg2'
    di _n

    * For g1
    mata: kqfe("`depvar'", "`touseg1'", "`w'", "", "at_range")
    di _n "> F`g1'`: di _dup(`decomp_level') "`g1'"'(`h',`c') "
    matrix Q`g1'`: di _dup(`decomp_level') "`g1'"' = Q[1...,2]
    matrix f`g1'`: di _dup(`decomp_level') "`g1'"' = f[1...,3]

    * For g2
    mata: kqfe("`depvar'", "`touseg2'", "`w'", "", "at_range")
    di _n "> F`g2'`: di _dup(`decomp_level') "`g2'"'(`h',`c') "
    matrix Q`g2'`: di _dup(`decomp_level') "`g2'"' = Q[1...,2]
    matrix f`g2'`: di _dup(`decomp_level') "`g2'"' = f[1...,3]


    * For counterfactual outcomes
    if "`dec'" == "dec" {
        foreach phi in `reweight_list' {
            local subscript = subinstr(subinstr("`phi'","_","",.),"phi","",.)
            mata: kqfe("`depvar'", "`touseg1'", "`w'", "`phi'", "at_range")
            di _n "> F`g1'`subscript'(`h',`c') "
            matrix Q`subscript' = Q[1...,2]
            matrix f`subscript' = f[1...,3]
        }
        qui AverageDecompOrder, decomp_level(`decomp_level') g1(`g1') g2(`g2') prefix(f Q)

        matrix b = r(result)
        matrix rownames b = b
        ereturn post b, esample(`touse')
    }


end




cap program drop decomp
program define decomp, rclass

syntax anything, depvar(varlist max=1) ///
       group(varlist max=1) g1(numlist max=1 integer) g2(numlist max=1 integer) ///
       [ points(numlist max=1 integer) ///
       weight(varlist max=1) ///
       reps(numlist max=1 integer) seed(numlist max=1 integer) save(string asis) ///
       ]

    di _n "*****************************"
    di    "* decomp (updated 2014/5/8) *"
    di    "*****************************" _n

    * Mark sample for finding range to evaluate kernel density
    local list `anything'
    while regexm("`list'", "(\()([A-Za-z0-9_ ]*)(\))(.*)")==1 {
        local ungrouped_varlist `ungrouped_varlist' `=regexs(2)'
        local list = regexs(4)
        local decomp_level = `decomp_level' + 1
    }

    tempvar touse
    mark `touse' if `group' == `g1' | `group' == `g2'
    markout `touse' `ungrouped_varlist' `depvar' `weight'



    * Grid
    qui count if `touse'== 1 & `group' == `g1'
    local n1 = r(N)
    qui count if `touse' == 1 & `group' == `g2'
    local n2 = r(N)

    * Restrict upon where we have support
    /*
    qui sum `depvar' if `touse' & `group' == `g1'
    local min = r(min)
    local max = r(max)
    qui sum `depvar' if `touse' & `group' == `g2'
    if r(min) > `min' local min = r(min)
    if r(max) < `max' local max = r(max)
    */

    * Should we restrict on the smallest sample?
    qui sum `depvar' if `touse'
    local min = r(min)
    local max = r(max)

    matrix at_range = J(`points',1,`min')
    foreach i of numlist 2/`points' {
        matrix at_range[`i',1] = at_range[`=`i'-1',1] + ((`max'-`min')/`points')
        local at_range_previous = at_range[`i',1]
    }    
    matrix colnames at_range = x

    * Bootstrap settings
    if "`reps'" != "" {
        local bs_setting bs _b, seed(`seed') reps(`reps') notable noheader nolegend :
        local se_
    }
    else {
        local bs_setting
        local se_
    }

    preserve
    keep `touse' `ungrouped_varlist' `depvar' `group' `weight'


    ****************************************************************************
    * Run decompositions                                                       *
    ****************************************************************************

    tempname result f Q save_mat fg1 fg2 Qg1 Qg2

    *** 1) Decomp
    `bs_setting' DecompEstimation `anything' , depvar(`depvar') group(`group') g1(`g1') g2(`g2') weight(`weight') dec
    matrix `result' = e(b)
    cap matrix `result' = e(b)\e(se)
    matrix `f' = `result'[1..rowsof(`result'),1..(`decomp_level'+2)*`points']'
    matrix `Q' = `result'[1..rowsof(`result'),((`decomp_level'+2)*`points')+1..colsof(`result')]'

    *** 2) Observed
    DecompEstimation `anything' , depvar(`depvar') group(`group') g1(`g1') g2(`g2') weight(`weight')

    foreach F in f Q {
        matrix colnames `F'`g1'`: di _dup(`decomp_level') "`g1'"' = "`F'`g1'`: di _dup(`decomp_level') "`g1'"'"
        matrix colnames `F'`g2'`: di _dup(`decomp_level') "`g2'"' = "`F'`g2'`: di _dup(`decomp_level') "`g2'"'"
        matrix colnames `F'`g1'`: di _dup(`decomp_level') "`g2'"' = "`F'`g1'`: di _dup(`decomp_level') "`g2'"'"
    }
   

    *** 3) Saving results

    clear

    qui {

    set obs `points'

    * Metadata for run settings
    gen str meta_parm = ""
    gen str meta_val = ""
    format %9s meta_val

    #delimit ;
    replace meta_parm = "decomp level" in 1      ; replace meta_val = "`decomp_level'" in 1 ;
    replace meta_parm = "bootstrap prefix" in 2  ; replace meta_val = "`bs_setting'" in 2 ;
    replace meta_parm = "functions" in 3         ; replace meta_val = "f Q" in 3 ;
    replace meta_parm = "group 1" in 4           ; replace meta_val = "`g1'" in 4 ;
    replace meta_parm = "group 2" in 5           ; replace meta_val = "`g2'" in 5 ;
    replace meta_parm = "obs 1" in 6             ; replace meta_val = "`n1'" in 6 ;
    replace meta_parm = "obs 2" in 7             ; replace meta_val = "`n2'" in 7 ;
    replace meta_parm = "depvar" in 8            ; replace meta_val = "`depvar'" in 8;
    #delimit cr

    * Observed functions
    foreach F in f Q {
        svmat `F'`g1'`: di _dup(`decomp_level') "`g1'"', names(col)
        svmat `F'`g2'`: di _dup(`decomp_level') "`g2'"', names(col)
        svmat `F'`g1'`: di _dup(`decomp_level') "`g2'"', names(col)
    }

    * Decomposition results
    foreach F in f Q {
        local n = rowsof(``F'')/(`decomp_level'+2)
        forvalues j = 1(1)`=`decomp_level'+2' {
            matrix `save_mat' = ``F''[((`j'-1)*`n')+1..(`j'*`n'),1..colsof(``F'')]
            matrix colname `save_mat' = b
            cap matrix colnames `save_mat' = b se
            local name : word 1 of `: roweq `save_mat' '
            matrix coleq `save_mat' = `name'_:
            svmat `save_mat', names(eqcol)

        }
    }

    * Estimated counterfactual outcome (f_g2 + effect)
    foreach F in f Q {
        forvalues j = 1(1)`decomp_level' {
            local E = char(`j'+64)
            gen `F'`g2'`: di _dup(`decomp_level') "`g2'"'_effect_`E' = ///
            `F'`g2'`: di _dup(`decomp_level') "`g2'"' + `F'_effect_`E'_b
        }
    }

    * X-axis (x,p)
    svmat at_range, names(col)
    order x, before(f`g1'`: di _dup(`decomp_level') "`g1'"')
    matrix p = Q[1...,1]
    matrix colname p = "p"
    svmat p, names(col)
    order p, before(Q_difference_b)
    order Q`g1'`: di _dup(`decomp_level') "`g1'"' Q`g2'`: di _dup(`decomp_level') "`g2'"', after(p)

    } /*qui*/




    if "`save'" == "" {
        local save decomp_`depvar'_`g1'_`g2'_lv`decomp_level'
    }

    qui save `save', replace
    di _n "Results saved as `save'"
    return local file = "`save'"


    restore

end
