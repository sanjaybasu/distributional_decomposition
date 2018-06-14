
* Given 
*


* The group conditional characteristic density f(z|G=g) can be itterated as
* f(z1|z2...zn,G=g) * f(z2|z3...zn,G=g) * ... (zn|G=g).  In our notation, the
* density of outcome y we can express as fg*ggg..., so if g = {1,2}, we 
* can say f1121

* In a sequential decomposition involving n characteristics, we have n! 
* unique decomposition orders.  Suppose n=3, one possible order is:
* f1-f2 = (f1111-f1112) + (f1112-f1122) + (f1122-f1222) + (f1222-f2222)     (1)

* A particular decomposition is characterized by the order in which we identify
* the effects.  In f1111-f1112 identifies the difference in density of y due
* to differences in conditional densities of z3.  See f1(11)1-f(11)2 because we're
* reweighing from group 1.  Thus, in (1) we're identifying the effects in the
* order of (z3,z2,z1) or (3,2,1) for short.  

* A) So first, we need to list all orders:

    * Johnson-Trotter Algorithm for listing all permutations 
    * http://www.cut-the-knot.org/Curriculum/Combinatorics/JohnsonTrotter.shtml

    * An directed integer <k(>) is called leftward (rightward) directed 
    * if it traverses to the left (right) in the set of permutations

    * A directed integer is said to be mobile if it is greater than its immediate 
    * neighbor in the direction it is looking at.

    * Initialize the first permutation with <1 <2 <3 ... <n 
    * (prog denotes < as -1, > as 1)
    * while there exists a largest mobile integer MMI                        (A)
    *    swap MMI and the adjacent integer it is looking at and              (B)
    *    reverse the direction of all integers larger than MMI

* B) Having the list, then we can play with substrings and positions of substrings
*    to create the sequence


cap prog drop MaxMobileInteger
prog define MaxMobileInteger, rclass 
syntax, matrix(name)

    * Given a sequence of mobile integers, find largest mobile integer
    tempname order M MM
    matrix `order' = `matrix'
    matrix `M' = J(rowsof(`order'),colsof(`order'),0)
    * Find list of mobile integers
    local mobile_integers 
    forvalues i = 1(1)`=colsof(`order')' {
        if ((`order'[1,`i'] > `order'[1,`=`i'-1']) & (`order'[2,`=`i''] == -1) & `i'-1 > 0 ) ///
        |  ((`order'[1,`i'] > `order'[1,`=`i'+1']) & (`order'[2,`=`i''] == 1) & `i'+1 <= `=colsof(`order')' ) {
            matrix `M'[1,`i'] = `order'[1,`i']
            matrix `M'[2,`i'] = `order'[2,`=`i'']
        }
        else {
            * Is not a mobile integer
        }
    }
    
    * Sort to find largest MMI 
    mata: st_matrix("`MM'", sort(st_matrix("`M'")',1)[rows(st_matrix("`M'")'),])

    if (`MM'[1,1] != 0) & (`MM'[1,2] != 0) {
        return local end = 0
        return local mmi_dir = `=`MM'[1,2]'
        return local mmi_val = `=`MM'[1,1]'
    }
    else {
        return local end = 1
    }

end

* SwapMMIwDNeighbour(current_order, position)
* Swaps entry for MMI in current_order with its nearest directional neighbour
capture mata mata drop SwapMMIwDNeighbour()
mata:
    void SwapMMIwDNeighbour(string current_order, scalar position)
    {
        c = st_matrix(current_order)
        // direction = {-1,1} = {left,right}
        MMI = c[.,position]
        direction = MMI[2,1]
        // Swap
        c[.,position] = c[.,position+direction]
        c[.,position+direction] = MMI[.,.]      
        st_matrix(current_order, c)     
    }
end


cap prog drop GetNewSequence
prog define GetNewSequence, rclass 
syntax, matrix(name) mmi_val(numlist max=1 integer)
    
    * Given a list of mobile integers and the value of the MMI, swap positions with
    * neighbour it is facing, then, for all other values greater than MMI, reverse
    * direction
    tempname orderI orderG
    matrix `orderI' = `matrix'
    matrix `orderG' = `matrix'

    *** Swap largest mobile integer with nearest neighbour in which it is facing
    forvalues i = 1(1)`=colsof(`orderI')' {
        * Go through all entries of matrix and check if it is MMI value.  If it is,
        * then [MMI value, direction]' = order_i[1,`i']
        * Index position identified at `i' and swap
        if (`orderI'[1,`i'] == 0) {
            exit
        }
        if ((`orderI'[1,`i']) == `mmi_val') {
            local mmi_dir = `orderI'[2,`i']
            mata: SwapMMIwDNeighbour("`orderG'",`i')
        }
    }

    *** Reverse directionality of all integers greater than largest mobile integer
    forvalues i = 1(1)`=colsof(`orderG')' {
        if (`orderG'[1,`i'] > `mmi_val') {
            matrix `orderG'[2,`i'] = (`orderG'[2,`i'])*(-1)
        }
    }

    return matrix order_new = `orderG'
end




cap prog drop AverageDecompOrder
prog define AverageDecompOrder, rclass 
syntax, decomp_level(integer) g1(integer) g2(integer) [prefix(namelist)]

    * Finds the average contribution due to each of the characteristics
    * in z using the matrices fijkl which are stored in memory.
    * (This will run even if matrices not defined, only, it will just
    * output the sequences without doing anything).

    * If calculating, can take list of matrix prefixes, ex: prefix(f Q)

    di _n
    di _n "`decomp_level' way decomposition with `=floor(exp(lnfactorial(`decomp_level')))' decomposition orders" _n


    if `decomp_level' > 7 {
        di as red "Maximum permitted is a 7 way decomposition due to matrix storage"
        di as red "limit of 11,000"
        exit
    }


    foreach pr in `prefix' {


        local f `pr'

        local prefix_c = `prefix_c'+1


        * Create matrix representing list of mobile integers
        tempname order_i Nlist
        mata : st_matrix("`order_i'", 1..strtoreal(st_local("decomp_level")) \ J(1,strtoreal(st_local("decomp_level")) ,-1))
        * Take integers without direction to compile list of permutations
        matrix `Nlist' = `order_i'[1,1...]

        * Create storage for calculated effects
        forvalues i = 1(1)`=colsof(`order_i')' {
            local compare1 : di _dup(`=colsof(`order_i')') "`g1'"

            matrix effect_`=char(`=`i'+64')' = J(rowsof(`f'`g1'`: di _dup(`=colsof(`order_i')') "`g1'"'),1,0)
            matrix colnames effect_`=char(`=`i'+64')' = "effect `=char(`=`i'+64')'"
        }

        *** Johnson-Trotter to find permutations
        MaxMobileInteger, matrix(`order_i')
        while `=r(end)' == 0 {                                         /* J-T step (A) */
            GetNewSequence, matrix(`order_i') mmi_val(`=r(mmi_val)')   /* J-T step (B) */
            matrix `order_i' = r(order_new)                            /* Resulting matrix */
            matrix `Nlist' = (`Nlist')\(`order_i'[1,1...])
            MaxMobileInteger, matrix(`order_i')    
        }

        * Given a decompositon order, we want to create the decompositon expression 
        * f1ijk = f1xyz = compare1 - compare2.  This is just playing with strings

        forvalues i = 1(1)`=rowsof(`Nlist')' { 
            matrix `order_i' = `Nlist'[`i',1...]
            matrix list `order_i', noheader nonames
            forvalues j = 1(1)`=colsof(`order_i')' {
                local effect = `order_i'[1,`j']
                local E = char(`=`order_i'[1,`j']+64')
                * Set first item of decomposition sequence
                if `j' == 1 {
                    local compare1 : di _dup(`=colsof(`order_i')') "`g1'"
                }
                * Set second item of decompositon sequence based on first item
                local compare2 "`=substr("`compare1'",1,`=`effect'-1')'`g2'`=substr("`compare1'",`=`effect'+1',.)'"
                * replace 
                * The expression to calculate effect is given by
                local quantity `f'`g1'`compare1' - `f'`g1'`compare2'

                * - Add to matrix of average effect
                cap matrix effect_`E' = effect_`E' + 1/`=(exp(lnfactorial(colsof(`order_i'))))'*(`quantity')
                cap matrix rownames effect_`E' = `f'_effect_`E':
                di "effect_`E' + (`quantity')" 

                * Update first item' as current second item and proced to next component
                * of decomposition expression
                local compare1 `compare2'
            }

        }


        * Compute raw difference and residual
        * Note, we're labelling the matrices because we're exploiting a trick for bootstraping
        * where we can bootstrap a set of results saved as the b matrix.  This relies on matrix 
        * column names (matrix') having unique names

        cap {
            matrix diff = `f'`g1'`: di _dup(`=colsof(`order_i')') "`g1'"' - `f'`g2'`: di _dup(`=colsof(`order_i')') "`g2'"'
            matrix rownames diff = `f'_difference:
            matrix resid = `f'`g1'`: di _dup(`=colsof(`order_i')') "`g2'"' - `f'`g2'`: di _dup(`=colsof(`order_i')') "`g2'"'
            matrix rownames resid = `f'_residual:

            if `prefix_c' == 1 {
                forvalues i = 1(1)`=colsof(`order_i')' {
                    local E = char(`i'+64)
                    if `i' == 1 {
                        matrix result = (effect_`E')'
                    }
                    else {
                        matrix result = result,(effect_`E')'
                    }
                }
                matrix result = diff', result, resid'
            }

            else {                
                matrix result = result,diff'
                forvalues i = 1(1)`=colsof(`order_i')' {
                    local E = char(`i'+64)
                    matrix result = result,(effect_`E')'
                }
                matrix result = result,resid'    
            }



 

        }






    }






    return matrix result = result



end


/*
set matsize 11000
* for testing... has fijkls 
*use testing.dta, clear
*mkmat f*, nomissing
AverageDecompOrder, decomp_level(4) g1(1) g2(2) f
*/



