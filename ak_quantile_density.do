/*
est:
/sum_{i=1}^{m}\phi_{i} \times I( \frac{\hat(q)-x_{i}}{h \times s_{i}}) =  r

where if 
\phi_{i} = 1/m \forall i 

it is unweighted kqfe

*/




*clear mata

mata:

    // bandwidth smoothing parameter                   <------------------------
    function h(M, |w) 
    {
        Vk = 3/(2*(pi()^4))
        sigma = cholesky(mm_var(M,w))
        h = ((Vk*(8/3)*sqrt(pi())*1/rows(M))^(1/5))*sigma
        return(h)
    }

    // kernel function and integrated kernel function   <-----------------------
    function IK(u) return(1:/(1:+exp(-u)))
    function K(u) return(exp(-u):/(1:+exp(-u)):^2)

    // adaptive adjustment term s                       <-----------------------
    function s(x,f,xx) {
        f_xi = spline3eval(spline3(x,f),xx)
        b = exp(sum(log(f_xi))*(1/rows(f_xi)))
        result = (f_xi/b):^(-1/2)
        return(result)
    }

    // constant c
    function cons(h,s,w) return(2*h/K(0)*(1/sum(s:*w)))

    // weighted variance
    function mm_var(X, |w) 
    {
        if (args()==1) w = 1
        CP = quadcross(w,0, X,1)
        n  = cols(CP)
        means = CP[|1\n-1|] :/ CP[n]
        if (missing(CP) | CP[n]==0) return(J(cols(X), cols(X), .))
        return(crossdev(X,0,means, w, X,0,means) :/ CP[n])
    }


    // calculate it all
    void kqfe(string scalar varname, 
                     string scalar touse,
                     | string scalar theta, 
                     string scalar phi,
                     string scalar at
                     )
    {
        real colvector M, survey_weight, phi_weight, percentiles, ww, x, s
        real colvector f, f2, f_k, Q_k
        real scalar n, h, c, p_ll, p_ii, p_ul
        numeric scalar p_rows




        ""
        "-----------------------------------------------------"
        ""
        "kqfe(y [, touse, svy_wgt, phi_wgt, at])"

        if (touse == "" ) {
            "  - No sample set"
            exit()
        }

        // M is data matrix
        st_view(M, ., (varname), touse)
        st_view(survey_weight, ., theta, touse)
        st_view(phi_weight, ., phi, touse)
        n = rows(M)


        // Reshape weights.  If no survey weight specified then 
        // theta_{i} =  1/N otherwise theta_{i} = 1/w(N)
        if (theta == "" ) {
            "  - No survey weights set"
            ww = J(rows(M),1,1/n)
        }
        else {
            ww = survey_weight :/ sum(survey_weight)
        }
        if (phi != "") ww = ww :* phi_weight                  /* <-------------------------- changed this */





        // Percentiles to evaluate
        p_ll = 0.025
        p_ii = 0.025
        p_ul = 0.975
        p_rows = round(((p_ul-p_ll)/p_ii)+1)

        percentiles = J(p_rows,1,p_ll)
        for (i=1; i<=p_rows; i++) {
            percentiles[i,1] = p_ll+(p_ii*(i-1))
        }

        
        // setting grid
        if (at == "") {
            "  - Use default grid"
            x = J(30,1,.)
            for (i=1; i<=30; i++) {
                if (i==1) x[i] = min(M)
                else x[i] = x[i-1]+(max(M)-min(M))/(30-1)
            }
        }
        else {
            "  - Custom grid set"
            x = st_matrix(at)
        }

        ""
        "KERNEL DENSITY ESTIMATION"
        ""
        /*
        Kernel density (f) and adaptive kernel density (f2) with
        parameters bandwidth (h) and adaptive (s)
        */

        "Stage 1: f(.) with fixed h bandwidth parameter"
        h = h(M,ww)
        z = ((M*J(1,rows(x),1)):-x')/h
        Khh = (1/h)*K(z)
        f = (ww'*Khh)'

        "Stage 2: f(.) with adaptive parameter s"
        s = s(x,f,sort(M,1))
        z = ((M*J(1,rows(x),1)):-x'):/(s*h)
        Khh = (1/h)*K(z)
        f2 = (ww'*editmissing(Khh,0))'
        f_k = x, f, f2


        c = cons(h,s,ww)

        ""
        "KERNEL QUANTILE ESTIMATION"
        ""
        "Estimate kernel quantiles with parameters:"
        "  h = "+strofreal(h)
        "  c = "+strofreal(c)

        // initialize q0 as a vector of means
        q0 = J(rows(percentiles),1,1)*mean(M,ww)
        // q_k+1 recurrently.  we can some convergence criterion instead

        max_iter = 1000
        history = percentiles, q0

        for (l=1;l<=max_iter;l++) {
            //"## itter "+strofreal(l)

            //for checking rows,col passed to next itteration
            //rows(q0),rows(q_kp),rows(percentiles)

            z= ((q0*J(1,rows(M),1)):-M'):/(s*h)'
            z = ww'*editmissing(IK(z'),0)
            q_kp = q0 + c*(percentiles - z')
            //" - calculations done"


            /* Convergence when qk-qk+1 = 0.  Where Q(p) becomes zero,
            omit for further calculations */
            check = q0 - q_kp
            check = (abs(check):>=0.00001):*check
            //"check: percentile, check, rounded check,value"
            //percentiles, q0-q_kp, check, q_kp

            /* Reduce percentiles to be calculate Q(p) where Q(p)!=0 */
            q0 = select(q_kp,check:!=0)
            percentiles = select(percentiles,check:!=0)
            check = percentiles,select(check,check:!=0)


            //"update history"
            // Current column from last
            history = history,history[1...,cols(history)]

            /* Update values of q_kp with new values 
               q_k = q_k-1 if check = 0
               q_k = q_k if check != 0 */
            for (ch=1;ch<=rows(history);ch++) {
                for (cp=1;cp<=rows(percentiles);cp++) {
                    // if ((history[ch,1] == check[cp,1]) & (check[cp,2] != 0)) {
                    if (history[ch,1] == check[cp,1]) {
                        history[ch,2+l] = q0[cp,1]
                    }
                }
            }

            q_kp = q0

            if (rows(check) == 0) {
                "Done!!! "+strofreal(l)+" itterations"
                break
            }

        }


        if (l==max_iter) {
            "<<< Warning: no convergence after "+strofreal(l)+" with state: >>>"
            "[p, check,,Q(p)]"
            percentiles,check, q_kp
        }

        ""
        "Data range: "+strofreal(min(M))+", "+strofreal(max(M))
        "Sample: "+strofreal(rows(M))
        "Returns matrix H in Stata for all itterations"

        // update listing for percentiles and grab q0 as last column of history

        Q_k = history[1...,1], history[1...,cols(history)] 


        st_matrix("H", history)

        st_matrix("f", f_k)
        st_matrix("Q", Q_k)

        st_local("h",strofreal(h))
        st_local("c",strofreal(c))

        ""
        "-----------------------------------------------------"
        ""

    }
end

/*

*  make sure to comment out this testing part when using
*




set more off
cd "$dropbox/Working Ado/DFL/.t/"
use "dataset_with_weights.dta", clear


mata: kqfe("mean_sbp", "touseg1", "newweight", "phi_3_3443", "at_range")






** sort out how you're handling  diagnostic outpts
** cleaner way to find mising values

clear

* show convergence of all quantiles
svmat H
rename H1 p
local incr = (colsof(H)-2)/(50-2)
forvalues i =2(`incr')`=colsof(H)' {

    local y `y' H`=int(`i')'
}

figure2, y(`y') x(p) title(q_`=colsof(H)' itterations; (h,c)=(`h', `c')) 


*/

/*
        // qk recurrently
        for (l=1;l<=100;l++) {
            z = (((q0*J(1,rows(M),1)):-M')':/(h*s))
            z = (J(1,rows(M),1)*IK(z))
            q_kp = q0 + c*(n*r - z')
            q0 = q_kp
        }
*/


*/
