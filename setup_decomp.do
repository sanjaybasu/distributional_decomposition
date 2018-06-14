clear all 
cap pr drop _all
set more off

do "decomp.do"
do "Dep/decomp_order.do"
do "Dep/ak_quantile_density.do"
do "Dep/graph_results.do"
do "Dep/figure2.do"


