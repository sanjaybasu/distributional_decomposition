clear matrix
set mem 4g
set more off
set matsize 11000

do setup_decomp.do

use [INSERT YOUR DATASET NAME HERE.DTA]


/*
Syntax:
decomp (varlist1)(varlist2)...(varlist7), depvar(varname) group(varname) g1(num1) g2(num2) 
       points(integer) 
       [
       weight(svyweight) reps(integer) seed(integer) save(name)
       ]

*/

decomp (smoking)(alcohol)(sodium)(bmi), depvar(bloodpressure) group(race) g1(white) g2(black) ///
    points(#ofpointsalongdistribution) weight(sampleweight) reps(bootstrap) seed(12345)
    *save(name)

/*
Syntax:
GraphResults (string1)(string2)...(string7), depvar(string) group(string) g1(string) g2(string)
    [
    xline(numlist) result_file(file name) save(file name)
    ]

result_file: tells GraphResult where saved results from decomp are stored, otherwise, it looks for
            the file name in rreturn memory
save: save graph name prefix
*/
GraphResults (smoking)(alcohol)(sodium)(bmi), depvar(bloodpressure) group(race) g1(white) g2(black) ///
    result_file(`=r(file)')
    *save(name)















/*
* another example
cd "$dropbox/Working Ado/DFL/example"
use NHANES_example.dta, clear
cd "$dropbox/Working Ado/DFL/example/Testing"
decomp (pad660 pad675 pad680)(dbd895 dbd900)(indhhin2), depvar(bmi) group(ridreth1) g1(3) g2(4) points(50) reps(10) seed(12345)
GraphResults (Activity level)(Diet)(Income), depvar(BMI) group(race) g1(white) g2(black) xline(18.5 25 30) result_file(`=r(file)')
*/
