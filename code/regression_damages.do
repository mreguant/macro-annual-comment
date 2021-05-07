* modify model to have kink at 1.2

global dirpath = "/Users/mreguant/Documents/git/macro-annual-comment"

clear all
set obs 100
gen temp = 1.2 + (_n-1) * 1.0/15.0
gen temp2 = temp^2/2.0
gen temp3 = (temp-2.0)^2 * (temp >= 2.0)
gen temp3_alt = (temp-1.2)^2

global g1 = 0.000177
global g2 = 0.0044
global g31 = 0.0
global g32 = 0.039
global g33 = 0.772
global thresh = 3.0

forvalues i=1(1)3 {
	gen d`i' =  temp * $g1 + temp^2 * $g2 / 2.0 + ${g3`i'} * temp3
	gen N`i' = exp(d`i')
	gen invN`i' = 1/N`i'
}

* partialled out
gen d2_po = d2 - $g1 * temp - $g2 * temp2
gen d3_po = d3 - $g1 * temp - $g2 * temp2


* approximate function
forvalues d=2(1)3 {
	reg d`d'_po temp3_alt if temp < $thresh, nocons
	gen d`d'_fit = $g1 * temp + $g2 * temp2 + _b[temp3_alt] * temp3_alt
	gen N`d'_fit = exp(d`d'_fit)
	gen invN`d'_fit = 1/N`d'_fit
}


twoway (scatter d3 temp) (line d3_fit temp) if temp < $thresh + 1.0
twoway (scatter N3 temp) (line N3_fit temp) if temp < $thresh + 1.0
twoway (scatter invN3 temp) (line invN3_fit temp) if temp < $thresh + 1.0
 

twoway line invN1 temp || line invN2 temp || line invN2_fit temp || line invN3 temp || line invN3_fit temp, ///
	ytitle("1/N") xtitle("Temp. increase in C") ///
	legend(order(1 "Case 1" 2 "Case 2" 4 "Case 3" 3 "Case Approx. 2" 5 "Case Approx. 3"))
graph export "/Users/mreguant/Documents/git/macro-annual-comment/output/damages.pdf" 
