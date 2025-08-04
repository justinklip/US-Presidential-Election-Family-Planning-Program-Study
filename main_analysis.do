/// Yet-To-Be-Treated Regression Specification

/// Setup
clear
eststo clear
use "/Users/klipjustin/Desktop/ECO499 Paper/Raw Data/analysis.dta", clear

rename v1 state
rename v2 county
replace dem = subinstr(dem, "%", "", .)
replace rep = subinstr(rep, "%", "", .)
replace other = subinstr(other, "%", "", .)
destring dem, replace
destring rep, replace

* Generate State Numbers
egen stateid = group(state)

* Fix the year measure
replace fp_year_p74_fed = fp_year_p74_fed + 1900

gen yeardist_fp = year - fp_year_p74_fed

/// Treatment Setup
* Narrow it down to just counties that do receive family planning programs

gen treatment_group = 0
replace treatment_group = 1 if !missing(fp_year_p74_fed)

gen postperiod = 0
replace postperiod = 1 if year >= fp_year_p74_fed

gen postperiodxtreatment = 0
replace postperiodxtreatment = 1 if postperiod == 1 & treatment_group == 1

/// Add Some Controls

* Define control variables
gen pop100 = pop*100000
global controls "land_area pop100 popdensity urbshare nonwhiteshare pop65 medage popforeign medfaminc mededuc unemployment_rate"

foreach y in 1948 1952 1956 1960 1964 1968 1972 {
    gen year`y' = (year == `y')
}
* Generate interactions of controls with year dummies
foreach v of global controls {
    foreach y in 1948 1952 1956 1960 1964 1968 1972 {
        gen `v'_yr`y' = `v' * year`y'
    }
}
gen dem_rep_diff = dem - rep
* Generate Previous Democratic Mean Variables Based on the average of 1952, 1956, and 1960
egen past_dem_mean = mean(cond(inlist(year, 1952, 1956, 1960), dem, .)), by(fips)
bysort fips (year): replace past_dem_mean = past_dem_mean[1]

*Generate the win margin variables
egen past_dem_rep_diff = mean(cond(inlist(year, 1948, 1952, 1956, 1960, 1964), dem_rep_diff, .)), by(fips)
bysort fips (year): replace past_dem_rep_diff = past_dem_rep_diff[1]

/// Regression Specifications - Compare treated units to not-yet-treated units

* Without controls, clustering by county
reg dem postperiod treatment_group postperiodxtreatment i.stateid ///
    if !missing(fp_year_p74_fed) & inlist(year, 1956, 1960, 1964, 1968), cluster(fips)
eststo

* With controls
reg dem postperiod treatment_group postperiodxtreatment i.stateid ///
    $controls if !missing(fp_year_p74_fed) & inlist(year, 1964, 1968), cluster(fips)
eststo

esttab using not_yet_treated_regression.tex, replace se starlevels(* 0.10 ** 0.05  *** 0.01)
eststo clear

/// Balance Test
reghdfe fp_year_p74_fed land_area pop100 popdensity urbshare nonwhiteshare pop65 medage popforeign ///
    medfaminc mededuc unemployment_rate past_dem_mean if year == 1964, cluster(fips)
eststo
esttab using balance_test_not_yet_treated.tex, replace se starlevels(* 0.10 ** 0.05 *** 0.01)
eststo clear



/// Full Regression w/ as many pre-periods as possible
* County and Year FE's
reghdfe dem postperiodxtreatment i.year if !missing(fp_year_p74_fed) & year <= 1972, ///
    cluster(fips) absorb(i.fips)
eststo

* Adding controls
reghdfe dem postperiodxtreatment i.year $controls if !missing(fp_year_p74_fed) & year <= 1972, ///
    cluster(fips) absorb(i.fips)
eststo


* Adding year differential effects of controls
local interactions
foreach y in 1948 1952 1956 1960 1964 1968 1972 {
    foreach v of global controls {
        local interactions `interactions' `v'_yr`y'
    }
}
gen postxtreatxpast_diff = postperiodxtreatment * past_dem_rep_diff 

*Interactions 
reghdfe dem postperiodxtreatment i.year `interactions' if !missing(fp_year_p74_fed) & year <= 1972, ///
    cluster(fips) absorb(i.fips)
	
eststo
esttab using regs_not_yet_treated.tex, replace se starlevels (* 0.10 ** 0.05 *** 0.01) keep(postperiodxtreatment *year )
eststo clear
	
*Heterogeneity Analysis
reghdfe dem postperiodxtreatment i.year `interactions' postxtreatxpast_diff if !missing(fp_year_p74_fed) & year <= 1972, ///
    cluster(fips) absorb(i.fips)

eststo
esttab using heterogeneity_not_yet_treated.tex, replace se starlevels (* 0.10 ** 0.05 *** 0.01) keep(postperiodxtreatment postxtreatxpast_diff *year )
eststo clear

*Binning
gen postxtreat = postperiodxtreatment
	
* Step 1: Create deciles of past_dem_rep_diff
xtile decile_past_diff = past_dem_rep_diff, n(10)


* Step 2: Run regression with interaction between decile and treatment
reghdfe dem i.decile_past_diff##i.postxtreat i.year i.year#i.decile_past_diff `interactions' ///
    if !missing(fp_year_p74_fed) & year <= 1972, ///
    cluster(fips) absorb(i.fips)
* Step 3: Plot interaction effects by decile
lincom 1.postxtreat + 1.decile_past_diff#1.postxtreat
lincom 1.postxtreat + 2.decile_past_diff#1.postxtreat
lincom 1.postxtreat + 3.decile_past_diff#1.postxtreat
lincom 1.postxtreat + 4.decile_past_diff#1.postxtreat
lincom 1.postxtreat + 5.decile_past_diff#1.postxtreat
lincom 1.postxtreat + 6.decile_past_diff#1.postxtreat
lincom 1.postxtreat + 7.decile_past_diff#1.postxtreat
lincom 1.postxtreat + 8.decile_past_diff#1.postxtreat
lincom 1.postxtreat + 9.decile_past_diff#1.postxtreat
lincom 1.postxtreat + 10.decile_past_diff#1.postxtreat

gen decile = .
gen effect = .
gen lower = .
gen upper = .

forvalues i = 1/10 {
    lincom 1.postxtreat + `i'.decile_past_diff#1.postxtreat
    replace decile = `i' in `i'
    replace effect = r(estimate) in `i'
    replace lower = r(lb) in `i'
    replace upper = r(ub) in `i'
}

twoway (scatter effect decile, msymbol(circle) msize(medium)) ///
       (rcap lower upper decile), ///
       ytitle("Treatment Effect (β₂ᵢ + β₃)") ///
       xtitle("Decile of Margin") ///
       title("Estimated Treatment Effects by Decile") ///
       xlabel(1(1)10, valuelabel angle(0)) ///
       legend(off) ///
       graphregion(color(white)) ///
       scheme(s2mono)
	   
* Step 3: Plot interaction effects by decile
coefplot, keep(*.postxtreat) xline(0) ///
    rename( ///
        1.decile_past_diff#i.postxtreat = "D1" ///
        2.decile_past_diff#i.postxtreat = "D2" ///
        3.decile_past_diff#i.postxtreat = "D3" ///
        4.decile_past_diff#i.postxtreat = "D4" ///
        5.decile_past_diff#i.postxtreat = "D5" ///
        6.decile_past_diff#i.postxtreat = "D6" ///
        7.decile_past_diff#i.postxtreat = "D7" ///
        8.decile_past_diff#i.postxtreat = "D8" ///
        9.decile_past_diff#i.postxtreat = "D9" ///
        10.decile_past_diff#i.postxtreat = "D10" ///
    ) ///
    title("Treatment Effects by Decile of Past Dem-Rep Diff") ///
    ylabel(, angle(horizontal)) xtitle("Estimated Effect on dem") ///
    xlabel(, labsize(small)) msymbol(O) mcolor(black) ciopts(lcolor(black))
