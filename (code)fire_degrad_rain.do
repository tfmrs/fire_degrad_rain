***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH
***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH
***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH***1st SUBPATH

*cd "G:\water_pp\Comet\"
*cd "F:\Dell_HD\water_pp"
cd "D:\Pesquisa\Pesquisa_2022\water_pp"
clear all
use "(dataset)fire_degrad_rain_1st_2nd_subpaths.dta", clear


global xlist_1st focos_up focos_nup ppt evapo area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area d2m t2m w s pib pop ts_dist_sede_100milhab ts_dist_roads d_ano2-d_ano12 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9

global xlist_1st_naive focos_agr ppt evapo area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area d2m t2m w s pib pop ts_dist_sede_100milhab ts_dist_roads d_ano2-d_ano12 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9

global xlist_1st_degrad focos_up focos_nup ppt evapo area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area degrad_area d2m t2m w s pib pop ts_dist_sede_100milhab ts_dist_roads d_ano2-d_ano9 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9

global xlist_1st_naive_degrad focos_agr ppt evapo area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area degrad_area d2m t2m w s pib pop ts_dist_sede_100milhab ts_dist_roads d_ano2-d_ano12 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9

*03/11 17h11: no SE in 2nd
global xlist_2nd focos_for area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area cloud d2m t2m w s d_ano2-d_ano12 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9

global xlist_2nd_degrad focos_for area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area degrad_area cloud d2m t2m w s d_ano2-d_ano9 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9

summ ppt evapo focos_for $xlist_1st if d_miss == 0
summ aod if d_miss_aod == 0
summ slope_pct if d_miss == 0
summ degrad if d_miss == 0 & ano <= 2016



*****(A) NAIVE FE ESTIMATION
**First stage
*xtreg focos_for $xlist_1st if d_miss == 0, fe vce(robust)
*estimates store d1st
*test

******WITHOUT DEGRAD******WITHOUT DEGRAD******WITHOUT DEGRAD******WITHOUT DEGRAD


**(A) First stage: 
***(A.1) naive
xtreg focos_for $xlist_1st_naive if d_miss == 0, fe vce(robust)
estimates store d1st_naive
test

***(A.2) R&V
xtreg focos_for $xlist_1st if d_miss == 0, fe vce(robust)
estimates store d1st_rv
test


**(B) Second stage
**(B.1) naive: ppt & evapo
foreach var of varlist ppt evapo {
xtreg `var' $xlist_2nd if d_miss == 0, fe vce(robust)
estimates store d2nd_`var'
test
}

**(B.2) naive: AOD
xtreg aod $xlist_2nd if d_miss_aod == 0, fe vce(robust)
estimates store d2nd_aod
test


******WITH DEGRAD******WITH DEGRAD******WITH DEGRAD******WITH DEGRAD******WITH DEGRAD******WITH DEGRAD******WITH DEGRAD******WITH DEGRAD

**(A) First stage: 
***(A.1) naive
xtreg focos_for $xlist_1st_naive_degrad if d_miss == 0 & ano <=2016, fe vce(robust)
estimates store d1st_naive_degrad
test

***(A.2) R&V
xtreg focos_for $xlist_1st_degrad if d_miss == 0 & ano <=2016, fe vce(robust)
estimates store d1st_rv_degrad
test


**(B) Second stage
**(B.1) naive: ppt & evapo
foreach var of varlist ppt evapo {
xtreg `var' $xlist_2nd_degrad if d_miss == 0 & ano <=2016, fe vce(robust)
estimates store d2nd_`var'_degrad
test
}

**(B.2) naive: AOD
xtreg aod $xlist_2nd_degrad if d_miss_aod == 0 & ano <=2016, fe vce(robust)
estimates store d2nd_aod_degrad
test


**(C) threat to first stage
xtreg aod $xlist_2nd_degrad if d_miss_aod == 0 & ano <=2016, fe vce(robust)
estimates store d2nd_aod_degrad
test


*****13/10/2022 non-linear
gen focos_up_sq = focos_up^2
gen focos_up_past = focos_up*area__pasture

global xlist_1st_nlin focos_up focos_nup focos_up_sq ppt evapo area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area d2m t2m w s pib pop ts_dist_sede_100milhab ts_dist_roads d_ano2-d_ano12 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9

global xlist_1st_degrad_nlin focos_up focos_nup focos_up_sq ppt evapo area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area degrad_area d2m t2m w s pib pop ts_dist_sede_100milhab ts_dist_roads d_ano2-d_ano9 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9

***(A.2) R&V full
xtreg focos_for $xlist_1st_nlin if d_miss == 0, fe vce(robust)
estimates store d1st_rv_nlin
test

***(A.2) R&V degrad
xtreg focos_for $xlist_1st_degrad_nlin if d_miss == 0 & ano <=2016, fe vce(robust)
estimates store d1st_rv_nlin_deg
test

**check first derivative > 0
summ focos_up
dis 0.0688 + 2*(-0.0000490)*2408


***Selected model was: squared AF_up, full sample
xtreg focos_for $xlist_1st_nlin if d_miss == 0, fe vce(robust)

***predicted change (fact(ual), count(erfactual))
capture drop ff_fact ff_count focos_up_store fx
predict ff_fact
gen focos_up_store = focos_up
tabstat t2m humidity_redown, s(p50)
replace focos_up = 0 if (t2m > 298.7657 & t2m !=.) & (humidity_redown < .0178116 & humidity_redown !=.)
predict ff_count

*Total avoided FFs
gen fx = ff_fact - ff_count
tabstat fx focos_for if ano == 2019, s(sum)
dis 14554.16/68524
*OBS: fx / focos_for is the reduction rate
codebook fx
tabstat fx, s(p90 p95 p99)
count if fx >= 1
*76,868
dis 76868 / _N
count if fx > 0
*189,507 (many fractionary predictions < 1)

*count by pixel (number of pixels with at least one FF avoided in at least one month)
gen fx_1 = 0
replace fx_1 = fx if fx >= 1
*help egen: total is the one as treats miss as zero
bysort fid_: egen m = total(fx_1) 
count if m >= 1 & fx !=.
dis 1348128 / 144
*9362
dis 1348128 / _N
*18%
*Other way to do the calc, but redundant as both num and denom divided by 144
*dis _N / 144
*dis 9362 / 51307
*2019 share
count if m >= 1 & fx !=. & ano == 2019
count if ano == 2019
*9362
dis 112344 / 615684
*total ppt
tabstat ppt if ano == 2019, s(sum) format(%20.5g)
help tabstat


gen ppt_fx = fx*(-9.908)/ppt
codebook ppt_fx ppt fx
count if ppt_fx < - 1
*13,668
count if ppt_fx < -0.0013 & ppt_fx !=. & ano == 2019
count if ppt_fx !=. & ano == 2019
dis 11068/611842
*1.8% with fx < 0.13%

count if ppt_fx < -0.5 & ppt_fx !=. & ano == 2019

*only pixels with AFs
bysort fid_: egen sum_focos_up_2019 = total(focos_up) if ano == 2019
tabstat fx focos_for if ano == 2019  & sum_focos_up_2019 >= 1, s(sum)
dis 14554.16/68524
*.21239507
tabstat ppt if ano == 2019  & sum_focos_up_2019 >= 1, s(sum) format(%20.5g)


codebook sum_focos_up_2019
*< 0.5%
count if ppt_fx < -0.005 & ppt_fx !=. & ano == 2019 & sum_focos_up_2019 >= 1 
count if ppt_fx !=. & ano == 2019 & sum_focos_up_2019 >= 1 
dis 1019/19464
*5%

*< 5%
count if ppt_fx < -0.05 & ppt_fx !=. & ano == 2019 & sum_focos_up_2019 >= 1 
count if ppt_fx !=. & ano == 2019 & sum_focos_up_2019 >= 1 
dis 580/19464
*5%

codebook ppt_fx if ano == 2019 & sum_focos_up_2019 >= 1 
codebook fx if ano == 2019 & sum_focos_up_2019 >= 1 
count if ano == 2019 & sum_focos_up_2019 >= 1 & (t2m > 298.7657 & t2m !=.) & (humidity_redown < .0178116 & humidity_redown !=.)

**Expo
outsheet fid_ fx ano mes using "sim_AFsq_13_10.txt" if ano == 2019, delimiter("#") replace
outsheet fid_ ppt_fx ano mes using "sim_2nd_stg_AFsq_13_10.txt" if ano == 2019, delimiter("#") replace

***05/11 count mun with ppt_fx < 0

drop c
bysort geocodig_m ano mes: gen c = _n
*bysort geocodig_m: egen m = total(c)
*count if ppt_fx < 0 and ano == 2019


**11h36: redoing the right way by summing both ff avoided and ppt
*totalizing at the key scale: mun-year
*remind: gen ppt_fx = fx*(-9.908)/ppt
bysort geocodig_m ano: egen mun_fx = total(fx)
bysort geocodig_m ano: egen mun_ppt = total(ppt)
drop mun_ppt_fx
gen mun_ppt_fx = mun_fx*(-9.908) / mun_ppt 
browse geocodig_m fid_ c ano mes mun_fx mun_ppt mun_ppt_fx if ano == 2019
*ok, changes only across geocodig_m (remains the same across months)

*(1) muns with increased ppt
count if mun_ppt_fx < 0 & ano == 2019 & mes==1 & c ==1
count if ano == 2019 & mes ==1 & c ==1
dis 623/771
*80% of muns encapsulated
*browse fid_ geocodig_m if ano == 2019 & mes == 1 & c ==1

*(2) muns with increased ppt in at least 1%
count if mun_ppt_fx <= -0.01 & ano == 2019 & mes==1 & c ==1
count if ano == 2019 & mes ==1 & c ==1
dis 17/771
*2%

*(3) muns with increased ppt in at least 0.1%
count if mun_ppt_fx < -0.001 & ano == 2019 & mes==1 & c ==1
count if ano == 2019 & mes ==1 & c ==1
dis 207/771
*27%
 

browse fid_ geocodig_m c if ano == 2019 & mes == 1


table c
duplicates report geocodig_m
help _n
browse fid_ geocodig_m if c == 1

tabstat c if ppt_fx < 0 and ano == 2019, s(sum)

*count of focos_up zeroed

gen focos_up_store = focos_up
tabstat t2m humidity_redown, s(p50)
tabstat focos_up_store if ano==2019 & (t2m > 298.7657 & t2m !=.) & (humidity_redown <.0178116 & humidity_redown !=.), s(sum)
tabstat focos_up_store if ano==2019, s(sum)
dis 211445/232068
*91% of fires up (that's surprising) [but due to high correl between low hum and high temp with fire occurrence]

count if ano==2019 & (t2m > 298.7657 & t2m !=.) & (humidity_redown <.0178116 & humidity_redown !=.)
count if ano==2019
dis 198677/615684
*32% of pixels



***Rob test: weather squared

foreach var of varlist ppt evapo humidity_redown d2m t2m w s {
gen `var'_sq = `var' ^2 
}

global xlist_1st_wsq focos_up focos_nup ppt evapo area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area d2m t2m w s pib pop ts_dist_sede_100milhab ts_dist_roads d_ano2-d_ano12 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9 ppt_sq evapo_sq humidity_redown_sq d2m_sq t2m_sq w_sq s_sq

global xlist_1st_degrad_wsq focos_up focos_nup ppt evapo area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area degrad_area d2m t2m w s pib pop ts_dist_sede_100milhab ts_dist_roads d_ano2-d_ano9 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9 ppt_sq evapo_sq humidity_redown_sq d2m_sq t2m_sq w_sq s_sq


***(A.2) R&V, no degrad
xtreg focos_for $xlist_1st_wsq if d_miss == 0, fe vce(robust)
estimates store d1st_rv_sq
test

***(A.2) R&V, degrad
xtreg focos_for $xlist_1st_degrad_wsq if d_miss == 0 & ano <=2016, fe vce(robust)
estimates store d1st_rv_degrad_sq
test



****ESTTAB
esttab d1st_naive_nb d1st_rv_nb d1st_naive_nb_deg d1st_rv_nb_deg using "test_grid_[put_date_here]_1st_nonlin_first_NB.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace

*13/10/22 no past int
esttab d1st_rv_nlin d1st_rv_nlin_deg using "test_grid_[put_date_here]1st_nonlin_first_onlyAFsq.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace


*01/11/22 squared weather
esttab d1st_rv_sq d1st_rv_degrad_sq using "test_grid_[put_date_here]_1st_nonlin_first_weather_sq_rob.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace

***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE
***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE
***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE***MATCHING 2nd STAGE

use "(dataset)fire_degrad_rain_1st_2nd_subpaths.dta", clear

count if focos_for > 0 & focos_for !=. & d_miss == 0
count if d_ff ==1 & d_miss == 0

*duplicates report 
*drop because oo
drop if ano > 2016 | d_miss > 0
collapse (sum) d_ff d_miss (mean) ppt focos_up focos_nup focos_agr area__for area__crops area__soy area__pasture area__grassland area__urban area__water defor_area degrad_area  slope_pct  humidity_redown d2m t2m w s cloud pib pop dist_sede_100milhab dist_roads d_uf2-d_uf9, by(fid_)


global xlist_1st focos_up focos_nup focos_agr area__for area__crops area__soy area__pasture area__grassland area__urban area__water defor_area degrad_area  slope_pct  humidity_redown d2m t2m w s cloud pib pop dist_sede_100milhab dist_roads  d_uf2-d_uf9

**without UFs
global xlist_1st_noufs focos_up focos_nup focos_agr area__for area__crops area__soy area__pasture area__grassland area__urban area__water defor_area degrad_area  slope_pct  humidity_redown d2m t2m w s cloud pib pop dist_sede_100milhab dist_roads

*redefine d_ff (Y)
*drop d_ff_bin
gen d_ff_bin = 0
replace d_ff_bin = 1 if d_ff > 0 & d_ff !=.

summ d_ff_bin $xlist_1st

*count of treated and with AFs
count if d_ff > 0
count if focos_agr > 0

**20/10
gen area__agric = area__crops + area__pasture
tabstat area__for area__agric if d_ff_bin == 0, s(sum)
tabstat area__for area__agric if d_ff_bin == 1, s(sum)
*0/1 ratios 1: for
dis 2145111/ 1722280
*1.24
*0/1 ratios 1: agr
dis 68876.22/632922.1
*0.1

*drop area__tot
*drop area__agric_s
gen area__tot = area__for + area__crops + area__pasture + area__grassland + area__urban + area__water
gen area__agric_s = area__agric / area__tot
gen area__for_s = area__for / area__tot
codebook area__agric_s area__for_s
*50% forest is below the 25% pctile, so will be the cutoff
count if area__for_s < 0.5 | area__for_s==.
drop if area__for_s <0.5 | area__for_s==.
*PS: ppt must be kept because mahalanobis matching requires outcome

*******(1) psmatch with Mahalanobis distance without caliper [matching first stage]*******(1) psmatch with Mahalanobis distance without caliper [matching first stage]
**(1.a) Estimation
psmatch2 d_ff_bin, outcome(ppt) mahalanobis($xlist_1st) ai(1) altvariance
*psmatch2 d_ff_bin, outcome(ppt) mahalanobis($xlist_1st_noufs) ai(1) altvariance

count if _n1 !=.


****(1.b) Rubin table [bef vs af]
pstest $xlist_1st, both t(d_ff_bin) mweight(_weight) label scatter

*pstest $xlist_1st_noufs, both t(d_ff_bin) mweight(_weight) label scatter


***(1.c) check of matching quality with support and balance graphs [bef vs af]***(5) check of matching quality with support and balance graphs [bef vs af]

**(1.c.A) gráfico suporte comum: before matching
regress d_ff_bin $xlist_1st, vce(robust)
*drop p_hat
predict double p_hat
codebook p_hat


*count if p_hat < 0 | (p_hat > 1 & p_hat !=.)
*count if d_miss == 0
*dis 6652/50894
*13%: this is small enough for probit/logit to be not needed

psgraph, pscore(p_hat)
graph export supp_bef_nocal_dist_20_10_22.tif, name(Graph) width(800) height(600) replace


**exporting to R
cd "D:\Pesquisa\Pesquisa_2022\water_pp"
outsheet fid_ d_ff_bin _id _n1 $xlist_1st using "matching_nocal_res_dist_20_10_22.txt", delimiter("#") replace
*drop _merge_filter
*drop matchres
*joinby _id using matching_cont_filter_dist_02_11, unmatched(both) _merge(_merge_filter)
joinby _id using matching_cont_filter_dist_20_10_22, unmatched(both) _merge(_merge_filter)
*note: in 12/11 generated other file but is the same as 02_11
*12/11 10h57: added count of times control was used in "count"

tab _merge_filter
browse d_ff_bin _id _n1 matchres ctr_count

*non-comparable controls are d_ff_bin == 0 & matchres == 0
*so it is only to use matchres ==. | matchres == 1

*count of controls
replace ctr_count ="." if ctr_count =="NA"
destring ctr_count, replace
browse d_ff_bin _id _n1 matchres ctr_count


*gráfico suporte comum: after matching
**p_hat_af
*drop p_hat_af
gen p_hat_af =. 
replace p_hat_af = p_hat if matchres==.|matchres==1
*codebook p_hat_af
count if p_hat_af < 0 | (p_hat_af > 1 & p_hat_af !=.)

psgraph, pscore(p_hat_af)
graph export supp_af_nocal_dist_20_10_22.tif, name(Graph) width(800) height(600) replace


***(1.c.B) gráfico balanceamento
// compare _pscores before matching & save graph to disk
//27/08: I replaced _pscore to p_hat estimated from probit [RETURN to check it is the same]
//also replaced, in xtitle, "propensity score" to "treatment probabilities"
twoway (kdensity p_hat if _treated==1) (kdensity p_hat if _treated==0, ///
lpattern(dash)), legend( label( 1 "treated") label( 2 "control" ) ) ///
xtitle("treatment probabilities BEFORE matching") saving(before_nocal_dist_12_11, replace)

twoway (kdensity p_hat_af if _treated==1) (kdensity p_hat_af if _treated==0 ///
& matchres==1, lpattern(dash)), legend( label( 1 "treated") label( 2 "control" )) ///
xtitle("treatment probabilities AFTER matching") saving(after_nocal_dist_12_11, replace)

graph combine before_nocal_dist_12_11.gph after_nocal_dist_12_11.gph, ycommon 
graph export bal_no_cal_dist_20_10_22.tif, name(Graph) width(800) height(600) replace



*******(2) psmatch with Mahalanobis distance with caliper [matching second stage]*******(2) psmatch with Mahalanobis distance with caliper [matching second stage]
**(2.a) Estimation
*[RETURN] one to one matching may not save pscore, in this case it there will be need to activate line below
*probit d_pa $xlist, vce(robust)
*predict double p_hat

*12/11: _pscore-based caliper now replaced by Arriagada and Ferraro sd-covar-diff-based

***Generate pair ids here

***check if _weight is not distorted by manual sd adjustment
codebook _weight
count if _weight > 1 & _weight !=. & _treated == 0
count if _weight > 1 & _weight !=. & _treated == 0 & matchres == 1

tabstat _weight if _treated == 0 & matchres == 0, s(min max)
*all miss
tabstat _weight if _treated == 0 & matchres == 1, s(min max)
*range = 1-421
tabstat _weight if _treated == 1, s(min max)
*all one

gen check = _weight - ctr_count
codebook check
codebook check if _treated == 0 & matchres == 1
*this is proof that _weight = number of times each control is used.


*if weight = 1, no distortion
 
drop matchres
drop _merge_filter
drop ctr_count
drop _merge_cont_filter

******1SD CALIPER

*treat filter for _support
joinby fid_ using "matching_treat_filter_cov1sd_20_10_22", unmatched(both) _merge(_merge_treat_filter)
table _merge_treat_filter
*ok, only treated
tab d_ff_bin

*control filter for graphs and _weight
joinby _id using "matching_cont_filter_cov1sd_cal_20_10_22", unmatched(both) _merge(_merge_cont_filter)
table _merge_cont_filter
*ok, only control

***Adjustment of key matching output variables: _support and _weight
*for pstest, psgraph and balance graph
**(1) _support
tab d_pair_out if _treated == 1

gen _support_str = _support
tab _support
*100% in support
*replace _support = "Off support" if d_pair_out == 1 & _treated == 1
*support = 1/0
replace _support = 0 if d_pair_out == 1 & _treated == 1
tab _support if _treated == 1
tab d_pair_out if _treated == 1
*ok


**(2) _weight
gen _weight_str = _weight
replace ctr_count ="." if ctr_count=="NA"
destring ctr_count, replace
replace _weight = ctr_count
codebook _weight if d_ff_bin == 0
*looks ok

****(2.b) Rubin table
*Note: with the option "both" it compares before and after
*teste de viés
*Doubt: after altering _support of treated, probably the test will ignore the excluded treated units
*but only in M rows
pstest $xlist_1st, both t(d_ff_bin) mweight(_weight) label scatter
*Doubt: if not, use:
*pstest $xlist_1st if d_pair_out == 0, both t(d_ff_bin) mweight(_weight) label scatter
*don't do that as _weights get distorted and so the mean-diff test

****(2.c.A) Common support graph
psgraph, pscore(p_hat)
graph export supp_bef_cal_cov1sd_dist_20_10_22.tif, name(Graph) width(800) height(600) replace

*gráfico suporte comum: after matching
**exporting to R
*outsheet fid_ d_ff_bin _id _n1 using "matching_cal_res_dist_02_11.txt", delimiter("#")

**p_hat_af
drop p_hat_af
gen p_hat_af =. 
replace p_hat_af = p_hat if matchres==.|matchres==1

psgraph, pscore(p_hat_af)
graph export supp_af_cal_cov1sd_dist_20_10_22.tif, name(Graph) width(800) height(600) replace


****(2.c.B) Common support graph
// compare _pscores before matching & save graph to disk
//27/08: I replaced _pscore to p_hat estimated from probit [RETURN to check it is the same]
//also replaced, in xtitle, "propensity score" to "treatment probabilities"

**Note: should count only treated on support
*Treated on support: 13,783
tab _support if _treated == 1
*Comparable controls: 5,105
tab matchres if _treated == 0

dis 13783+5105
count if (matchres==. & _support==1)|matchres==1
*16,112 ok

*01/11: important note: this graph includes all treated, even out of support so must subsel

twoway (kdensity p_hat if _treated==1) (kdensity p_hat if _treated==0, ///
lpattern(dash)), legend( label( 1 "treated") label( 2 "control" ) ) ///
xtitle("treatment probabilities BEFORE matching") saving(before_cal_cov1sd_dist_15_11, replace)
// compare _pscores *after* matching & save graph to disk

*important: when there are treated out of support, must apply _support == 1 below
*as generation of p_hat_af did not account for that (in this sense matchres == 1 is redundant below)

twoway (kdensity p_hat_af if _treated==1 & _support==1) (kdensity p_hat_af if _treated==0 ///
& matchres==1, lpattern(dash)), legend( label( 1 "treated") label( 2 "control" )) ///
xtitle("treatment probabilities AFTER matching") saving(after_cal_cov1sd_dist_15_11, replace)


// combine these two graphs that were saved to disk
// put both graphs on y axes with common scales
graph combine before_cal_cov1sd_dist_15_11.gph after_cal_cov1sd_dist_15_11.gph, ycommon 
graph export bal_cov1sd_cal_dist_20_10_22.tif, name(Graph) width(800) height(600) replace

**saving whole information on matching
*outsheet fid_ _id _n1 _support _weight d_ff_bin matchres using "match_half_cal__cov1sd_dist_info_15_11.txt", delimiter("#") replace



******0.5SD CALIPER
*deleting all vars from 1sd cal
drop _merge_treat_filter-ctr_control p_hat_af

*treat filter for _support
*cd "E:\Dell_HD\water_pp"
joinby fid_ using "matching_treat_filter_cov0_5sd_20_10_22", unmatched(both) _merge(_merge_treat_filter)
table _merge_treat_filter
*ok, only treated
tab d_ff_bin

*control filter for graphs and _weight
joinby _id using "matching_cont_filter_cov0_5sd_cal_20_10_22", unmatched(both) _merge(_merge_cont_filter)
table _merge_cont_filter
*ok, only control

***Adjustment of key matching output variables: _support and _weight
*for pstest, psgraph and balance graph
**(1) _support
tab d_pair_out if _treated == 1

*resuming original _support, will all treated matched
table _support
replace _support = 1
tab _support
*OBS: here I deleted by mistake the _support_str the first time I run 0.5sd

*100% in support
*replace _support = "Off support" if d_pair_out == 1 & _treated == 1
*support = 1/0
replace _support = 0 if d_pair_out == 1 & _treated == 1
tab _support if _treated == 1
tab d_pair_out if _treated == 1
*ok


**(2) _weight
*gen _weight_str = _weight
replace ctr_count ="." if ctr_count=="NA"
destring ctr_count, replace
replace _weight = ctr_count
codebook _weight if d_ff_bin == 0
*looks ok

****(2.b) Rubin table
*Note: with the option "both" it compares before and after
*teste de viés
*Doubt: after altering _support of treated, probably the test will ignore the excluded treated units
*but only in M rows
pstest $xlist_1st, both t(d_ff_bin) mweight(_weight) label scatter
*Doubt: if not, use:
*pstest $xlist_1st if d_pair_out == 0, both t(d_ff_bin) mweight(_weight) label scatter
*don't do that as _weights get distorted and so the mean-diff test

****(2.c.A) Common support graph
psgraph, pscore(p_hat)
graph export supp_bef_cal_covhalfsd_dist_20_10_22.tif, name(Graph) width(800) height(600) replace

*gráfico suporte comum: after matching
**exporting to R
*outsheet fid_ d_ff_bin _id _n1 using "matching_cal_res_dist_02_11.txt", delimiter("#")

**p_hat_af
drop p_hat_af
gen p_hat_af =. 
replace p_hat_af = p_hat if matchres==.|matchres==1

psgraph, pscore(p_hat_af)
graph export supp_af_cal_covhalfsd_dist_20_10_22.tif, name(Graph) width(800) height(600) replace


****(2.c.B) Common support graph
// compare _pscores before matching & save graph to disk
//27/08: I replaced _pscore to p_hat estimated from probit [RETURN to check it is the same]
//also replaced, in xtitle, "propensity score" to "treatment probabilities"

**Note: should count only treated on support
*Treated on support: 6,316
tab _support if _treated == 1
*Comparable controls: 3,769
tab matchres if _treated == 0

dis 6316+3769
count if (matchres==. & _support==1)|matchres==1
*16,112 ok

*01/11: important note: this graph includes all treated, even out of support so must subsel

twoway (kdensity p_hat if _treated==1) (kdensity p_hat if _treated==0, ///
lpattern(dash)), legend( label( 1 "treated") label( 2 "control" ) ) ///
xtitle("treatment probabilities BEFORE matching") saving(before_cal_covhalfsd_dist_15_11, replace)
// compare _pscores *after* matching & save graph to disk

*important: when there are treated out of support, must apply _support == 1 below
*as generation of p_hat_af did not account for that (in this sense matchres == 1 is redundant below)

twoway (kdensity p_hat_af if _treated==1 & _support==1) (kdensity p_hat_af if _treated==0 ///
& matchres==1, lpattern(dash)), legend( label( 1 "treated") label( 2 "control" )) ///
xtitle("treatment probabilities AFTER matching") saving(after_cal_covhalfsd_dist_15_11, replace)


// combine these two graphs that were saved to disk
// put both graphs on y axes with common scales
graph combine before_cal_covhalfsd_dist_15_11.gph after_cal_covhalfsd_dist_15_11.gph, ycommon 
graph export bal_covhalfsd_cal_dist_20_10_22.tif, name(Graph) width(800) height(600) replace

**saving whole information on matching
outsheet fid_ _id _n1 _support _weight d_ff_bin matchres d_pair_out ctr_count using "match_half_cal__covhalfsd_dist_info_20_10_22.txt", delimiter("#") replace



***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING
***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING
***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING***2nd SUBPATH-POSTMATCHING
clear all
cd "your directory here"
insheet using "(dataset)fire_degrad_rain_1st_2nd_subpaths.txt", delimiter("#") clear


global xlist_2nd d_ff area__for area__crops area__soy area__pasture area__grassland area__urban area__water defor_area degrad_area humidity_redown d2m t2m w s cloud d_ano2-d_ano9 d_mes2-d_mes12

global xlist_2nd_nodeg d_ff area__for area__crops area__soy area__pasture area__grassland area__urban area__water defor_area humidity_redown d2m t2m w s cloud d_ano2-d_ano12 d_mes2-d_mes12


***(A) WITH DEGRAD
**Second stage: ppt and evapo
foreach var of varlist ppt evapo {
xtreg `var' $xlist_2nd if d_samp == 1 & ano <= 2016, fe vce(robust)
estimates store d2nd_`var'
test
}

**Second stage: AOD
xtreg aod $xlist_2nd if d_samp == 1 & ano <= 2016, fe vce(robust)
estimates store d2nd_aod
test


***(B) WITHOUT DEGRAD
**Second stage: ppt and evapo
foreach var of varlist ppt evapo {
xtreg `var' $xlist_2nd_nodeg if d_samp == 1, fe vce(robust)
estimates store d2nd_`var'_nodeg
test
}

**Second stage: AOD
xtreg aod $xlist_2nd_nodeg if d_samp == 1, fe vce(robust)
estimates store d2nd_aod_nodeg
test

esttab d2nd_ppt_nodeg d2nd_evapo_nodeg d2nd_aod_nodeg using "test_grid_[put_date_here]_postmatch_nopibpop_nodeg_calhalfsd.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace


****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage
****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage
****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage****THREAT TEST 1st stage

gen d_for_congruent = 0
replace d_for_congruent = 1 if class_angle_for =="congruent"

gen d_for_congruent_no5km = 0
replace d_for_congruent_no5km = 1 if class_angle_for_no5km =="congruent"

Set of controls: (i) precipitation to account for wet/dry season, 
*(ii) distance to roads, (iii) crop area, (iv) deforested area, (v) pib, (vi) pop, (vii) slope, 
*(viii) year, month, state dummies.

global thr_list ppt dist_roads area__crops area__pasture defor_area pib pop ts_slope_pct d_ano2-d_ano12 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9
global thr_list_nodum ppt dist_roads area__crops area__pasture defor_area pib pop ts_slope_pct

***(1) regress
regress focos_agr d_for_congruent if d_miss ==0, vce(robust)
estimates store thr

regress focos_agr d_for_congruent_no5km if d_miss ==0, vce(robust)
estimates store thr_no5km

ttest focos_agr if d_miss ==0, by(class_angle_for)

ttest focos_agr if d_miss ==0, by(class_angle_for_no5km)

***(2) FE
xtreg focos_agr d_for_congruent $thr_list, fe vce(robust)
estimates store thr_fe

xtreg focos_agr d_for_congruent_no5km $thr_list, fe vce(robust)
estimates store thr_fe_no5km

xtreg focos_agr d_for_congruent $thr_list_nodum, fe vce(robust)
estimates store thr_fe_nd

xtreg focos_agr d_for_congruent_no5km $thr_list_nodum, fe vce(robust)
estimates store thr_fe_no5km_nd



esttab thr thr_no5km thr_fe thr_fe_no5km thr_fe_nd thr_fe_no5km_nd using "test_grid_[put_date_here]_threat_test.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps mtitles star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace



**** robustness**** robustness**** robustness**** robustness**** robustness**** robustness**** robustness**** robustness
**** robustness**** robustness**** robustness**** robustness**** robustness**** robustness**** robustness**** robustness


*****Climate change 1
global xlist_1st focos_up focos_nup ppt evapo area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area d2m t2m w s pib pop ts_dist_sede_100milhab ts_dist_roads d_ano2-d_ano12 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9

***(R.0) R&V full sample, all years
xtreg focos_for $xlist_1st if d_miss == 0, fe vce(robust)
estimates store rob_eyal0
*mat b_full = e(b)
*matlist b_full
scalar b_full = e(b)[1,1]
*mat s_full = e(V)
*matlist s_full
scalar s_full = (e(V)[1,1])^(1/2)
scalar n_full = e(N)
*scalar list

test

***(R.1) R&V full sample, first 4 years
xtreg focos_for $xlist_1st if d_miss == 0 & ano <= 2011, fe vce(robust)
estimates store rob_eyal1
*mat b_f4 = e(b)
*matlist b_f4
scalar b_f4 = e(b)[1,1]
*mat s_f4 = e(V)
*matlist s_f4
scalar s_f4 = (e(V)[1,1])^(1/2)
scalar n_f4 = e(N)
scalar list

test


***(A.2) R&V full sample, last 4 years
xtreg focos_for $xlist_1st if d_miss == 0 & ano >= 2016, fe vce(robust)
estimates store rob_eyal1_last

*mat b_l4 = e(b)
*matlist b_l4
scalar b_l4 = e(b)[1,1]
*mat s_l4 = e(V)
*matlist s_l4
scalar s_l4 = (e(V)[1,1])^(1/2)
scalar n_l4 = e(N)
scalar list
help xtreg


test

**f4 vs full
scalar s_a_f4_full = ((n_f4-1)*((s_f4)^2) + (n_full-1)*((s_full)^2))/(n_f4+n_full-2)
scalar estat_f4_full = (b_f4 - b_full)/((s_a_f4_full^(1/2))*((1/n_f4+1/n_full)^(1/2)))
scalar pval_f4_full = t(n_f4+n_full-2,-abs(estat_f4_full))
dis estat_f4_full
dis pval_f4_full

**l4 vs full
scalar s_a_l4_full = ((n_l4-1)*((s_l4)^2) + (n_full-1)*((s_full)^2))/(n_l4+n_full-2)
scalar estat_l4_full = (b_l4 - b_full)/((s_a_l4_full^(1/2))*((1/n_l4+1/n_full)^(1/2)))
scalar pval_l4_full = t(n_l4+n_full-2,-abs(estat_l4_full))
dis estat_l4_full
dis pval_l4_full

**f4 vs l4
scalar s_a_f4_l4 = ((n_f4-1)*((s_f4)^2) + (n_l4-1)*((s_l4)^2))/(n_f4+n_l4-2)
scalar estat_f4_l4 = (b_f4 - b_l4)/((s_a_f4_l4^(1/2))*((1/n_f4+1/n_l4)^(1/2)))
scalar pval_f4_l4 = t(n_f4+n_l4-2,-abs(estat_f4_l4))
dis estat_f4_l4
dis pval_f4_l4

dis ((n_f4-1) / (n_f4+n_l4-2))*((s_f4)^2) + ((n_l4-1) / (n_f4+n_l4-2))*((s_l4)^2)
mat matstat = [estat_f4_full, estat_l4_full, estat_f4_l4] \ [pval_f4_full, pval_l4_full, pval_f4_l4]
matlist matstat

scalar list
*check OMs of s_a
dis s_a_f4_l4 ^(1/2)
dis s_a_f4_full ^(1/2)
dis s_a_l4_full ^(1/2)
*ok, 10^-3, same as ses
dis b_f4/b_full-1
dis b_l4/b_full-1
dis b_f4/b_l4-1



*13/10/22 no past int
esttab rob_eyal0 rob_eyal1 rob_eyal1_last using "bioecon_rob_CC1.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace



*********Climate change 2

collapse (mean) t2m humidity_redown, by(fid_ ano)
reshape wide t2m humidity_redown, i(fid_) j(ano)
gen t2m_r = t2m2019 / t2m2008
gen hum_r = humidity_redown2019 / humidity_redown2008
tabstat t2m_r hum_r, s(p25 p75)

*t2m d
gen d_t2m_big =0 
*> p75
replace d_t2m_big = 1 if t2m_r >  1.002861
gen d_t2m_small =0 
*< p25
replace d_t2m_small = 1 if t2m_r <  1.001528

*hum_redown d (p25 < 1, so fall)
gen d_hum_big =0 
*<p25
replace d_hum_big = 1 if hum_r < .9918307
gen d_hum_small =0 
*>p75
replace d_hum_small = 1 if hum_r >  1.037275

sum d_t2m_big d_t2m_small d_hum_big d_hum_small
*ok, all near 25%
codebook d_t2m_big d_t2m_small d_hum_big d_hum_small

keep fid_ d_t2m_big d_t2m_small d_hum_big d_hum_small
save class_cc_experiences_24_10, replace

use "(dataset)fire_degrad_rain_1st_2nd_subpaths.dta", clear
joinby fid_ using class_cc_experiences_24_10, unmatched(both) _merge(_merge_class_cc_exp)
tab _merge_class_cc_exp


*****estimations
global xlist_1st focos_up focos_nup ppt evapo area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area d2m t2m w s pib pop ts_dist_sede_100milhab ts_dist_roads d_ano2-d_ano12 d_mes2-d_mes12 ts_d_uf2-ts_d_uf9

***(R.0) R&V full sample, all years
xtreg focos_for $xlist_1st if d_miss == 0, fe vce(robust)
estimates store rob_eyal0
*mat b_full = e(b)
*matlist b_full
scalar b_full = e(b)[1,1]
*mat s_full = e(V)
*matlist s_full
scalar s_full = (e(V)[1,1])^(1/2)
scalar n_full = e(N)
*scalar list

test

***(R.1) R&V full sample, big change
xtreg focos_for $xlist_1st if d_miss == 0 & d_t2m_big == 1 & d_hum_big == 1, fe vce(robust)
estimates store rob_eyal2_big
*mat b_big = e(b)
*matlist b_big
scalar b_big = e(b)[1,1]
*mat s_big = e(V)
*matlist s_big
scalar s_big = (e(V)[1,1])^(1/2)
scalar n_big = e(N)
scalar list

test


***(A.2) R&V full sample, small change
xtreg focos_for $xlist_1st if d_miss == 0 & d_t2m_small == 1 & d_hum_small == 1, fe vce(robust)
estimates store rob_eyal2_small

*mat b_small = e(b)
*matlist b_small
scalar b_small = e(b)[1,1]
*mat s_small = e(V)
*matlist s_small
scalar s_small = (e(V)[1,1])^(1/2)
scalar n_small = e(N)
scalar list
help xtreg


test

**f4 vs full
scalar s_a_big_full = ((n_big-1)*((s_big)^2) + (n_full-1)*((s_full)^2))/(n_big+n_full-2)
scalar estat_big_full = (b_big - b_full)/((s_a_big_full^(1/2))*((1/n_big+1/n_full)^(1/2)))
scalar pval_big_full = t(n_big+n_full-2,-abs(estat_big_full))
dis estat_big_full
dis pval_big_full

**l4 vs full
scalar s_a_small_full = ((n_small-1)*((s_small)^2) + (n_full-1)*((s_full)^2))/(n_small+n_full-2)
scalar estat_small_full = (b_small - b_full)/((s_a_small_full^(1/2))*((1/n_small+1/n_full)^(1/2)))
scalar pval_small_full = t(n_small+n_full-2,-abs(estat_small_full))
dis estat_small_full
dis pval_small_full

**f4 vs l4
scalar s_a_big_small = ((n_big-1)*((s_big)^2) + (n_small-1)*((s_small)^2))/(n_big+n_small-2)
scalar estat_big_small = (b_big - b_small)/((s_a_big_small^(1/2))*((1/n_big+1/n_small)^(1/2)))
scalar pval_big_small = t(n_big+n_small-2,-abs(estat_big_small))
dis estat_big_small
dis pval_big_small

dis ((n_big-1) / (n_big+n_small-2))*((s_big)^2) + ((n_small-1) / (n_big+n_small-2))*((s_small)^2)

capture drop matstat
mat matstat = [estat_big_full, estat_small_full, estat_big_small] \ [pval_big_full, pval_small_full, pval_big_small]
matlist matstat

scalar list
*check OMs of s_a
dis s_a_f4_l4 ^(1/2)
dis s_a_f4_full ^(1/2)
dis s_a_l4_full ^(1/2)
*ok, 10^-3, same as ses
dis b_f4/b_full-1
dis b_l4/b_full-1
dis b_f4/b_l4-1
*cleary there's sig diff for f4 vs full and f4 vs l4 (coeff %change >= 13%)


*13/10/22 no past int
esttab rob_eyal0 rob_eyal2_big rob_eyal2_small using "bioecon_rob_Eyal2_24_10.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace


********Robustness 2: squared weather
clear all
*cd "G:\water_pp\Comet\"
cd "D:\Pesquisa\Pesquisa_2022\water_pp"
*log using "log_stata_17_08_11_8h.txt", text
use "(dataset)fire_degrad_rain_1st_2nd_subpaths.dta", clear

foreach var of varlist humidity_redown d2m t2m w s {
gen `var'_sq = `var' ^2 
}

*without pib and pop [01/11/22]
global xlist_2nd d_ff area__for area__crops area__soy area__pasture area__grassland area__urban area__water defor_area degrad_area humidity_redown d2m t2m w s cloud d_ano2-d_ano9 d_mes2-d_mes12 humidity_redown_sq d2m_sq t2m_sq w_sq s_sq


*without pib and pop
global xlist_2nd_nodeg d_ff area__for area__crops area__soy area__pasture area__grassland area__urban area__water defor_area humidity_redown d2m t2m w s cloud d_ano2-d_ano12 d_mes2-d_mes12 humidity_redown_sq d2m_sq t2m_sq w_sq s_sq


***(A) WITH DEGRAD
**Second stage: ppt and evapo
foreach var of varlist ppt evapo {
xtreg `var' $xlist_2nd if d_samp == 1 & ano <= 2016, fe vce(robust)
estimates store d2nd_`var'_wsq
test
}

**Second stage: AOD
xtreg aod $xlist_2nd if d_samp == 1 & ano <= 2016, fe vce(robust)
estimates store d2nd_aod_wsq
test


***(B) WITHOUT DEGRAD
**Second stage: ppt and evapo
foreach var of varlist ppt evapo {
xtreg `var' $xlist_2nd_nodeg if d_samp == 1, fe vce(robust)
estimates store d2nd_`var'_wsq_nodeg
test
}

**Second stage: AOD
xtreg aod $xlist_2nd_nodeg if d_samp == 1, fe vce(robust)
estimates store d2nd_aod_wsq_nodeg
test

esttab d2nd_ppt_wsq d2nd_evapo_wsq d2nd_aod_wsq d2nd_ppt_wsq_nodeg d2nd_evapo_wsq_nodeg d2nd_aod_wsq_nodeg using "test_grid_[put_date_here]_postmatch_weather_sq_rob.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace

***03/11 SDs
summ focos_up if d_miss == 0

*how likely the shock is: half sd extreme count
summ focos_up if d_miss == 0
return list
scalar mb = r(mean) -  r(sd)
scalar ma = r(mean) +  r(sd)

count if d_miss ==0 & focos_up !=. & (focos_up < 0.5*mb | focos_up > 0.5*ma)
count if d_miss ==0 
dis 92482/7328736
*1.26%

*how likely the shock is: half sd extreme count: deg
summ focos_up if d_miss == 0 & ano <=2016
return list
scalar mbd = r(mean) -  r(sd)
scalar mad = r(mean) +  r(sd)

count if d_miss ==0 & ano <=2016 & focos_up !=. & (focos_up < 0.5*mbd | focos_up > 0.5*mad)
count if d_miss ==0 & ano <=2016

dis 70241/5496552
*1.23%
dis r(mean) - r(sd)

*ext drought year dummy
summ d_ano9 if d_miss == 0

summ d_mes6 if d_miss == 0

***
***(B) WITHOUT DEGRAD
*03/11 without pib and pop
*Base = may
global xlist_2nd_nodeg d_ff area__for area__crops area__soy area__pasture area__grassland area__urban area__water defor_area humidity_redown d2m t2m w s cloud d_ano2-d_ano12 d_mes1-d_mes4 d_mes6-d_mes12

*Base = Aug
*global xlist_2nd_nodeg d_ff area__for area__crops area__soy area__pasture area__grassland area__urban area__water defor_area humidity_redown d2m t2m w s cloud d_ano2-d_ano12 d_mes1-d_mes7 d_mes9-d_mes12


**Second stage: ppt and evapo
foreach var of varlist ppt {
xtreg `var' $xlist_2nd_nodeg if d_samp == 1, fe vce(robust)
estimates store d2nd_`var'_nodeg
test
}

 
***(A.2) R&V
global xlist_1st focos_up focos_nup ppt evapo area__for area__crops area__soy area__pasture area__grassland area__urban area__water humidity_redown ts_slope_pct defor_area d2m t2m w s pib pop ts_dist_sede_100milhab ts_dist_roads d_ano2-d_ano12 d_mes1-d_mes7 d_mes9-d_mes12 ts_d_uf2-ts_d_uf9

xtreg focos_for $xlist_1st if d_miss == 0, fe vce(robust)
estimates store d1st_rv
test

table mes if d_miss ==0, stat(mean ppt)
table mes if d_miss ==0, stat(mean humidity_redown)
table mes if d_miss ==0, stat(mean t2m)



***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH
***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH
***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH***3rd SUBPATH

clear all
cd "your directory here"
insheet using "(dataset)fire_degrad_rain_3rd.dta"

***panel set
xtset geocodig_m time


***Variable lists
**no log vs
global xlist_month price wage_all mean_ppt_1-mean_ppt_12 sur_ppt_1-sur_ppt_12 t2m_wet t2m_dry solar ts_slope_pct25 ts_slope_pct50 ts_slope_pct75 ts_pca_soil ts_dist_roads ts_dist_sede_100milhab ts_fertilizer ts_tractor ts_pesticide ts_urban ts_area_tot area_for area_pasture area_crops evapo humidity ts_cot ts_csv ts_mze ts_soy ts_d_uf_*


global xlist_month_wagr price wage_agr mean_ppt_1-mean_ppt_12 sur_ppt_1-sur_ppt_12 t2m_wet t2m_dry solar ts_slope_pct25 ts_slope_pct50 ts_slope_pct75 ts_pca_soil ts_dist_roads ts_dist_sede_100milhab ts_fertilizer ts_tractor ts_pesticide ts_urban ts_area_tot area_for area_pasture area_crops evapo humidity ts_cot ts_csv ts_mze ts_soy ts_d_uf_*


*nolog
sum vaa $xlist_month i.time if d_miss == 0
tab d_miss uf
*RO and TO out of sample

sum vaa $xlist_month_wagr i.time if d_miss_wagr == 0
tab d_miss_wagr uf
*browse


****MODELS: GDP
*FE Cobb: wage of all sectors
xtreg vaa $xlist_month i.time if d_miss == 0, fe vce(cluster geocodig_m)
estimates store FE_Cobb_all
test
*both nonsig
margin, dydx(sur_ppt_8) atmeans

*FE Cobb: wage only agriculture
xtreg vaa $xlist_month_wagr i.time if d_miss_wagr == 0, fe vce(cluster geocodig_m)
estimates store FE_Cobb_wagr
*both nonsig
margin, dydx(sur_ppt_8) atmeans
test

*GDP effect
dis 493/43758.37 

esttab FE_Cobb_all FE_Cobb_wagr using "water_3rd_[put_date_here].txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace


*(01/11) check if double miss leaves no attrition
count if d_miss == 0 & d_miss_temp == 0
dis 8640/12
*720, ok

****MODELS: output (wage of all sectors only) [re-run 01/11 without attrition]
**(a) Temp

foreach y in  "q_cott" "q_cass" "q_maiz" "q_soyb" {         
xtreg `y' $xlist_month i.time if d_miss ==0 & d_miss_temp == 0, fe vce(cluster geocodig_m)
estimates store temp_`y'
test

***storing year-total surprise rain: start
mat beta = e(b)
mat var = vecdiag(e(V))
scalar df = e(N) - 68
*71 vars (with intercept) - last ts*uf -time08 - time19
*CHECK DF
capture mat drop pval
forvalues i = 15(1)26 {
matrix pval = nullmat(pval) \ 2*t(df,-abs(beta[1,`i']/sqrt(var[1,`i'])))
}
scalar sum_sur = 0
forvalues i = 1(1)12 {
if (pval[`i',1] <= 0.05) {
*matrix ex = nullmat(ex) \ `i'
scalar sum_sur = sum_sur + beta[1,14+`i']
}
}
matrix ex_temp = nullmat(ex_temp) \ sum_sur
mat drop beta var pval
scalar drop df sum_sur
***storing year-total surprise rain: end
}
*sum checked for the 4 crops in excel: ok


*(01/11) check if double miss leaves no attrition
count if d_miss == 0 & d_miss_pere == 0
dis 7140/12
*595, ok


**(a) Pere
foreach y in  "q_bana" "q_coco" "q_coff" "q_palm" "q_pepp" {         
xtreg `y' $xlist_month i.time if d_miss ==0 & d_miss_pere == 0, fe vce(cluster geocodig_m)
estimates store pere_`y'
test

***storing year-total surprise rain: start
mat beta = e(b)
mat var = vecdiag(e(V))
scalar df = e(N) - 68
*71 vars (with intercept) - last ts*uf -time08 - time19
*CHECK DF
capture mat drop pval
forvalues i = 15(1)26 {
matrix pval = nullmat(pval) \ 2*t(df,-abs(beta[1,`i']/sqrt(var[1,`i'])))
}
scalar sum_sur = 0
forvalues i = 1(1)12 {
if (pval[`i',1] <= 0.05) {
*matrix ex = nullmat(ex) \ `i'
scalar sum_sur = sum_sur + beta[1,14+`i']
}
}
matrix ex_pere = nullmat(ex_pere) \ sum_sur
mat drop beta var pval
scalar drop df sum_sur
***storing year-total surprise rain: end
}
*sum checked for the 5 crops in excel: ok


esttab temp_* using "water_3rd_[put_date_here]_output_temp.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace

esttab pere_* using "water_3rd_[put_date_here]_output_pere.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace

*OBS: the produc value models were excluded from here in 01/11

***milk, 27/10
*(01/11) check if double miss leaves no attrition
count if d_miss == 0 & d_miss_milk == 0
dis 8712/12
*726, ok

sum $xlist_month
**(a) milk, no herd
capture mat drop beta var
capture mat drop ex_milk
capture scalar drop df sum_sur

foreach y in  "q_milk" {         
xtreg `y' $xlist_month i.time if d_miss == 0 & d_miss_milk == 0, fe vce(cluster geocodig_m)
estimates store milk_`y'
test

***storing year-total surprise rain: start
mat beta = e(b)
mat var = vecdiag(e(V))
scalar df = e(N) - 68
*71 vars (with intercept) - last ts*uf -time08 - time19
*CHECK DF
capture mat drop pval
forvalues i = 15(1)26 {
matrix pval = nullmat(pval) \ 2*t(df,-abs(beta[1,`i']/sqrt(var[1,`i'])))
}
scalar sum_sur = 0
forvalues i = 1(1)12 {
if (pval[`i',1] <= 0.05) {
*matrix ex = nullmat(ex) \ `i'
scalar sum_sur = sum_sur + beta[1,14+`i']
}
}
matrix ex_milk = nullmat(ex_milk) \ sum_sur
mat drop beta var pval
scalar drop df sum_sur
***storing year-total surprise rain: end
}
matlist ex_milk


capture mat drop beta var
capture mat drop ex_milk
capture scalar drop df sum_sur
**(b) milk, herd
foreach y in  "q_milk" {         
xtreg `y' $xlist_month herd i.time if d_miss ==0 & d_miss_milk == 0, fe vce(cluster geocodig_m)
estimates store milk_herd_`y'
test

***storing year-total surprise rain: start
mat beta = e(b)
mat var = vecdiag(e(V))
scalar df = e(N) - 68
*71 vars (with intercept) - last ts*uf -time08 - time19
*CHECK DF
capture mat drop pval
forvalues i = 15(1)26 {
matrix pval = nullmat(pval) \ 2*t(df,-abs(beta[1,`i']/sqrt(var[1,`i'])))
}
scalar sum_sur = 0
forvalues i = 1(1)12 {
if (pval[`i',1] <= 0.05) {
*matrix ex = nullmat(ex) \ `i'
scalar sum_sur = sum_sur + beta[1,14+`i']
}
}
matrix ex_milk = nullmat(ex_milk) \ sum_sur
mat drop beta var pval
scalar drop df sum_sur
***storing year-total surprise rain: end
}
matlist ex_milk


*esttab milk_* using "water_3rd_27_10_22_output_milk.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace

esttab milk_* using "water_3rd_[put_date_here]_output_milk.txt", se scalars(N chi2 F ll ll_0 p r2_a r2_o r2_w r2_b N_clust) nolines nogaps star(+ 0.10 * 0.05 ** 0.01 *** 0.001) varwidth(41) replace

*herd and milk summarys 02/11
global cropmilk q_cott q_cass q_maiz q_soyb q_bana q_coco q_coff q_palm q_pepp q_milk herd 

summ $cropmilk if d_miss == 0 & d_miss_temp == 0 & d_miss_pere == 0 & d_miss_milk == 0

**05/11 likelihood of 1SD in suprise rainfall
codebook sur_ppt_3 sur_ppt_4 if d_miss == 0
count if sur_ppt_3 >= 1 & sur_ppt_4 >= 1 
count if sur_ppt_3 <= -1 & sur_ppt_4 <= -1 




