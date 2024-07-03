********************************************************************************
* Program: master_paper_final
* Author: Emma Riley
* Created: JUL 2018
* Modified: 		
* Purpose: Replication of "Resisting Social Pressure in the Household Using Mobile Money: Experimental Evidence on Microenterprise Investment in Uganda" by Emma Riley, 2023
********************************************************************************

***************************
****** MODIFICATIONS ******
* Modified on June 28, 2024 to get the code to run through in one go
* Two main changes required:
* (1) Setting up an environment, in order to install the Stata packages without conflicts.
/* (2) The "multe" package has been updated, and no longer runs for this code due to a change in how it identifies the treatment. Accordingly, you need to run the code using the old version of the package.
	a. Go to the GitHub page for multe, navigate to the commit 27fe818b881d4bae5b844e1e73f9eda4f4288d71 (it should be from about three months ago), and download the full package to your computer.
	b. Point the corresponding "net install…" command to where you have stored this version of the package locally.
***************************

clear all
set scheme    s2color  //s1mono //sets colour scheme
set matsize 1000
set seed 99999
set sortseed 99999
version 18

* 1. Set up your user specific root directory to Replication Materials folder
global analysis "/Users/victoriamooers/Library/CloudStorage/GoogleDrive-vrm2120@columbia.edu/My Drive/replicationgame2024/194886-V1/replication"

*2. controls whether some more time consuming parts of the code are run. Set all to 1 to run everything 
global graph=1 // set to 1 to generate graphs 
global perm=1 // set to 1 to run the permutation analysis
global forest=1 // set to 1 to run the causal forest analysis
global install 1 // set to 1 if first time running 

global clean 0 // set to 1 to clean and construct the datasets
*cleaning globals
global bottom_w=0
global bottom_tr=0
global top_w=99 
global top_tr=99
global winsor=1
global trim=0
global shillings_conversion=3600 //usd to ugx

***
/* install any packages locally */
di "=== Redirecting where Stata searches for ado files ==="
capture mkdir "$analysis/ado_add"
sysdir
adopath 
*adopath - PERSONAL
*adopath - OLDPLACE
*adopath - SITE
sysdir set PLUS     "$analysis/ado_add/plus"
sysdir set PERSONAL "$analysis/ado_add"       // may be needed for some packages
sysdir

local author_adopath "ado"
if "`author_adopath'" != "" {             // The author adopath variable is filled out
    adopath ++ "$analysis/`author_adopath'"
}
***


if $install==1 {	
ssc install ciplot
ssc install coefplot
ssc install estout
ssc install texsave, replace
ssc install parmest
ssc install dsconcat
ssc install outreg2
ssc install qqvalue
ssc install icw_index
*move the ado file _gweightave from the ado folder to your personal ado directory

local github "https://raw.githubusercontent.com"
cap noi net uninstall multe
net install multe, from("$analysis/stata-multe-27fe818b881d4bae5b844e1e73f9eda4f4288d71") // you must go to the package's GitHub page, navigate to this commit, download it to your own computer, and call it from there -- because the old version is not available to download from github 

*1) ensure R is installed on your device - download from https://cran.r-project.org/ 
*2) download github to your device from https://desktop.github.com/
*3) install rcall
net install github, from("https://haghish.github.io/github/")
github install haghish/github
github install haghish/rcall, stable
*4) install and setup mlrtime
net install mlrtime, from("https://raw.githubusercontent.com/NickCH-K/mlrtime/master/")
mlrtimesetup, go

}


// outcomes
*primary
global main_results     earn_business    much_saved    capital   

*secondary
global profit t_sales sales monthly_profit weekly_profit 
global savings  saving net_saving use_saving_5 saving_amount_5 saving_goal_6 
global assets asset_ent_index ent_asset_value ent_asset_dummy ent_asset_count inventory_value hh_asset_value
global labour total_hoursbusiness   hours_week adult_hoursbusiness child_hoursbusiness non_hhemployee employee_hours 
global empower_all  switch_m  own_decision equal_decision  control_money remittance_share womans_income_share empower_1 empower_2
global happy happiness life_satisfaction  worry_money worry_money_dum
global household hh_income consumption_total consump_food  consump_nonfood_exsch consump_school 
global earnings other_work_earning spouse_earning  otherhh_earning spouse_bus_earn otherhh_bus_earn  spouse_allearning otherhh_allearning
global give_spouse give_spouse gave_spouse receive_spouse 
global records records_1 records_2 records_3 records_4 
global loan_use loan_use_business loan_use_family  loan_use_school loan_use_home  loan_use_exp loan_use_sav  loan_use_loan
global group group_talk group_receivehelp group_givehelp
global remittance remittance_samount remittance_ramount net_remittancerec  remittance_mm receive_remittance sent_remittance

*hetero
global hetero high_profits_base above_m_median_base abovem_inventory_base current_loan_base hyperbolic_base impatient_base abovem_risk_base abovem_sav_base  above_med_basset_base  married_base above_med_emp_base sent_fam_dummy_base spouse_fam_takes_base saving_goal_6_base spouse_bus_base hh_bus_base mm_close

global hetero_index hetero_perf_median hetero_selfc_median hetero_family_median

global index_vars empower 

*ids 
global baseline_ids  id_key_new treatment  loan_amount_dis today2_dis  disburse_noncompliance account_noncompliance strata_fixed_base branch_name_base 


*strata
global strata3 strata_fixed_base

*days controls
global days1 days days2

***************************************************************************
****************Data cleaning and constructione*****************************
**************************************************************************
if $clean==1 {
do "${analysis}/do-files/baseline_cleaning"
do "${analysis}/do-files/endline_cleaning"
do "${analysis}/do-files/combine_datasets"
do "${analysis}/do-files/brac_admin_cleaning"
do 	"${analysis}/do-files/clean_mno"
do "${analysis}/do-files/BRAC admin MNO merge"

}


*load data
use "${analysis}/input/survey_data.dta", clear

***************************************************************************
****************Tables using baseline sample*****************************
**************************************************************************

*take-up - Figure 1 and Table A3
preserve

do "${analysis}/do-files/analysis_takeup.do"

restore


//attrition - Table A4

tab treatment if consent==1
tab treatment if consent!=1

gen attrition=0
replace attrition=1 if consent!=1

eststo clear
eststo: areg attrition treatment2 treatment3, absorb($strata3) vce(robust)
sum attrition if treatment1==1
local control: di %9.3fc r(mean)
estadd scalar control `control'
test treatment2=treatment3
return list        
local test: di %9.3fc r(p)
estadd scalar test `test'

esttab  using "${analysis}/output/appendix/table_attrition.tex" , ///
					cells(b(fmt(%12.3fc)) se(fmt(%12.3fc) par)) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label ///
					stats(N r2 control  test, fmt(%9.0fc 3 3  ) label ("Observations" "R-squared"  "Control mean" "p-value MA=MD")) /// 
					title(Attrition \label{attrition}) postfoot(\bottomrule \multicolumn{2}{p{10cm}}{Linear regression of treatment indicators on a variable equal to one if the woman was not surveyed at endline. Regression controls for strata fixed effects. Robust standard errors in parentheses. }  \\ \end{tabular} \end{table}) replace


					
*Attrition on baseline correlates - Table A5
preserve 
local income loan_amount_dis weekly_profit   earn_business   much_saved hh_income
foreach x of local income {
replace `x'_base=`x'_base/1000
}

ren treatment2 treatment2_base
ren treatment3 treatment3_base

eststo clear

reg attrition  respondent_age_base married_base hh_size_base completed_primary_base  completed_secondary_base work_occupation_base ///
  loan_amount_dis_base weekly_profit_base  high_profits_base current_loan_base  much_saved_base mobile_account_base ///
  hyperbolic_base impatient_base womans_income_share_base   spouse_fam_takes_base 
testparm respondent_age_base married_base hh_size_base completed_primary_base  completed_secondary_base work_occupation_base ///
  loan_amount_dis_base weekly_profit_base  high_profits_base current_loan_base  much_saved_base mobile_account_base ///
  hyperbolic_base impatient_base womans_income_share_base   spouse_fam_takes_base 
local F1 = `r(p)'

local variables  "treatment2 treatment3 respondent_age married hh_size completed_primary  completed_secondary work_occupation  loan_amount_dis weekly_profit  high_profits current_loan  much_saved mobile_account hyperbolic impatient womans_income_share above_m_median  spouse_fam_takes " 

local count=0
foreach x of local variables {
loc count = `count' + 1	
eststo attrition: reg  attrition `x'_base, vce(robust) 
estadd 	scalar Ftest `F1'

if `count'==1  {
esttab  using "${analysis}/output/appendix/table_attrition2.tex", ///
					cells(b(fmt(%9.3f)) se(fmt(%9.3f)  par)) ///
					nostar ///
					keep (`x'_base)   noobs nonotes  noline nocons collabels(none)  replace f eqlabels(none) label ///
					booktabs    ///
					prehead(\begin{table}[hbt] \caption{Correlates of attrition } \label{attrition2}  \centering \resizebox{0.47\textheight}{!}{ \begin {tabular}{lc} \hline ) posthead(\hline)
eststo clear 
}
if `count'>1 & `count'!=19 {
esttab  using "${analysis}/output/appendix/table_attrition2.tex" , ///
					append f cells(b(fmt(%9.3f)) se(fmt(%9.3f)  par)) ///
					nostar ///
					keep (`x'_base)  nomtitle noline nonumber  collabels(none) plain nocons noobs eqlabels(none) label ///
						booktabs nodep   
eststo clear
}

if `count'==19 {
esttab  using "${analysis}/output/appendix/table_attrition2.tex" , ///
					append f cells(b(fmt(%9.3f)) se(fmt(%9.3f)  par)) ///
					nostar ///
					keep (`x'_base) nomtitles noline nonumber  collabels(none)  nocons noobs eqlabels(none) booktabs label ///
					stats(Ftest N, fmt(2 0) label ( "\hline F-test p-value" "Observations")) ///
					postfoot(\bottomrule \multicolumn{2}{p{11cm}}{ Linear regression of baseline characteristics on a variable equal to one if the woman was not surveyed at endline. Each row represents a separate regression. Monetary amounts in '000 USD and winsorized at the 99\% level. The F-test p-value comes from regressing the attrition variable on all the characteristics and testing if they are jointly zero. Robust standard errors in parentheses. } \\ \end{tabular} } \end{table})
						 
eststo clear
}
}

restore



//balance test - Table A2
preserve
do "${analysis}/do-files/analysis_balancecheck2.do"
restore



//final sample for subsequent analysis
drop if consent!=1


********************************************************************
***************************Main tables******************************
********************************************************************

//1. primary outcomes - Table 1 
*q-values

eststo clear
preserve 
local i=1
foreach x of global main_results {
tempfile tf`i' 
eststo: areg  `x' `x'_base treatment2 treatment3 , absorb($strata3) vce(robust)
parmest, format(p %12.9f) saving(`tf`i'',replace) idn(`i')
sum `x' if treatment1==1
local control: di %9.2fc r(mean)
estadd scalar control `control'
sum `x'_base if treatment1==1
local control_base: di %9.2fc r(mean)
estadd scalar control_base `control_base'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd scalar test `test'
local i=`i'+1

}


foreach x of numlist 1/3 {	
	matrix define Tq`x' = (. , .)
	matrix colnames Tq`x' = treatment2 treatment3
}


dsconcat `tf1' `tf2' `tf3'
keep if parm=="treatment2"
keep p idnum
qqvalue p, method(yekutieli) qvalue(q)
format q   %12.8fc
*qqvalue p, method(bonferroni) qvalue(q2)
foreach x of numlist 1/`i' {
global T2q`x'=q[`x']
display ${T2q`x'}
}

foreach x of numlist 1/3 {
matrix Tq`x'[1,1] = (q[`x'])
}


dsconcat `tf1' `tf2' `tf3'
keep if parm=="treatment3"
keep p idnum
qqvalue p, method(yekutieli) qvalue(q)
format q %9.6fc
*qqvalue p, method(bonferroni) qvalue(q2)
foreach x of numlist 1/`i' {
global T1q`x'=q[`x']
display ${T1q`x'}
}

foreach x of numlist 1/3 {
matrix Tq`x'[1,2] = (q[`x'])
}

foreach x of numlist 1/3{
estadd matrix q = Tq`x' : est`x'
}

restore


esttab  using "${analysis}/output/table_primary_results.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par) ///
				p( par([ ])) ///
				q( par(\{ \}) pvalue(q))) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label ///
					stats(N r2 control control_base test, fmt(%9.0fc 2 2 2 ) label ("Observations" "R-squared"  "Control mean" "Control mean baseline" "p-value MA=MD")) /// 
					title(Treatment effects on woman's business profits, savings and business capital \label{results}) postfoot(\bottomrule \multicolumn{4}{p{12cm}}{ Intent-to-treat estimates. All outcomes are winsorized at the 99\% level. USD. All regressions include strata dummies and include the baseline value of the outcome. Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account.  Profits refers to the woman's self-reported monthly business profit. Savings is individual savings held by the woman. Capital is the value of all assets the woman uses in her business plus the value of inventory held for her business. Control mean endline is the mean value of the outcome in the control group at endline. Control mean baseline is the mean value of the outcome in the control group at baseline.   False discovery rate (FDR) adjusted p-values, also known as q-values, were used to correct for multiple hypothesis testing. They are shown in curly brackets. These were calculated following the method of \cite{Benjamini2006}. Robust p-values in square brackets. Robust standard errors in parentheses. }  \\ \end{tabular} \end{table}) replace




//Table 2 heterogeneity by sc and fp
*hetero by index
	foreach y of global hetero_index  {
	gen treatment2_`y'=treatment2*`y'
	gen treatment3_`y'=treatment3*`y'
}
la var  treatment2_hetero_selfc_median "MA*self control"
la var treatment3_hetero_selfc_median "MD*self control"
la var treatment2_hetero_family_median "MA*family pressure"
la var treatment3_hetero_family_median "MD*family pressure"

eststo clear
preserve 
local i=1
foreach x of global main_results {
tempfile tf`i' 
eststo: 	areg `x'  `x'_base treatment2 treatment3 treatment2_hetero_selfc_median treatment3_hetero_selfc_median  treatment2_hetero_family_median treatment3_hetero_family_median hetero_selfc_median hetero_family_median , absorb($strata3) vce(robust)
parmest, format(p %12.9f) saving(`tf`i'',replace) idn(`i')
sum `x' if treatment==0 & hetero_selfc_median>0
local c1: di %9.2fc r(mean)
estadd scalar control_sc `c1'
sum `x' if treatment==0 & hetero_family_median>0
local c2: di %9.2fc r(mean)
estadd scalar control_fp `c2'
sum `x'_base if treatment==0 & hetero_selfc_median>0
local c3: di %9.2fc r(mean)
estadd scalar controlb_sc `c3'
sum `x'_base if treatment==0 & hetero_family_median>0
local c4: di %9.2fc r(mean)
estadd scalar controlb_fp `c4'

test treatment3_hetero_family_median=treatment3_hetero_selfc_median
local test_temp: di %9.2fc r(p)
estadd scalar test_md `test_temp'
test treatment2_hetero_family_median=treatment2_hetero_selfc_median
local test_temp2: di %9.2fc r(p)
estadd scalar test_ma `test_temp2'

local i=`i'+1

}

foreach x of numlist 1/3 {	
	matrix define Tq`x' = (. , . , . , ., . , .)
	matrix colnames Tq`x' = treatment2 treatment3 treatment2_hetero_selfc_median treatment3_hetero_selfc_median treatment2_hetero_family_median  treatment3_hetero_family_median 
}

dsconcat `tf1' `tf2' `tf3'
keep if parm=="treatment2"
keep p idnum
qqvalue p, method(yekutieli) qvalue(q) 
replace q=0.99 if q==1
format q %12.0g
*qqvalue p, method(bonferroni) qvalue(q2)
foreach x of numlist 1/3 {
matrix Tq`x'[1,1] = (q[`x'])
}


dsconcat `tf1' `tf2' `tf3'
keep if parm=="treatment3"
keep p idnum
qqvalue p, method(yekutieli) qvalue(q)
replace q=0.99 if q==1
format q %12.0g
*qqvalue p, method(bonferroni) qvalue(q2)
foreach x of numlist 1/3 {
matrix Tq`x'[1,2] = (q[`x'])
}

dsconcat `tf1' `tf2' `tf3'
keep if parm=="treatment2_hetero_selfc_median"
keep p idnum
qqvalue p, method(yekutieli) qvalue(q)
replace q=0.99 if q==1
format q %12.0g
*qqvalue p, method(bonferroni) qvalue(q2)
foreach x of numlist 1/3 {
matrix Tq`x'[1,3] = (q[`x'])
}

dsconcat `tf1' `tf2' `tf3'
keep if parm=="treatment3_hetero_selfc_median"
keep p idnum
qqvalue p, method(yekutieli) qvalue(q)
replace q=0.99 if q==1
format q %12.0g
*qqvalue p, method(bonferroni) qvalue(q2)
foreach x of numlist 1/3 {
matrix Tq`x'[1,4] = (q[`x'])
}



dsconcat `tf1' `tf2' `tf3'
keep if parm=="treatment2_hetero_family_median"
keep p idnum
qqvalue p, method(yekutieli) qvalue(q)
replace q=0.99 if q==1
format q %12.0g
*qqvalue p, method(bonferroni) qvalue(q2)
foreach x of numlist 1/3 {
matrix Tq`x'[1,5] = (q[`x'])
}


dsconcat `tf1' `tf2' `tf3'
keep if parm=="treatment3_hetero_family_median"
keep p idnum
qqvalue p, method(yekutieli) qvalue(q)
replace q=0.99 if q==1
format q %12.0g
foreach x of numlist 1/3 {
matrix Tq`x'[1,6] = (q[`x'])
}

foreach x of numlist 1/3{
estadd matrix q = Tq`x' : est`x'
}

restore


esttab  using "${analysis}/output/table_hetero_index_family.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par) ///
				q(par(\{ \}) pvalue(q) keep(treatment2_hetero_selfc_median treatment3_hetero_selfc_median  treatment2_hetero_family_median treatment3_hetero_family_median treatment2 treatment3))) ///
					nostar ///
					keep (treatment2_hetero_selfc_median treatment3_hetero_selfc_median  treatment2_hetero_family_median treatment3_hetero_family_median treatment2 treatment3 hetero_family_median)     collabels(none)  nocons  eqlabels(none) booktabs label ///
					stats(N r2 control_sc control_fp controlb_sc controlb_fp test_md test_ma, fmt(%9.0fc 2 2 2 ) label ("Observations" "R-squared"  "Control mean self control" "Control mean family pressure" "Control mean baseline self-control" "Control mean baseline family pressure" "p-val MD self control=MD family pressure" "p-val MA self control=MA family pressure")) /// 
					title(Heterogeneous treatment effects by baseline self control and family pressure index \label{hetero_index_family}) postfoot(\bottomrule \multicolumn{4}{p{14cm}}{ Intent-to-treat estimates. Monetary outcomes are winsorized at the 99\% level and in USD. All regressions include strata dummies. Mobile Account (MA) is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse (MD) is the treatment where a mobile money account was provided and the loan also disbursed onto this account. Heterogeneous indexes are defined in section \ref{mechanisms} and the construction is shown in Appendix Table \ref{index}. The interaction is for someone who is above the median in the index. Profit is self-reported monthly profit.  Savings is total savings in each form of saving used. Capital is composed of business assets and inventories.    False discovery rate (FDR) adjusted p-values, also known as q-values, were used to correct for multiple hypothesis testing. They are shown in square brackets. These were calculated following the method of \cite{Benjamini2006}. Robust standard errors in parentheses. }  \\ \end{tabular} \end{table}) substitute(\_ _)  replace



																
//Tables 3 (loan use) 
eststo clear
foreach x of global loan_use {
		eststo: areg  `x'  treatment2 treatment3, absorb($strata3) vce(robust)
		sum `x' if treatment==0
		local control: di %9.4fc  r(mean)
		estadd scalar control `control'
		test treatment2=treatment3
		local test: di %9.4fc r(p)
		estadd scalar test `test'

}

esttab  using "${analysis}/output/table_second_loan_use.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control test , fmt(%9.0fc 2 2 2) label ("Observations" "R-squared"  "Control mean"  "p-value MA=MD" )) /// 
					 prehead(\begin{table}[htb] \centering \caption{Treatment effects on use of the loan in the first week after disbursement} \label{second_loan}  \begin{tabulary}{1\textwidth}{lCcccCCc} \hline & (1) & (2) & (3) & (4) & (5) & (6) & (7) \\ & Business & Sharing & School & Home & Expenditure & Saving & Loan \\) postfoot(\bottomrule \multicolumn{8}{p{15.5cm}}{Not specified in pre-analysis plan.  Intent-to-treat estimates. All outcomes are winsorized at the 99\% level. USD. All regressions include strata dummies. Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account.  Amount of loan spent on each category 1 week after receiving loan. Business is business inventory and assets, sharing is money given to the spouse, friends or other family members, both at home and elsewhere, school is money spent on school fees and related expenditures, home is money spent on items for the home or home improvements, expenditure is money spent on food, clothes, transport etc. and loan is money spent paying back other loans. Recall 8 months later. Average loan amount was USD 380. Robust standard errors in parentheses. } \\ \end{tabulary} \end{table}) replace substitute(\_ _) 

//Tables  4 (give spouse) 
eststo clear
foreach x of global give_spouse {
		eststo: areg  `x' `x'_base treatment2 treatment3, absorb($strata3) vce(robust)
		sum `x' if treatment==0
		local control: di %9.4fc  r(mean)
		estadd scalar control `control'
		sum `x'_base if treatment==0
		local control_base: di %9.4fc  r(mean)
		estadd scalar control_base `control_base'
		test treatment2=treatment3
		local test: di %9.4fc r(p)
		estadd scalar test `test'

}

esttab  using "${analysis}/output/table_second_give_spouse.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control control_base test , fmt(%9.0fc 2 2 2 2) label ("Observations" "R-squared"  "Control mean" "Control mean baseline" "p-value MA=MD" )) /// 
					 prehead(\begin{table}[hbt] \centering \caption{Treatment effects on amount and whether the woman gave money to her spouse and amount received from her spouse in the last month} \label{give_spouse}  \begin{tabulary}{1\textwidth}{lCCC} \hline  & (1) & (2) & (3) \\  & Amount Given Spouse & Dummy Gave Spouse Money & Amount Received Spouse \\) postfoot(\bottomrule \multicolumn{4}{p{15cm}}{Not in pre-analysis plan. Intent-to-treat estimates. All outcomes are winsorized at the 99\% level. '000 Ugandan Shillings.  All regressions include strata dummies and include the baseline value of the outcome. Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account. Amount give/received spouse is the monthly transfer to/from the spouse. Robust standard errors in parentheses.  } \\ \end{tabulary} \end{table}) replace substitute(\_ _) 

					 
					 
//mm offered - Table 5 first 2 columns
eststo clear
eststo:areg mm_offered treatment2 treatment3, absorb($strata3) vce(robust)
sum  mm_offered if treatment1==1
local control: di %9.4fc r(mean)
estadd scalar control `control'
test treatment2=treatment3
return list        
local test: di %9.4fc r(p)
estadd scalar test `test'

eststo:areg mm_offered treatment2  treatment3 treatment2_hetero_family_median treatment3_hetero_family_median hetero_family_median, absorb($strata3) vce(robust)
sum  mm_offered if treatment1==1 & hetero_family_median==1
local control: di %9.4fc r(mean)
estadd scalar control `control'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd scalar test `test'


*subsequent loan - columns 3 onwards Table 5
preserve
use "${analysis}/input/BRAC_MM_merged", clear

	foreach y of varlist hetero_family_median  {
	gen treatment2_`y'=treatment2*`y'
	gen treatment3_`y'=treatment3*`y'
}
la var treatment2_hetero_family_median "MA*family pressure"
la var treatment3_hetero_family_median "MD*family pressure"

foreach x of varlist deposit_any deposit_max deposit_share {
eststo: areg `x'   treatment3 , absorb($strata3) vce(robust)
sum `x' if treatment2==1
local control: di %9.4fc r(mean)
estadd scalar control `control'
eststo:areg `x'   treatment3  treatment3_hetero_family_median hetero_family_median, absorb($strata3) vce(robust)
sum  `x' if treatment2==1 & hetero_family_median==1
local control: di %9.4fc r(mean)
estadd scalar control `control'
}


restore


esttab  using "${analysis}/output/table_subsequent.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment2 treatment3 treatment2_hetero_family_median treatment3_hetero_family_median )     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control  test , fmt(%9.0fc 2  2 2) label ("Observations" "R-squared"  "Control mean" "p-value MA=MD" )) /// 
					 prehead(\begin{table}[hbt]\centering\caption{Treatment effects on stated and actual use of the mobile money account for subsequent loans} \label{tab:subsequent_loan}\resizebox{1\textwidth}{!}{   \begin{tabulary}{1.4\textwidth}{lCCCCCCCC} \hline&  \multicolumn{2}{c}{Stated preferences}  & \multicolumn{6}{c}{Own deposit behaviour} \\ \cmidrule(lr){2-3}    \cmidrule(lr){4-9} & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) \\ & Want Mobile Deposit & Want Mobile Deposit  & Any Own Deposit & Any Own Deposit  & Own Deposit Amount & Own Deposit Amount & Own Deposit Share &Own Deposit Share \\) postfoot(\bottomrule \multicolumn{9}{p{22cm}}{Intent-to-treat estimates.  All regressions include strata dummies. Stated preferences refer to the responses women gave when asked hypothetical questions about receiving future loans on a mobile money account. Want mobile deposit means that the woman indicated in the endline survey that she would like to receive a subsequent loan on a mobile money account. Own deposit variables capture if the woman made any despoit of her subsequent loan herself to the mobile money account. A subsequent loan is any loan disbursed in 2017 after the loan that disbursement was randomised for in this study. This is only examined within the treatment groups. Any own deposit any means a deposit was made to the mobile money account in the 2 week period after the subsequent loan was disbursed. Own deposit amount is the maximum single deposited amount in the 2 weeks after the subsequent loan was disbursed, in USD. Own deposit share  is the value of the deposited amount as a share of the subsequent loan amount. Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account. Family pressure is a dummy variable if the woman experienced sharing pressure from her family that was above the median at baseline.  Control mean is the mean value of that outcome in the control group for stated preferences and the Mobile Account treatment group for own deposit behaviour. Robust standard error in parentheses. } \\ \end{tabulary}} \end{table}) replace substitute(\_ _) 


//Tables  6 (household)  
eststo clear
foreach x of global household {
		eststo: areg  `x' `x'_base  treatment2 treatment3, absorb($strata3) vce(robust)
		sum `x' if treatment==0
		local control: di %9.4fc  r(mean)
		estadd scalar control `control'
		sum `x'_base if treatment==0
		local control_base: di %9.4fc  r(mean)
		estadd scalar control_base `control_base'
		test treatment2=treatment3
		local test: di %9.4fc r(p)
		estadd scalar test `test'

}

esttab  using "${analysis}/output/table_second_household.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control control_base test , fmt(%9.0fc 2 2 2 2) label ("Observations" "R-squared"  "Control mean" "Control mean baseline" "p-value MA=MD" )) /// 
					 prehead(\begin{table}[hbt] \centering \caption{Treatment effects on household outcomes} \label{second_consum}  \begin{tabulary}{1\textwidth}{lCCCCC} \hline & (1) & (2) & (3) & (4) & (5)  \\ & Total Household Income & Total Household Consumption & Food Consumption & Non-food Consumption exl school & School Expenditure \\) postfoot(\bottomrule \multicolumn{6}{p{16cm}}{Intent-to-treat estimates. All outcomes are winsorized at the 99\% level. USD. All regressions include strata dummies and include the baseline value of the outcome.  Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account.  All values are monthly for the entire household. Non-food consumption excludes temptation spending and transfers. Control mean endline is the mean value of the outcome in the control group at endline. Control mean baseline is the mean value of the outcome in the control group at baseline.  Robust standard errors in parentheses. } \\ \end{tabulary} \end{table}) replace substitute(\_ _) 
					 
										 
********************************************************					 
***********************APPENDIX**************************
********************************************************
*loan size distribution - Figure A1
if ${graph}==1 {

histogram loan_amount_dis_base, color(blue%60)  ylabel(,grid) xtitle("Loan Size (USD)") ///
title("Loan Size Distribution") graphregion(color(white) lc(white) lw(med)) bgcolor(white) xlabel(200 400 600 800  1000 1200)
graph export "$analysis/output/appendix/graphs/hist_disbursed.png", replace
}


//business types - figure A2 
if ${graph}==1 {
gen one = 1
graph hbar (sum) one, over(business_type_base) ytitle(frequency) bar(1, color(blue%60)) plotregion(ilstyle(none) lcolor(white))  ylabel(,grid) graphregion(color(white) lc(white) lw(med)) bgcolor(white)
 graph export "${analysis}/output/appendix/graphs/business_types.png", replace
}


//figure A3 and A4
if $graph==1 {
foreach x of global main_results { 
capture cumul 	`x' if treatment==0, 	gen(cumul1)
capture cumul 	`x' if treatment==1, 	gen(cumul2)
capture cumul 	`x' if treatment==2	, 	gen(cumul3)

global 	StarSize 	= 1.5
twoway 	(line cumul1 `x'	, lpattern(solid) lcolor(green%60) sort) ///
(line cumul2 `x'	, lpattern(shortdash) lcolor(red%60) sort) ///
(line cumul3 `x'	, lpattern(longdash) lcolor(blue%60) sort), /// 
 ytitle("Percentile") ylabel(0.05 0.10 0.25 0.50 0.75 0.90 0.95)  legend(label(1 "Control") label(2 "Account") label(3 "Disburse")  rows(1)) graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white) icolor(white) ilstyle(none) lcolor(white))  ylabel(,grid) 
graph export "${analysis}//output/appendix/graphs/CDF_`x'.png", replace
drop cumul*

capture drop ln_`x'
gen ln_`x'=ln(`x')
cumul 	ln_`x' if treatment==0, 	gen(cumul1)
cumul 	ln_`x' if treatment==1, 	gen(cumul2)
cumul 	ln_`x' if treatment==2	, 	gen(cumul3)

global 	StarSize 	= 1.5
twoway 	(line cumul1 ln_`x'	, lpattern(solid) lcolor(green%60) sort) ///
(line cumul2 ln_`x'	, lpattern(shortdash) lcolor(red%60)  sort) ///
(line cumul3 ln_`x'	, lpattern(longdash) lcolor(blue%60) sort), /// 
 ytitle("Percentile") ylabel(0.05 0.10 0.25 0.50 0.75 0.90 0.95)  legend(label(1 "Control") label(2 "Account") label(3 "Disburse")  rows(1)) graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white) ilstyle(none) icolor(white) lcolor(white))  ylabel(,grid)
graph export "${analysis}//output/appendix/graphs/CDF_ln_`x'.png", replace
drop cumul*
}

gen treat23=0 if treatment2==1 | treatment3==1
replace treat23=1 if treatment3==1
ksmirnov(earn_business), by(treat23)
ksmirnov(capital), by(treat23)
ksmirnov(much_saved), by(treat23)
gen treat13=0 if treatment1==1 | treatment3==1
replace treat13=1 if treatment3==1
ksmirnov(earn_business), by(treat13)
ksmirnov(capital), by(treat13)
ksmirnov(much_saved), by(treat13)
gen treat12=0 if treatment2==1 | treatment1==1
replace treat12=1 if treatment2==1
ksmirnov(earn_business), by(treat12)
ksmirnov(capital), by(treat12)
ksmirnov(much_saved), by(treat12)
}

//Table A1
tab switch_m //note numbers directly copied into table, no output generated. 

// Tables A6 (profit) 
eststo clear
foreach x of global profit {
		eststo: areg  `x' `x'_base  treatment2 treatment3, absorb($strata3) vce(robust)
		sum `x' if treatment==0
		local control: di %9.4fc  r(mean)
		estadd scalar control `control'
		sum `x'_base if treatment==0
		local control_base: di %9.4fc  r(mean)
		estadd scalar control_base `control_base'
		test treatment2=treatment3
		local test: di %9.4fc r(p)
		estadd scalar test `test'

}

esttab  using "${analysis}/output/appendix/table_second_profit.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control control_base test , fmt(%9.0fc 2 2 2 2) label ("Observations" "R-squared"  "Control mean" "Control mean baseline" "p-value MA=MD" )) /// 
					 prehead(\begin{table}[htb] \centering \caption{Treatment effects on additional business outcomes} \label{second_profit} \begin{tabulary}{1\textwidth}{lCCCC} \hline  & (1) & (2) & (3) & (4)  \\ & Monthly Sales & Weekly Sales & Monthly Profit & Weekly Profit \\) postfoot(\bottomrule \multicolumn{5}{p{15cm}}{Intent-to-treat estimates. All outcomes are winsorized at the 99\% level. USD.  All regressions include strata dummies and include the baseline value of the outcome. Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account.  Monthly and weekly profit are calculated by subtracting the corresponding expenditures from  sales. All outcomes refer to the woman's business. Robust standard errors in parentheses. } \\ \end{tabulary} \end{table}) replace substitute(\_ _) 

//A7 (savings) 					 
eststo clear
foreach x of global savings {
	capture confirm variable `x'_base 
	if !_rc {
		eststo: areg  `x' `x'_base  treatment2 treatment3, absorb($strata3) vce(robust)
		sum `x' if treatment==0
		local control: di %9.4fc  r(mean)
		estadd scalar control `control'
		sum `x'_base if treatment==0
		local control_base: di %9.4fc  r(mean)
		estadd scalar control_base `control_base'
		test treatment2=treatment3
		local test: di %9.4fc r(p)
		estadd scalar test `test'
	}
	else {
		eststo: areg  `x'  treatment2 treatment3, absorb($strata3) vce(robust)
		sum `x' if treatment==0
		local control: di %9.4fc  r(mean)
		estadd scalar control `control'
		test treatment2=treatment3
		local test: di %9.4fc r(p)
		estadd scalar test `test'		
		
		
	}

}

esttab  using "${analysis}/output/appendix/table_second_savings.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control control_base test , fmt(%9.0fc 2 2 2 2) label ("Observations" "R-squared"  "Control mean" "Control mean baseline" "p-value MA=MD" )) /// 
					 prehead( \begin{table}[h!tb] \centering \caption{Treatment effects on saving outcomes} \label{second_saving}  \begin{tabulary}{1\textwidth}{lCCCCC} \hline & (1) & (2) & (3) & (4) & (5) \\ & Calculated Saving & Net Saving & Saves Mobile Money & Amount Mobile Money & Saving Goal Business \\) postfoot(\bottomrule \multicolumn{6}{p{15cm}}{Intent-to-treat estimates.  All outcomes are winsorized at the 99\% level. USD. All regressions include strata dummies. Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account.  All outcomes reported here were only collected at endline. Calculated savings is the sum of savings in each form of saving. Net savings is additions-withdrawals from savings in the last month. Saves mobile money is a dummy equal to one if the the respondent reported saving on a mobile money account. Amount mobile money is the value of savings on a mobile money account. Saving goal business is a dummy if the reported goal of saving is to use it for the business.  Robust standard errors in parentheses.} \\ \end{tabulary} \end{table}) replace substitute(\_ _) 

					 
// A8 (assets) 
eststo clear
foreach x of global assets {
		eststo: areg  `x' `x'_base  treatment2 treatment3, absorb($strata3) vce(robust)
		sum `x' if treatment==0
		local control: di %9.4fc  r(mean)
		estadd scalar control `control'
		sum `x'_base if treatment==0
		local control_base: di %9.4fc  r(mean)
		estadd scalar control_base `control_base'
		test treatment2=treatment3
		local test: di %9.4fc r(p)
		estadd scalar test `test'

}

esttab  using "${analysis}/output/appendix/table_second_assets.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control control_base test , fmt(%9.0fc 2 2 2 2) label ("Observations" "R-squared"  "Control mean" "Control mean baseline" "p-value MA=MD" )) /// 
					 prehead(\begin{table}[hbt] \centering \caption{Treatment effects on business and household asset outcomes} \label{second_asset}  \resizebox{1\textwidth}{!}{\begin{tabulary}{1.4\textwidth}{lCCCCCC} \hline & \multicolumn{5}{c}{Business} & Household \\  \cmidrule(lr){2-6}   \cmidrule(lr){7-7}  & (1) & (2) & (3) & (4) & (5) & (6)\\ & PCA Index Bus. Assets & Value Bus. Assets & Unique Bus. Assets & Count Bus. Assets & Inventory Value & Household Wealth \\) postfoot(\bottomrule \multicolumn{7}{p{22cm}}{Intent-to-treat estimates. All outcomes are winsorized at the 99\% level. USD. All regressions include strata dummies and include the baseline value of the outcome.  Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account. Principal component analysis of assets used in the business. Higher values mean a larger number of different assets are used in the business. Household wealth includes the value of all assets used only in the household (and not the business). Control mean endline is the mean value of the outcome in the control group at endline. Control mean baseline is the mean value of the outcome in the control group at baseline. Robust standard errors in parentheses.} \\ \end{tabulary}} \end{table}) replace substitute(\_ _) 

					 
					 
// A9 (labour)
eststo clear
foreach x of global labour {
		eststo: areg  `x' `x'_base  treatment2 treatment3, absorb($strata3) vce(robust)
		sum `x' if treatment==0
		local control: di %9.4fc  r(mean)
		estadd scalar control `control'
		sum `x'_base if treatment==0
		local control_base: di %9.4fc  r(mean)
		estadd scalar control_base `control_base'
		test treatment2=treatment3
		local test: di %9.4fc r(p)
		estadd scalar test `test'

}

esttab  using "${analysis}/output/appendix/table_second_labour.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control control_base test , fmt(%9.0fc 2 2 2 2) label ("Observations" "R-squared"  "Control mean" "Control mean baseline" "p-value MA=MD" )) /// 
					 prehead(\begin{table}[hbt] \centering \caption{Treatment effects on household labour outcomes} \label{second_labour}  \begin{tabulary}{1\textwidth}{lCCCCCC} \hline & (1) & (2) & (3) & (4) & (5) & (6) \\ & All Hours & Woman's Hours & Adult Family Hours & Child Family Hours & No. Employees & Employee Hours \\) postfoot(\bottomrule \multicolumn{7}{p{15cm}}{Intent-to-treat estimates.  All outcomes are winsorized at the 99\% level. Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account. All regressions include strata dummies and include the baseline value of the outcome. All variables refer to hours or employment in the woman's business. All hours is composed of columns (2), (3), (4) and (6). Robust standard errors in parentheses.} \\ \end{tabulary} \end{table}) replace substitute(\_ _) 
					 
		

*do impacts differ by baseline business type - Table A10
eststo clear 
foreach x of varlist bus1-bus17  {
eststo: areg  `x' `x'_base treatment2 treatment3, absorb($strata3) vce(robust)
sum `x' if treatment1==1
local control: di %9.2fc r(mean)
estadd scalar control `control'
sum `x'_base if treatment1==1
local control_base: di %9.2fc r(mean)
estadd scalar control_base `control_base'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd scalar test `test'
}

eststo: areg  change_business treatment2 treatment3, absorb($strata3) vce(robust)
sum change_business if treatment1==1
local control: di %9.2fc r(mean)
estadd scalar control `control'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd scalar test `test'

esttab  using "${analysis}/output/appendix/table_second_bus_type.tex" , ///
					cells(b(fmt(%12.3fc)) se(fmt(%12.3fc) par)) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control control_base  test , fmt(%9.0fc 2  2 2) label ("Observations" "R-squared"  "Control mean" "Control mean baseline" "p-value MA=MD" )) /// 
					 prehead(\afterpage{ \begin{landscape} \begin{table}[H] \centering \caption{Treatment effects on business type} \label{tab:bus_type} \resizebox{1.5\textwidth}{!}{ \begin{tabulary}{2.4\textwidth}{LCCCCCCCCCCCCCCCCCC} \hline & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) & (13) & (14) & (15) & (16) & (17) & (18) \\ & Agriculture & Beauty \& Hairdressing & Boda Boda & Brick laying & Charcoal seller & Cook  & Food stall & Hawker & Landlord & Mobile money agent & Other  & Restaurant/bar & Shop & Seamstress & Laundry & Clothes resale & Drug store & Change business \\) postfoot(\bottomrule \multicolumn{19}{p{38cm}}{Intention-to-treat estimates. Each column shows a dummy variable for whether the woman reported that industry as a her primary business at endline. All regressions control for whether the woman was also doing that business at baseline, as well as strata dummies. Change business is a dummay variable capturing if the  business is different at endline than baseline. Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account.  Control mean endline is the mean value of the outcome in the control group at endline. Control mean baseline is the mean value of the outcome in the control group at baseline. Robust standard errors in parentheses.} \\ \end{tabulary}} \end{table} \end{landscape}}) replace substitute(\_ _) 



//winsorizing robustness Table A11
eststo clear
foreach x of global main_results {
eststo: areg  `x'_100 `x'_100_base treatment2 treatment3 , absorb($strata3) vce(robust)
sum `x'_100 if treatment1==1
local control: di %9.2fc r(mean)
display "`control'"
estadd local control `control'
sum `x'_100_base if treatment1==1
local control_base: di %9.2fc r(mean)
estadd local control_base `control_base'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd local test `test'
}
esttab 	using "$analysis/output/appendix/table_robust_winsor.tex", replace f ///
								prehead(\begin{table}[htbp]\centering \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi} \caption{Robustness Checks - winsorizing \label{winsor}} \begin{tabular}{l*{@M}{c}} \toprule) ///
								posthead(\hline \multicolumn{4}{l}{\textit{No winsorizing}} \\) ///
								b(2) se(2) nostar  /// 
								keep(treatment2 treatment3 ) ///
								varlabels(treatment2 "Mobile Account" treatment3 "Mobile Disburse") ///
								label  booktabs noobs nonotes  collabels(none) ///
								mtitle("Profit" "Saving" "Capital" ) ///
								stats( control control_base test, fmt(2 2 2) ///
								label ("Control Mean" "Control mean baseline" "p-value MA=MD" )) 

eststo clear

foreach x of global main_results {
eststo:areg  `x'_995 `x'_995_base treatment2 treatment3 , absorb($strata3) vce(robust)
sum `x'_995 if treatment1==1
local control: di %9.2fc r(mean)
display "`control'"
estadd local control `control'
sum `x'_995_base if treatment1==1
local control_base: di %9.2fc r(mean)
estadd local control_base `control_base'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd local test `test'
}

esttab 	using "$analysis/output/appendix/table_robust_winsor.tex", append f ///
								posthead(\hline \multicolumn{4}{l}{\textit{Winsorizing 99.5\%}} \\) ///
								b(2) se(2) nostar  /// 
								keep(treatment2 treatment3 ) ///
								varlabels(treatment2 "Mobile Account" treatment3 "Mobile Disburse") ///
								label  booktabs nonum nomtitles noobs nonotes  collabels(none) ///
								stats( control control_base test, fmt(2 2 2) ///
								label ("Control Mean" "Control mean baseline" "p-value MA=MD" )) 


eststo clear


foreach x of global main_results {
eststo:areg  `x'_98 `x'_98_base treatment2 treatment3 , absorb($strata3) vce(robust)
sum `x'_98 if treatment1==1
local control: di %9.2fc r(mean)
display "`control'"
estadd local control `control'
sum `x'_98_base if treatment1==1
local control_base: di %9.2fc r(mean)
estadd local control_base `control_base'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd local test `test'
}

esttab 	using "$analysis/output/appendix/table_robust_winsor.tex", append f ///
								posthead(\hline \multicolumn{4}{l}{\textit{Winsorizing 98\%}} \\) ///
								b(2) se(2) nostar  /// 
								keep(treatment2 treatment3 ) ///
								varlabels(treatment2 "Mobile Account" treatment3 "Mobile Disburse") ///
								label  booktabs nonum nomtitles noobs nonotes  collabels(none) ///
								stats( control control_base test, fmt(2 2 2) ///
								label ("Control Mean" "Control mean baseline" "p-value MA=MD" )) 


eststo clear


foreach x of global main_results {
eststo:areg  `x'_95 `x'_95_base treatment2 treatment3 , absorb($strata3) vce(robust)
sum `x'_95 if treatment1==1
local control: di %9.2fc r(mean)
display "`control'"
estadd local control `control'
sum `x'_95_base if treatment1==1
local control_base: di %9.2fc r(mean)
estadd local control_base `control_base'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd local test `test'
}


esttab 					using "$analysis/output/appendix/table_robust_winsor.tex", append f ///
						posthead(\hline \multicolumn{4}{l}{\textit{Winsorizing 95\%}} \\) ///
						postfoot(\bottomrule ///
						\multicolumn{4}{p{0.7\textwidth}}{\footnotesize  Intent-to-treat estimates. All outcomes are unwinsorized. Values in USD. All regressions include strata dummies and include the baseline value of the outcome.  Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account.  Profits refers to the woman's self-reported monthly business profit. Savings is individual savings held by the woman. Capital is the value of all assets the woman uses in her business plus the value of inventory held for her business. Control mean endline is the mean value of the outcome in the control group at endline. Control mean baseline is the mean value of the outcome in the control group at baseline. Robust standard errors in parentheses. } ///
						\\ \end{tabular} \\ \end{table}) b(2) se(2) nostar ///
						keep(treatment2 treatment3 ) ///
						varlabels(treatment2 "Mobile Account" treatment3 "Mobile Disburse") ///
						booktabs nodep nonum nomtitles  ///
						noobs  collabels(none)  ///
						stats( control control_base test N, fmt(2 2 2 %9.0fc) ///
						label ("Control Mean" "Control mean baseline" "p-value MA=MD" "Observations")) 




//robustness
*control for unbalanced variables at baseline - Table A12
eststo clear
foreach x of global main_results {
eststo:  areg  `x' `x'_base treatment2 treatment3 mobile_account_base completed_secondary_base hyperbolic_base   , absorb($strata3) vce(robust)
sum `x' if treatment1==1
local control: di %9.4fc r(mean)
display "`control'"
sum `x'_base if treatment1==1
local control_base: di %9.4fc r(mean)
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd local test `test'

}

esttab 	using "$analysis/output/appendix/table_robust.tex", replace f ///
								prehead(\begin{table}[htbp]\centering \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi} \caption{Robustness Checks - controls \label{robustness}} \begin{tabular}{l*{@M}{c}} \toprule) ///
								posthead(\hline \multicolumn{4}{l}{\textit{Controlling for imbalanced variables at baseline}} \\) ///
								b(2) se(2) nostar  /// 
								keep(treatment2 treatment3 ) ///
								varlabels(treatment2 "Mobile Account" treatment3 "Mobile Disburse") ///
								label  booktabs noobs nonotes  collabels(none) ///
								mtitle("Profit" "Saving" "Capital" ) 


eststo clear								
//Including days 
foreach x of global main_results {
eststo: areg  `x' `x'_base treatment2 treatment3 $days1, absorb($strata3) vce(robust)
sum `x' if treatment1==1
local control: di %9.2fc r(mean)
display "`control'"
sum `x'_base if treatment1==1
local control_base: di %9.2fc r(mean)
test treatment2=treatment3
return list        
local test: di %9.4fc r(p)
}


esttab 	using "$analysis/output/appendix/table_robust.tex", append f ///
								posthead(\hline \multicolumn{4}{l}{\textit{Controlling for linear and quadratic time trend}} \\) ///
								b(2) se(2) nostar  /// 
								keep(treatment2 treatment3 ) ///
								varlabels(treatment2 "Mobile Account" treatment3 "Mobile Disburse") ///
								label  booktabs noobs nonum nomtitle nonotes  collabels(none) ///
		



eststo clear
*controling for take up correlates Table A12
foreach x of global main_results {
eststo: areg  `x' `x'_base treatment2 treatment3 married_base own_decision_base  , absorb($strata3) vce(robust)
sum `x' if treatment1==1
local control: di %9.2fc r(mean)
estadd local control `control'
sum `x'_base if treatment1==1
local control_base: di %9.2fc r(mean)
estadd local control_base `control_base'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd local test `test'
}

esttab 					using "$analysis/output/appendix/table_robust.tex", append f ///
						posthead(\hline \multicolumn{4}{l}{\textit{Controlling for correlates of takeup}} \\) ///
						postfoot(\bottomrule ///
						\multicolumn{4}{p{0.7\textwidth}}{\footnotesize  Intent-to-treat estimates. All outcomes are winsorized at the 99\% level. USD.  All regressions include strata dummies and include the baseline value of the outcome. The first panel controls for those variables imbalanced in Table \ref{balance}. The second panel controls for linear and quadratics of the number of days between loan disbursement and endline. The third panel controls for correlates of takeup in Table \ref{takeup2} column (2). Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account. Profits refers to the woman's self-reported monthly business profit. Savings is individual savings held by the woman. Capital is the value of all assets the woman uses in her business plus the value of inventory held for her business. Control mean endline is the mean value of the outcome in the control group at endline. Control mean baseline is the mean value of the outcome in the control group at baseline.  Robust standard errors in parentheses.  } ///
						\\ \end{tabular} \\ \end{table}) b(2) se(2) nostar ///
						keep(treatment2 treatment3 ) ///
						varlabels(treatment2 "Mobile Account" treatment3 "Mobile Disburse") ///
						booktabs nodep nonum nomtitles  ///
						noobs  collabels(none)  ///
						stats( control control_base N, fmt(2 2 %9.0fc) ///
						label ("Control Mean" "Control mean baseline" "Observations")) 


*logs Table A13
eststo clear 
foreach x of global main_results {
eststo est`x': areg  ln_`x' ln_`x'_base treatment2 treatment3 , absorb($strata3) vce(robust)
sum `x' if treatment1==1
local control: di %9.2fc r(mean)
estadd scalar control `control'
sum `x'_base if treatment1==1
local control_base: di %9.2fc r(mean)
estadd scalar control_base `control_base'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd scalar test `test'
}

esttab  using "${analysis}/output/appendix/table_primary_logs.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label ///
					stats(N r2 control control_base test, fmt(%9.0fc 2 2 2 ) label ("Observations" "R-squared"  "Control mean" "Control mean baseline" "p-value MA=MD")) /// 
					title(Treatment effects on woman's log business profits,  log savings and log business capital \label{results_logs}) postfoot(\bottomrule \multicolumn{4}{p{12cm}}{Intent-to-treat estimates. All outcomes are winsorized at the 99\% level and reported in logs. Observations vary due to the presence of zeros. All regressions include strata dummies and include the baseline value of the outcome. Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account.  Profits refers to the woman's self-reported monthly business profit. Savings is individual savings held by the woman. Capital is the value of all assets the woman uses in her business plus the value of inventory held for her business. Control mean endline is the mean value of the outcome in the control group at endline. Control mean baseline is the mean value of the outcome in the control group at baseline.   Robust standard errors in parentheses. }  \\ \end{tabular} \end{table}) replace substitute(\_ _) 


*Table A14 permutation test
if $perm==1{
clear matrix
*permutation tests
foreach x of global main_results {

	matrix define b`x' = (. , .)
	matrix colnames b`x' = treatment2 treatment3


eststo est`x': areg  `x' `x'_base treatment2 treatment3 , absorb($strata3) vce(robust)
sum `x' if treatment1==1
local control: di %9.2fc r(mean)
estadd scalar control `control'
sum `x'_base if treatment1==1
local control_base: di %9.2fc r(mean)
estadd scalar control_base `control_base'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd scalar test `test'


permute treatment2 _b[treatment2], rep(1000) strata(${strata3}) rseed(9999) : reg `x' treatment2 treatment3  
mat b`x'[1,1] = r(p_twosided) 
permute treatment3 _b[treatment3], rep(1000) strata(${strata3}) rseed(9999) : reg `x' treatment2 treatment3 
mat b`x'[1,2] = r(p_twosided) 

estadd matrix bp = b`x' : est`x'
}



esttab  using "${analysis}/output/appendix/table_primary_perm.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par) ///
				p( par([ ])) ///
				bp( par(\{ \}) pvalue(bp))) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label ///
					stats(N r2 control control_base test, fmt(%9.0fc 2 2 2 ) label ("Observations" "R-squared"  "Control mean" "Control mean baseline" "p-value MA=MD")) /// 
					title(Treatment effects on woman's business profits, savings and business capital - permutation test \label{results_perm}) postfoot(\bottomrule \multicolumn{4}{p{12cm}}{Intent-to-treat estimates. All outcomes are winsorized at the 99\% level. '000 Ugandan Shillings. All regressions include strata dummies and include the baseline value of the outcome. Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account.  Profits refers to the woman's self-reported monthly business profit. Savings is individual savings held by the woman. Capital is the value of all assets the woman uses in her business plus the value of inventory held for her business. Control mean endline is the mean value of the outcome in the control group at endline. Control mean baseline is the mean value of the outcome in the control group at baseline.  Permutation p-values are shown in curly brackets. These used the permute command in Stata and 1000 repetitions. Robust p-values in square brackets. Robust standard errors in parentheses. }  \\ \end{tabular} \end{table}) replace substitute(\_ _) 


}


  




*contamination bias
multe earn_business treatment , control($strata3) 
multe, decomp // small bias
multe, gen(lambda(M_) tau(tauhat_))
multe capital treatment , control($strata3 )
multe, decomp // small bias



//Table A21
if $forest==1 {
preserve
do "$analysis/do-files/causal_forest.do"
restore
}

//spouse present - Table A22
 gen treatment2_no_spouse_home=treatment2*no_spouse_home_base
 la var treatment2_no_spouse_home "Mobile Account * no spouse at home"
 gen treatment3_no_spouse_home=treatment3*no_spouse_home_base
 la var treatment3_no_spouse_home "Mobile Disburse * no spouse at home"

eststo clear 
foreach x of varlist earn_business capital  {
eststo: areg  `x' `x'_base treatment2 treatment3 treatment2_no_spouse_home treatment3_no_spouse_home no_spouse_home_base, absorb($strata3) vce(robust)
sum `x' if treatment1==1
local control: di %9.2fc r(mean)
estadd scalar control `control'
sum `x'_base if treatment1==1
local control_base: di %9.2fc r(mean)
estadd scalar control_base `control_base'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd scalar test `test'
test treatment2_no_spouse_home=treatment3_no_spouse_home
return list        
local test2: di %9.2fc r(p)
estadd scalar test2 `test2'

}

esttab  using "${analysis}/output/appendix/table_spouse_present.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment2 treatment3 treatment2_no_spouse_home treatment3_no_spouse_home no_spouse_home_base)     collabels(none)  nocons  eqlabels(none) booktabs label  ///
					stats(N r2 control control_base  test test2, fmt(%9.0fc 2  2 ) label ("Observations" "R-squared"  "Control mean" "Control mean baseline" "p-value MA=MD" "p-value MA=MD interaction")) /// 
					 title(Heterogeneous treatment effects by spouse presence \label{tab:spouse_presence}) postfoot(\bottomrule \multicolumn{3}{p{15cm}}{ Intent-to-treat estimates. Monetary outcomes are winsorized at the 99\% and in USD.  Mobile Account (MA) is the treatment where only a mobile money account was provided and the loan was disbursed as cash. Mobile Disburse (MD) is the treatment where a mobile money account was provided and the loan also disbursed onto this account. No spouse at home is a dummy variable equal to 1 if the woman was either not married or her spouse lived away from home at baseline. Profits refers to the woman's self-reported monthly business profit. Capital is the value of all assets the woman uses in her business plus the value of inventory held for her business. Means are shown at both baseline and endline in the control group. Robust standard errors in parentheses. } \\ \end{tabular} \end{table}) replace substitute(\_ _) 

					 
//Tables 23 and 24
preserve
do "$analysis/do-files/k-means.do"
restore


// Tables A25 (earnings)  A26 (happiness) A27 (empowerment) A28 (records)  A29 (remittance) A30 (group)
foreach root in   earnings happy empower_all records remittance           group   {
local i=1
foreach x of global `root' {
	capture confirm variable `x'_base 
	if !_rc {
		areg  `x' `x'_base treatment2 treatment3, absorb($strata3) vce(robust)
		sum `x' if treatment==0
		local control: di %9.4fc  r(mean)
		sum `x'_base if treatment==0
		local control_base: di %9.4fc  r(mean)
		test treatment2=treatment3
		local test=r(p)
	if `i'==1 {
outreg2 using "${analysis}/output/appendix/table_second_`root'.tex", keep(`x' treatment2 treatment3)  addstat(Control mean endline, `control', Control mean baseline, `control_base', p-value MA=MD, `test')  nocons dec(2)  replace label adec(2) noaster
}
else {
outreg2 using "${analysis}/output/appendix/table_second_`root'.tex", keep(`x' treatment2 treatment3)  addstat(Control mean endline, `control', Control mean baseline, `control_base', p-value MA=MD, `test') nocons dec(2)  label adec(2) noaster
}

	}
	
	else {
		areg  `x'  treatment2 treatment3, absorb($strata3) vce(robust)
		sum `x' if treatment==0
		local control=r(mean)
		test treatment2=treatment3
		local test=r(p)
if `i'==1 {
outreg2 using "${analysis}/output/appendix/table_second_`root'.tex", keep(`x' treatment2 treatment3)  addstat(Control mean endline, `control', p-value MA=MD, `test')  nocons dec(2)  replace label adec(2) noaster
}
else {
outreg2 using "${analysis}/output/appendix/table_second_`root'.tex", keep(`x' treatment2 treatment3)  addstat(Control mean endline, `control', p-value MA=MD, `test') nocons dec(2)  label adec(2) noaster
}

	}

local i=`i'+1
}
}





//loan outcomes Table A31
preserve
use "${analysis}/input/BRAC admin data", clear

eststo clear 
foreach x of varlist  missed_payment missed_days_brac principal_outstanding_brac interest_outstanding_brac savings_amt_brac overdue_amount_brac  {
eststo: areg  `x' treatment2 treatment3 , absorb($strata3) vce(robust)
sum `x' if treatment1==1
local control: di %9.2fc r(mean)
estadd scalar control `control'
test treatment2=treatment3
return list        
local test: di %9.2fc r(p)
estadd scalar test `test'
}

esttab  using "${analysis}/output/appendix/table_brac_admin.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment2 treatment3)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control  test, fmt(%9.0fc 2  2 ) label ("Observations" "R-squared"  "Control mean"  "p-value MA=MD")) /// 
					 prehead(\begin{table}[hbt] \centering \caption{Treatment effects on loan repayment}\label{tab:it_data} \begin{tabulary}{1\textwidth}{lCCCCCC} \hline  & (1) & (2) & (3) & (4) & (5) & (6) \\ & missed payment & missed days & principal outstanding & interest outstanding & savings amt & overdue amount \\ ) postfoot(\bottomrule \multicolumn{7}{p{16cm}}{Data from BRAC administrative records, hence sample is all 2959 baseline women. Intent-to-treat estimates. All regressions include strata dummies. Mobile Account is the treatment where only a mobile money account was provided and the loan was disbursed as cash.  Mobile Disburse is the treatment where a mo-bile money account was provided and the loan also disbursed onto this account. Missed payment is a dummy variable if a payment was not made the week it was due. Missed days is the number of days a payment is overdue, 0 for those without an overdue payment. Principal outstanding is the amount of loan still remaining to be paid, interest outstanding is the amount of interest remaining to be paid. Saving amount is the saving balance held by BRAC. Overdue amount is the amount  due for overdue payments, or 0 otherwise. Columns (3)-(6) are in USD. No winsorizing is applied to this data. Robust standard errors in parentheses.}   \\ \end{tabulary} \end{table}) replace substitute(\_ _) 

restore


	//Tables A32-A34 all heterogeneity
foreach x of global hetero   {
	gen treatment2_`x'=treatment2*`x'
		gen treatment3_`x'=treatment3*`x'
}

	*dimensions of hetero for main outcomes
	foreach x of global main_results {
	local i=1
	foreach y of global hetero {
	capture drop temp *temp
	gen temp=`y'
	gen treatment2_temp=treatment2*`y'
	label var treatment2_temp "MA*interaction"
	gen treatment3_temp=treatment3*`y'
	label var treatment3_temp "MD*interaction"
	sum `x' if treatment==0 & `y'==1 
	local c1: di %9.2f r(mean)
	display `c1'
	sum `x'_base if treatment==0 & `y'==1
	local c2: di %9.2f r(mean)
	display `c2'
	sum `y'
	local int: di %9.3f r(mean)
	display `int'
	areg `x'  `x'_base treatment2_temp treatment3_temp temp treatment2 treatment3, absorb($strata3) vce(robust)
	*reg `x'  `x'_base treatment2_temp treatment3_temp temp treatment2 treatment3 $strata1, robust
	test treatment2=treatment3
	local 	test: display %9.3f `r(p)' 
	display `test'
	test treatment2+treatment2_temp=treatment3+treatment3_temp
	local 	test_temp: display %9.3f `r(p)' 
	display `test_temp'

	if `i'==1 {
		outreg2 using "${analysis}/output/appendix/table_hetero_`x'.tex" , keep(`x' treatment2 treatment3 treatment2_temp treatment3_temp ) label ctitle(`y') ///
	addstat("Control mean", `c1' ,  "Interaction mean", `int',  "MA=MD", `test', "MA=MD interaction", `test_temp' )  nocons replace  dec(2)  adec(2) noaster
	 }
	else {
	 	outreg2 using "${analysis}/output/appendix/table_hetero_`x'.tex" , keep(`x' treatment2 treatment3 treatment2_temp treatment3_temp ) label ctitle(`y') ///
	addstat("Control mean", `c1' ,   "Interaction mean", `int',  "MA=MD", `test', "MA=MD interaction", `test_temp' )  nocons  dec(2) adec(2) noaster
	}
	local i=`i'+1	
	drop temp *temp
	}
	}

					 
					 
***********************************************************
***********************	TRANSACTION DATA*******************
**********************************************************


*MM transaction data  - Figure 2 

use "${analysis}/input/MM_balances.dta", clear
if $graph==1 {
graph twoway scatter balance_loan days if treatment==1, mcolor( blue%60) msize(vsmall) || scatter balance_loan days if treatment==2,  mcolor(red%60 ) msize(vsmall)  legend( order(1 "Mobile account" 2 "Mobile disbursement" ) )  plotregion(ilstyle(none) lcolor(white)) graphregion(color(white) lc(white) lw(med)) bgcolor(white) ylabel(,grid) title("Average mobile money account daily balance")  subtitle("% of initial loan value")
graph export "$analysis/output/graphs/balances.png", replace
}


*Tables A17 and A20, Table A15 and Table A16 useage tables and Figure A6 
use "${analysis}/input/MM_useage", clear

//summary - Table A15
eststo clear
bysort treatment: eststo: ///
estpost sum  ever_deposit ever_withdrawal number_deposit number_withdrawal av_deposit av_withdrawal total_deposit total_withdrawal withdrew_perc withdrew_disburse, detail
esttab using "$analysis/output/appendix/table_int_summary.tex", replace nonum nomtitle collabels("mean" "sd" "median")  cells(" mean(fmt(2)) sd(fmt(2))   p50(fmt(2))")  label prehead(\begin{table}[htbp] \centering \caption{Summary statistics of mobile money account usage - compliers only}  \label{int_sum} \resizebox{1\textwidth}{!}{   \begin{tabulary}{1.45\textwidth}{Lcccccc} \hline  & \multicolumn{3}{c}{Mobile account} & \multicolumn{3}{c}{Mobile disburse} \\ \cmidrule(r){2-4}  \cmidrule(r){5-7} \hline  ) postfoot(\bottomrule \multicolumn{7}{p{17cm}}{Monetary outcomes are in USD. Includes all women who received a sim card as part of the study in Mobile Account and all women who received their loan on the mobile money account in Mobile Disbursement - excludes women in each treatment arm who did not receive the assigned treatment (non-compliers). All variables are defined over the first 180 days after the account was provided. I cap transactions at 180 since the last mobile money accounts were given out in June 2017 and the administrative data ends in January 2018. Deposits always excludes the loan disbursement for the mobile disbursement treatment group. Ever deposit and withdraw are dummy variables if at least one transaction of that type occurred. Number of deposits and withdrawals is the count of each transaction for an account.  Deposit amount and withdrawal amount summarises the mean transaction amount if that type of transaction occurred. Total deposits and withdrawals are cumulative transactions on an account.  Withdrew day 1 and \% loan withdrew loan day 1 are only captured for the Mobile Disbursement group and capture whether the woman withdrew any of the loan the day it was disbursed and what percentage of the loan she withdrew the day the loan was disbursed.} \end{tabulary}} \end{table}) substitute(| $|$)
eststo clear

//Figure A6
if $graph==1 {
reg total_deposit   month2-month6  if treatment==1,  noconstant
estimates store model_1 
coefplot model_1 , ytitle("Total deposits to account by month of opening USD") xtitle("Month of loan disbursement") vertical      mfcolor(blue%60) ylabel(,grid) graphregion(color(white)) bgcolor(white)
graph export "$analysis/output/appendix/graphs/MA_time2.png", replace

}


//Table A16
eststo clear 
foreach x of varlist  ever_deposit number_deposit av_deposit total_deposit ever_withdrawal  number_withdrawal av_withdrawal total_withdrawal  {
eststo: areg  `x' treatment_2  , absorb($strata3) vce(robust)
sum `x' if treatment==1
local control: di %9.2fc r(mean)
estadd scalar control `control'
}

esttab  using "${analysis}/output/appendix/table_int3_results.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment_2)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control , fmt(%9.0fc 2 2 ) label ("Observations" "R-squared"  "Mobile account mean" )) /// 
					 prehead(\begin{table}[htb] \centering \caption{Treatment effects on mobile money usage outcomes} \label{int3_results} \resizebox{1\textwidth}{!}{ \begin{tabulary}{1.4\textwidth}{lCCCCCCCC} \hline  & (1) & (2) & (3) & (4) & (5) & (6)  & (7) & (8) \\  & Ever deposit & Number deposit & Average deposit & Total deposit  & Ever withdraw  & Number withdrawals  & Average withdrawal & Total withdrawals \\ ) postfoot(\bottomrule \multicolumn{9}{p{22cm}}{Impacts amongst those who received sim cards (Mobile Account) and received the loan on the account (Mobile Disburse). All regressions include strata dummies.  Monetary outcomes in USD. All variables are defined over the first 180 days after the account was provided. MD (Mobile Disburse) is the treatment where a mobile money account was provided and the loan also disbursed onto this account.  Mobile account mean refers to the mean in the mobile account group. Robust standard errors in parentheses. }   \\ \end{tabulary}} \end{table}) replace substitute(\_ _) 

					 
//Table A17
eststo clear 
foreach x of varlist  av_balance7  av_balance15 av_balance30  av_balance45 av_balance60 av_balance90 av_balance180 final_balance  {
eststo: areg  `x' treatment_2  , absorb($strata3) vce(robust)
sum `x' if treatment==1
local control: di %9.2fc r(mean)
estadd scalar control `control'
}

esttab  using "${analysis}/output/appendix/table_int2_results.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment_2)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control , fmt(%9.0fc 2 2 ) label ("Observations" "R-squared"  "Mobile account mean" )) /// 
					 prehead(\begin{table}[htb] \centering \caption{Treatment effects on mobile money balances} \label{int2_results} \resizebox{1\textwidth}{!}{\begin{tabulary}{1.3\textwidth}{LCCCCCCCC} \hline  & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) \\ & Average balance 0-7 & Average balance 8-15 & Average balance 15-30 & Average balance 30-45 & Average balance 45-60 & Average balance 60-90 & Average balance 90-180 & Final balance \\) postfoot(\bottomrule \multicolumn{9}{p{20cm}}{Impacts amongst those who received sim cards (Mobile Account) and the loan on the account (Mobile Disburse) - excludes those who didn't receive the treatment as assigned (non-compliers). Average balance in USD. All regressions include strata dummies. Mobile Disburse is the treatment where a mobile money account was provided and the loan also disbursed onto this account. Average balance is the average end of day balance on the account the specified number of days after the account was given to the client. Final balance is the balance at the last transaction made within 180 days of account opening. Mobile account mean refers to the mean in the mobile account group. Robust standard errors in parentheses. }   \\ \end{tabulary}} \end{table}) replace substitute(\_ _) 


//Table A20
local hetero hetero_family_median
foreach x of local hetero {
gen `x'_treat_2=`x'*treatment_2
}
la var hetero_family_median_treat_2 "MD*family pressure"

eststo clear 
foreach x of varlist  av_balance7  av_balance15 av_balance30  av_balance45 av_balance60 av_balance90 av_balance180 final_balance  {
eststo: areg  `x' treatment_2  hetero_family_median_treat_2 hetero_family_median , absorb($strata3) vce(robust)
sum `x' if treatment==1
local control: di %9.2fc r(mean)
estadd scalar control `control'
}

esttab  using "${analysis}/output/appendix/table_int_hetero.tex" , ///
					cells(b(fmt(%12.2fc)) se(fmt(%12.2fc) par)) ///
					nostar ///
					keep (treatment_2 hetero_family_median_treat_2 hetero_family_median)     collabels(none)  nocons  eqlabels(none) booktabs label nonum nomtitles ///
					stats(N r2 control , fmt(%9.0fc 2 2 ) label ("Observations" "R-squared"  "Mobile account mean" )) /// 
					 prehead(\begin{table}[htb] \centering \caption{Treatment effect heterogeneity by sharing pressure on mobile money account usage outcomes} \label{transaction_hetero} \resizebox{1\textwidth}{!}{\begin{tabulary}{1.4\textwidth}{LCCCCCCCC} \hline & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) \\ & Average balance 0-7 & Average balance 8-15 & Average balance 15-30 & Average balance 30-45 & Average balance 45-60 & Average balance 60-90 & Average balance 90-180 & Final balance \\) postfoot(\bottomrule \multicolumn{9}{p{22cm}}{Impacts amongst those who received sim cards (Mobile Account) and received the loan on the account (Mobile Disburse). All regressions include strata dummies.  Monetary outcomes in USD. All variables are defined over the first 180days after the account was provided. MD (Mobile Disburse) is the treatment where a mobile money account was provided and the loan also disbursed onto this account.  Control mean refers to the mean in the mobile account group. Robust standard errors in parentheses. }   \\ \end{tabulary}} \end{table}) replace substitute(\_ _) 




//Table A18 and Figure A5
use "${analysis}/input/MM_transactions.dta", clear

*figure A5
if $graph==1 {
graph twoway (hist transfer_type if treatment==1, percent xlabel(1 "Cash in"  2 "Cash out" 3 "Debit" 4 "payment" 5 "reversal" 6 "Transfer In" 7 "Transfer Out") color(blue%30)  xtitle("Transfer type") ) (hist transfer_type if treatment==2 , percent xlabel(1 "Cash in"  2 "Cash out" 3 "Debit" 4 "payment" 5 "reversal" 6 "Transfer In"  7 "Transfer Out")  xtitle("Transfer type") color(red%30)), legend( order(1 "Mobile account" 2 "Mobile disbursement" ) ) graphregion(margin(5 10 5 5) fcolor(white) lcolor(white)) plotregion(fcolor(white) ilstyle(none) icolor(white) lcolor(white))
graph export "$analysis/output/appendix/graphs/trans_type.png", replace
}

*Table A18
eststo clear
bysort treatment transfer_type:  sum  AMOUNT 
eststo ma: estpost tabstat AMOUNT if  treatment==1, by(transfer_type) statistics(mean sd n) nototal
eststo md: estpost tabstat AMOUNT if treatment==2, by(transfer_type) statistics(mean sd n) nototal

esttab ma md using "$analysis/output/appendix/table_trans_type.tex", replace cells("mean(pattern(1 1) fmt(2)) sd(pattern(1 1) fmt(2)) count(pattern(1 1 1) fmt(%9.0fc))") label mtitles("Mobile Account" "Mobile Disburse") noobs booktabs stats(N, fmt(%9.0gc) label("Observations")) title( Mobile money account transaction value by transaction type and treatment group \label{tab:transfer_type}) substitute(\_ _ {l} {p{11cm}}) addnotes(Transaction value USD. Excludes loan deposit. Cash in refers to depositing cash using a mobile money agent. Cash out is withdrawing cash through a mobile money agent. Debit is a transfer from another mobile money account or bank. Payment is principally buying airtime or data. Transfer refers to sending/receiving money to/from another mobile money account or bank account.)



