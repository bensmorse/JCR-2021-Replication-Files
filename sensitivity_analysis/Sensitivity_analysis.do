clear
//cd "\\ad.ucl.ac.uk\homee\uctqswe\Documents\REPLICATION\REPLICATION"
cd "C:\Users\Ben Morse\Dropbox (MIT)\IDP_Hosting_Paper\Syria\Data\REPLICATION"
use "stata_file.dta",clear
cap do "psacalc.do"

gl outcome1 currently_hosting hosted_number_duration
gl ctrls i.education_before_conflict i.work_before_conflict prev_res_land_type i.prev_res_who_own i.prev_res_struc_type ethnic respondent_gender age
cap matrix drop A


foreach y of varlist $outcome1 {


// naive -- no FE, no ctrls
// obtain naive beta, SE, r-squared 
qui xi: reg `y' violence_index_std , cl(current_village)
local rsqnaive = `e(r2)'
matrix coef = e(b)
local betanaive = coef[1,1]
matrix coefvar = e(V)
local senaive = sqrt(coefvar[1,1])

// obtain fully controlled r-squared and se
qui xi: reg `y' violence_index_std $ctrls, cl(current_village)
local rsqfull = `e(r2)'
matrix coef = e(b)
local betafull = coef[1,1]
matrix coefvar = e(V)
local sefull = sqrt(coefvar[1,1])


//local rmax0 = `rsqfull' + `rsqfull'-`rsqnaive' // sets rmax rsqfull plus the influence of the observed after conditioning on villages effects
local rmax1 = `rsqfull'*1.3  // sets rmax to experimental benchmark recommended by Oster
local rmax2 = `rsqfull'*1.6  // sets rmax to twice the experimental benchmark recommended by Oster
local rmax3 = `rsqfull'*2 // sets rmax to twice the fully controlled regression


// Delta, the effect of included control variables on treatment relative to omitted variables on treatment, is set to 1, following Oster 2014. The assumption is that the unobservables allowed to be as important as the observables in explaining the treatment.

// mcontrol are controls that should be maintained throughout-- i.e. those after which treatment should be exogenous
psacalc  violence_index_std   beta, rmax(`rmax1') delta(1) model(xi: reg `y' violence_index_std $ctrls) 
	local betarobust1 = `r(output)'
psacalc  violence_index_std  beta, rmax(`rmax2') delta(1) model(xi: reg `y' violence_index_std $ctrls) 
	local betarobust2 = `r(output)'
psacalc  violence_index_std  beta, rmax(`rmax3') delta(1) model(xi: reg `y' violence_index_std $ctrls) 
	local betarobust3 = `r(output)'
	
// store results

	matrix m`y'  = (`betanaive',`senaive',`rsqnaive',`betafull',`sefull',`rsqfull',`betarobust1',`betarobust2',`betarobust3')
	matrix list m`y'
	matrix A = nullmat(A) \ m`y'

}
	* Set column names 
	local cnames `" "B naive" "SE naive" "RSQ naive" "B full" "SE full" "RSQ full" "B Rmax=1.3 x RSQ full" "B Rmax=1.6 x RSQ full" "B Rmax=2 x RSQ full" "' 
	
	* Set row names
	local rowname
	foreach y of varlist $outcome1 {

			local rowname `" `rowname' "`y'" "' //"	
	}	

	
matrix list A
esttab matrix(A), tex























	