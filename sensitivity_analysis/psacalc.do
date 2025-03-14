	   
program psacalc, rclass
	version 11.2
	syntax namelist(min=2 max=2) [, mcontrol(string) rmax(real 1.0) delta(real 1.0) beta(real 0) model(string) weight(string)]
	tokenize `namelist'
	local treatment "`1'"
	local type "`2'"
	
	tempvar yhat4566 yvar4533 A4533 B4533 C4533 D4533 theta4533 sample xvar4533 tau4533 xhat4566 aa4533 bb4533 cc4533 R4533 Q4533 sol_4533_1 cubicyhat cubicwhat sol1_poss sol2_poss sol3_poss dist_1 dist_2 dist_3 mindist M4533 xhat4567
	
	if "`model'"=="" {
		quietly `e(cmdline)'
	}
	
	else {
		quietly `model'
	}
			
	
	/************ Error cases -> exit program *******************/
		
	// No previous regression, and failed to specify linear model
	if e(cmd)~="regress" {
		di as error _n "This command only works after linear regression, or else must specify full linear regression model in options"
		exit 5001				
	}
	
	// rmax out of bounds
	else if `rmax'>1 {
		di as err _n "The maximum possible R-squared is 1"
		exit 5001
	} 
	
	// rmax less than controlled R-squared
	else if `rmax'<e(r2) {
		di as err _n "Maximum R-squared provided is less than controlled R-squared"
		exit 5001
	}
	
	// bad type
	else if "`type'"~="set" & "`type'"~="beta" & "`type'"~="delta" & "`type'"~="rmax" {
		di as err _n "Invalid type: must be set, beta, delta, or rmax"
		exit 5001
	}
	
	// specified treatment not in the model
	local names: colnames e(b)	
	else if strpos("`names'","`treatment'")==0 {

		// treatment is not in the model at all
		if strpos(e(cmdline),"`treatment'")==0 {
			di as err _n "Specified treatment variable missing from model"
			exit 5001
		}
		
		// treatment is not an independent variable in the model
		else {
			di as err _n "Specified treatment variable is not an indep var in the model"
			exit 5001
		}
	}
	
	
		/************* Calculations (main program arc) ***************/
	
	else {
		preserve
		
		local command=e(cmdline)
		tokenize `command'
		macro shift 3
		local rest `*'
		tokenize `command'	
		macro shift 1
		local all `*'
		tokenize `command'
	
		
		
		quietly `command' 
		quietly predict `sample' if e(sample)
		
			
		if "`weight'"~="" {
			quietly reg `2' `treatment' `mcontrol' [`weight'] if `sample'~=.
		}
		
		else {
	
			 quietly reg `2' `treatment' `mcontrol' if `sample'~=.
		}
		
		local hat_beta=_b[`treatment']
		local hat_r=e(r2)
		
		quietly `command'
		
		local tilde_beta=_b[`treatment']
		local tilde_r=e(r2)
		
	
		
			// Define variance Terms
		
		quietly sum `e(depvar)'
		quietly gen `yvar4533'=r(Var)
		quietly reg `treatment' `mcontrol'
		quietly predict `xhat4567', resid
		
		quietly sum `xhat4567'
		quietly gen `xvar4533'=r(Var)
		
		quietly reg `treatment' `rest'
		quietly predict `xhat4566', resid
		quietly sum `xhat4566'
		quietly {
		gen `tau4533'=r(Var) 
		
		
		gen `A4533'=(`hat_beta'-`tilde_beta')
		gen `B4533'=(`tilde_r'-`hat_r')*`yvar4533'
		gen `C4533'=(`rmax'-`tilde_r')*`yvar4533'
		}

		
		// Solve for Delta for Beta = beta-hat
		
			
		local boundx=((`tilde_beta'-`beta')*`B4533'*`tau4533'+(`tilde_beta'-`beta')*`xvar4533'*`tau4533'*(`A4533')^2+2*(`tilde_beta'-`beta')^2*(`tau4533'*`A4533'*`xvar4533')+(`tilde_beta'-`beta')^3*(`tau4533'*`xvar4533'-`tau4533'^2))/(`C4533'*`A4533'*`xvar4533'+(`tilde_beta'-`beta')*`C4533'*(`xvar4533'-`tau4533')+(`tilde_beta'-`beta')^2*(`tau4533'*`A4533'*`xvar4533')+(`tilde_beta'-`beta')^3*(`tau4533'*`xvar4533'-`tau4533'^2))
	


	if `delta'==1 {
	
	quietly {
		gen `aa4533' = `C4533'*(`xvar4533'-`tau4533')-`B4533'*`tau4533'-`xvar4533'*`tau4533'*(`A4533')^2
		gen `bb4533'=4*`C4533'*(`A4533')^2*(`xvar4533')^2*`tau4533'
		gen `cc4533'=-2*`tau4533'*(`A4533')*`xvar4533'
		
		local solution_1=(-1*`aa4533'-sqrt((`aa4533')^2+`bb4533'))/(`cc4533')
		local solution_2=(-1*`aa4533'+sqrt((`aa4533')^2+`bb4533'))/(`cc4533')
			forvalues i=1/2 {
		
			local beta_`i' = `tilde_beta'-`solution_`i''
			gen `dist_`i'' = (`beta_`i''-`tilde_beta')^2
			replace `dist_`i''=999999 if sign(`beta_`i''-`tilde_beta')~=sign(`tilde_beta'-`hat_beta')
		}
		
	egen `mindist'=rowmin(`dist_1' `dist_2')
	local j = 2
		gen temp45041=.
		gen temp45042=.
	
			forvalues i=1/2 {
			if `dist_`i''==`mindist' {
				local betax=`beta_`i''
			}
				else if `dist_`i''~=`mindist' {
				replace temp4504`i'=`beta_`i''
			}
			
		}
		
		if temp45041~=. {
		
		local altroot1=temp45041
		}
		if temp45041==. & temp45042~=. {
		local altroot1=temp45042 
		}
		
		drop temp4504*
	
	
	}
	}
	
	else if `delta'~=1 {
	


	quietly {
	gen `aa4533'=(`tau4533'*`A4533'*`xvar4533'*(`delta'-2))/((`delta'-1)*(`tau4533'*`xvar4533'-`tau4533'^2))
	gen `bb4533'=(`delta'*`C4533'*(`xvar4533'-`tau4533')-`B4533'*`tau4533'-`xvar4533'*`tau4533'*`A4533'^2)/((`delta'-1)*(`tau4533'*`xvar4533'-`tau4533'^2))
	gen `cc4533'= (`C4533'*`delta'*`A4533'*`xvar4533')/((`delta'-1)*(`tau4533'*`xvar4533'-`tau4533'^2))
	
	
	
	gen `Q4533' = (`aa4533'^2-3*`bb4533')/9
	gen `R4533' = (2*`aa4533'^3-9*`aa4533'*`bb4533'+27*`cc4533')/54
	gen `M4533' =  `R4533'^2-`Q4533'^3
	local discrim = `R4533'^2-`Q4533'^3
	}
		

		// Solve For Cubic Roots if delta neq 1	
		
		
		
	if `discrim' < 0 {
	
	
	
	gen `theta4533' = acos(`R4533'/sqrt(`Q4533'^3))
	
	local solution_1 = -2*sqrt(`Q4533')*cos(`theta4533'/3)-(`aa4533'/3)
	local solution_2= -2*sqrt(`Q4533')*cos((`theta4533'+2*_pi)/3)-(`aa4533'/3)
	local solution_3= -2*sqrt(`Q4533')*cos((`theta4533'-2*_pi)/3)-(`aa4533'/3)		
	
		forvalues i=1/3 {
		
			local beta_`i' = `tilde_beta'-`solution_`i''
			local cov_`i' = `solution_`i''*(`xvar4533'-`tau4533')+`A4533'*`xvar4533'
			gen `dist_`i'' = (`beta_`i''-`tilde_beta')^2
			
			replace `dist_`i''=999999 if sign(`beta_`i''-`tilde_beta')~=sign(`tilde_beta'-`hat_beta')
		}
		
	egen `mindist'=rowmin(`dist_1' `dist_2' `dist_3')
	local j = 3
			gen temp45041=.
			gen temp45042=.
			gen temp45043=.
			
			forvalues i=1/3 {
			if `dist_`i''==`mindist' {
				local betax=`beta_`i''
			}
				else if `dist_`i''~=`mindist' {
				replace temp4504`i'=`beta_`i''
			}
			
		}
		
		
		
		if temp45041~=. {
		local altroot1=temp45041
		}
		
		if temp45041==. & temp45042~=. {
		local altroot1=temp45042 
		}
		
		if  temp45041~=. & temp45042~=. {
		local altroot2=temp45042 
		}
		
		if temp45041~=. & temp45042==. {
		local altroot2 = temp45043
		}
		
		drop temp4504*
		
		
		

	}

	


 // Case with 1 Real Root 
	else if `discrim'>0 {

	
	local t1=-1*`R4533'+sqrt(`M4533')
	local t2=-1*`R4533'-sqrt(`M4533')
	

		local solution_1 = `t1'^(1/3)+`t2'^(1/3)-(`aa4533'/3)
		
	local betax=`tilde_beta'-`solution_1'
	local j=1
	}
}			


// OUTPUT
	
			// Case 2: beta
		if "`type'"=="beta" {
		
			di _n as txt ///
			_col(18) "{hline 4} Treatment Effect Estimate {hline 4}" _n ///
			_col(1) "{hline 13}{c +}{hline 64}" _n ///
			_col(1) "beta" _col(14) "{c |}" _col(18) as result %11.5f `betax' _n ///
			as txt _col(1) "{hline 13}{c +}{hline 64}" 	
		
			di _n as txt ///
			_col(18) "{hline 4} Inputs from Regressions {hline 4}" _n ///
			_col(14) "{c |}" _col(21) "Coeff." _col(49) "R-Squared" _n ///
			_col(1) "{hline 13}{c +}{hline 64}" _n ///
			_col(1) "Uncontrolled" _col(14) "{c |}" _col(18) as res %12.5f `hat_beta' _col(49) %5.3f `hat_r' _n ///
			_col(1) as txt "Controlled" _col(14) "{c |}" _col(18) as res %12.5f `tilde_beta' _col(49) %5.3f `tilde_r' _n ///
			_col(1) as txt "{hline 13}{c +}{hline 64}" 

			di _n as txt ///
			_col(18) "{hline 4} Other Inputs {hline 4}" _n ///
			_col(1) "{hline 13}{c +}{hline 64}" _n ///
			_col(1) "R_max" _col(14) "{c |}" _col(18) %5.3f `rmax' _n ///
			_col(1) "Delta" _col(14) "{c |}" _col(18) %5.3f `delta' _n ///
			_col(1) "M Controls" _col(14) "{c |}" _col(18) "`mcontrol'"  _n ///
			_col(1) "{hline 13}{c +}{hline 64}"  

				return scalar output=`betax'
				return scalar root_count=`j'
			

				
				
			if `j'==2 {
				di as txt _col(5) "Note: There is an alterantive solution, which is: "   `altroot1' 
			}
			
			if `j'==3 {
				di as txt _col(5) "Note: There are alterantive solutions, which are: "   `altroot1' " , " `altroot2'
			}
			
		}
	
			// Case 3: delta
		if "`type'"=="delta" {
		
			di _n as txt ///
			_col(18) "{hline 4} Bound Estimate {hline 4}" _n ///
			_col(1) "{hline 13}{c +}{hline 64}" _n ///
			_col(1) "delta" _col(14) "{c |}" _col(18) as result %11.5f `boundx' _n ///
			as txt _col(1) "{hline 13}{c +}{hline 64}" 	
		
			di _n as txt ///
			_col(18) "{hline 4} Inputs from Regressions {hline 4}" _n ///
			_col(14) "{c |}" _col(21) "Coeff." _col(49) "R-Squared" _n ///
			_col(1) "{hline 13}{c +}{hline 64}" _n ///
			_col(1) "Uncontrolled" _col(14) "{c |}" _col(18) as res %12.5f `hat_beta' _col(49) %5.3f `hat_r' _n ///
			_col(1) as txt "Controlled" _col(14) "{c |}" _col(18) as res %12.5f `tilde_beta' _col(49) %5.3f `tilde_r' _n ///
			_col(1) as txt "{hline 13}{c +}{hline 64}" 

			di _n as txt ///
			_col(18) "{hline 4} Other Inputs {hline 4}" _n ///
			_col(1) "{hline 13}{c +}{hline 64}" _n ///
			_col(1) "R_max" _col(14) "{c |}" _col(18) %5.3f `rmax' _n ///
			_col(1) "Beta" _col(14) "{c |}" _col(18) %9.6f `beta' _n ///
			_col(1) "M Controls" _col(14) "{c |}" _col(18) "`mcontrol'"  _n ///
			_col(1) "{hline 13}{c +}{hline 64}"  
			
			di as txt _col(5) "Reported delta matches a treatment effect of " as result `beta' 
			
		
			return scalar output=`boundx'
		}

	restore 
	
}		
	

 	
end

		