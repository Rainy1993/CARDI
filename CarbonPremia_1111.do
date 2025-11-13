// remove the concern of autoregression
// block bootstrap
// oos RME


clear
clear matrix
clear mata
global date = "20251015"
global root = "/Users/ruting/Documents/macbook/PcBack/23.Topic_wkh_cooperation/FRM/Output/EfficiencyTest/"
cd $root

set maxvar  120000

/**************** 1. data generating ****************/
use InsampleBeta_Premia_1014_month.dta, clear

tsset Date

* 控制变量
global control "Change_TY3M Slope TED RealEstate_excess MKreturn MKvol CarbonVol_Shenzhen CarbonVol_Guangdong CarbonVol_Hubei"

/********************
  2. 生成滞后变量
*********************/
local laglist 5 21 63 84 105 126
foreach lag of local laglist {
    gen L`lag'_FRM_5_High_Low = L`lag'.FRM_5_High_Low
    gen L`lag'_FRM_1_High_Low = L`lag'.FRM_1_High_Low
    gen L`lag'_FRM_10_High_Low = L`lag'.FRM_10_High_Low
}

save InsampleBeta_Premia_1014_lag.dta, replace

* 删除滞后缺失

/********************
  2. 双重循环 block bootstrap
*********************/

* Step 1. 创建总结果文件
tempfile allres
postfile allhandle lag frm mean_beta se_beta t_stat p_value using `allres', replace

* Step 2. 定义回归程序
capture program drop myreg
program define myreg, rclass
    args lag frm
    quietly regress Return L`lag'_FRM_`frm'_High_Low $control
    return scalar beta_CARDI = _b[L`lag'_FRM_`frm'_High_Low]
end

* Step 3. 双循环
foreach lag in 63 84 105 {
    foreach frm in 5 {
        use InsampleBeta_Premia_1014_lag.dta, clear
        tsset Date

        local frmvar L`lag'_FRM_`frm'_High_Low
        local varlist `frmvar' $control

        * 删除缺失观测
        gen dropflag = 0
        foreach var of local varlist {
            replace dropflag = 1 if missing(`var')
        }
        drop if dropflag==1
        drop dropflag

        * block分组
        sort Date
        local block = 3
        gen block_id = ceil((_n-1)/`block') + 1
        count
        local nblock = ceil(r(N)/`block')

        * Bootstrap采样
        tempfile bootres
        postfile handle beta using `bootres', replace

        set seed 12345
        forvalues i = 1/500 {
            preserve
            levelsof block_id, local(blocks)
            local sample_blocks ""
            forvalues j=1/`nblock' {
                local blk = word("`blocks'", ceil(runiform()*`=wordcount("`blocks'")'))
                local sample_blocks "`sample_blocks' `blk'"
            }
            gen keep_row = 0
            foreach b of local sample_blocks {
                replace keep_row = 1 if block_id==`b'
            }
            keep if keep_row==1
            drop keep_row

            myreg `lag' `frm'
            post handle (r(beta_CARDI))
            restore
        }

        postclose handle

        * 计算统计量
        use `bootres', clear
	    local filename = "boot_lag`lag'_frm`frm'.dta"
		save "`filename'", replace
		
        quietly summarize beta
        local beta_mean = r(mean)
        local beta_se   = r(sd)
        local tstat     = `beta_mean'/`beta_se'
        local pval      = 2*ttail(`=_N-1', abs(`tstat'))

        display "Lag `lag'  FRM_`frm' | β = " %6.4f `beta_mean' "  SE = " %6.4f `beta_se' ///
                "  t = " %6.3f `tstat' "  p = " %6.3f `pval'

        * 保存到总结果文件
        post allhandle (`lag') (`frm') (`beta_mean') (`beta_se') (`tstat') (`pval')
    }
}

postclose allhandle

/********************
  3. 导出结果表
*********************/
use `allres', clear

order lag frm mean_beta se_beta t_stat p_value
format mean_beta se_beta t_stat p_value %9.4f
export excel using "Bootstrap_FRM_results.xlsx", replace firstrow(variables)











/********************
  OOS 
*********************/
*------------------------------------------------------------*
* 0. 初始化
*------------------------------------------------------------*
/**************************************************
  Sample-out-of-sample + Diebold-Mariano test
  with multiple lags and FRM indicators
**************************************************/

use InsampleBeta_Premia_1014_lag.dta, clear
tsset Date

* 控制变量
global control "Change_TY3M Slope TED RealEstate_excess MKreturn MKvol CarbonVol_Shenzhen CarbonVol_Guangdong CarbonVol_Hubei"

* 样本外预测起点
local T0 = 120
local N_obs = _N
local start = `T0' + 1

* 创建 frame 保存结果
frame drop DMres
frame create DMres lag mse_base mse_model diff tstat pval

foreach lag in 63 84 105 {
    di as result "==============================="
    di as result "Running lag = `lag'..."
    di as result "==============================="

    preserve
    * 生成存放预测值的变量
    gen double yhat_model = .
    gen double yhat_base  = .

    *-----------------------------
    * 2. 样本外滚动预测
    *-----------------------------
    forvalues t = `start'/`N_obs' {
        quietly regress Return L63_FRM_5_High_Low $control if _n < `t'
        quietly predict double yhat_temp, xb
        replace yhat_model = yhat_temp if _n == `t'
        drop yhat_temp

        quietly regress Return $control if _n < `t'
        quietly predict double yhat_base_temp, xb
        replace yhat_base = yhat_base_temp if _n == `t'
        drop yhat_base_temp
    }

    *-----------------------------
    * 3. 计算预测误差平方
    *-----------------------------
    gen double e_model2 = (Return - yhat_model)^2
    gen double e_base2  = (Return - yhat_base)^2

    quietly summarize e_model2
    local rmse_model = sqrt(r(mean))
    quietly summarize e_base2
    local rmse_base = sqrt(r(mean))

    di as txt "RMSE(model) = " %6.4f `rmse_model'
    di as txt "RMSE(base)  = " %6.4f `rmse_base'
    di as txt "RMSE reduction = " %6.4f 1-(`rmse_model'/`rmse_base')

    *-----------------------------
    * 4. Diebold-Mariano 测试
    *-----------------------------
    dmariano Return yhat_base yhat_model, maxlag(12) kernel(bartlett)
	matrix list r(table)
	
	
    matrix res = r(table)

    local mse_base  = res[1,1]
    local mse_model = res[2,1]
    local diff      = res[3,1]
    local tstat     = res[4,1]
    local pval      = res[5,1]

    di as txt "DM results: mse_base=" %6.4f `mse_base' ///
        ", mse_model=" %6.4f `mse_model' ///
        ", diff=" %6.4f `diff' ///
        ", t=" %6.3f `tstat' ///
        ", p=" %6.3f `pval'

    *-----------------------------
    * 5. 保存到 DMres frame
    *-----------------------------
    frame DMres {
        set obs `=_N' + 1
        local last = _N
        replace lag       = `lag'       in `last'
        replace mse_base  = `mse_base'  in `last'
        replace mse_model = `mse_model' in `last'
        replace diff      = `diff'      in `last'
        replace tstat     = `tstat'     in `last'
        replace pval      = `pval'      in `last'
    }

    restore
}

*-----------------------------
* 6. 导出结果
*-----------------------------
frame DMres
list
save DM_results.dta, replace
export excel DM_results.xlsx, firstrow(variables) replace




// 定义变量组
local varlist "VaRRMAY_95 VaREMAY_95 VaREMAY_99 VaRRMAY_99 VaRRMAY_90 VaREMAY_90 VaRRMAY_5 VaREMAY_5 VaRRMAY_1 VaREMAY_1 VaRRMAY_10 VaREMAY_10 CResidual"

foreach ivar of local varlist {
	replace `ivar' = `ivar'*100
}
	
//local varlist "Return"

// 外层循环：遍历所有风险度量变量
foreach ivar of local varlist {
	
	
	// 内层循环：遍历所有滞后期
	foreach lag in 5 21 63 84 105 126 {

		// 定义不同滞后期的 FRM 指标
		local FRM_1  L`lag'_FRM_1_High_Low
		local FRM_5  L`lag'_FRM_5_High_Low
		local FRM_10 L`lag'_FRM_10_High_Low

		// 定义输出路径（确保 $save/$date 已定义）
		local outreg_file0 "$date/`ivar'_1015_`lag'.xls"

		// ---------- 回归 1 ----------
		reg `ivar' `FRM_1' $control
		outreg2 using "`outreg_file0'", ///
			stat(coef se) bdec(3) sdec(3) addstat(Adj. R-squared, e(r2_a)) ///
			title("Lag `lag' - FRM_1") replace drop(_I** _est* o.**) ///
			addtext(Controls, YES, Firm FE, YES, Y-var, `ivar', Cluster, Firm)

		// ---------- 回归 2 ----------
		reg `ivar' `FRM_5' $control
		outreg2 using "`outreg_file0'", ///
			stat(coef se) bdec(3) sdec(3) addstat(Adj. R-squared, e(r2_a)) ///
			title("Lag `lag' - FRM_5") append drop(_I** _est* o.**) ///
			addtext(Controls, YES, Firm FE, YES, Y-var, `ivar', Cluster, Firm)

		// ---------- 回归 3 ----------
		reg `ivar' `FRM_10' $control
		outreg2 using "`outreg_file0'", ///
			stat(coef se) bdec(3) sdec(3) addstat(Adj. R-squared, e(r2_a)) ///
			title("Lag `lag' - FRM_10") append drop(_I** _est* o.**) ///
			addtext(Controls, YES, Firm FE, YES, Y-var, `ivar', Cluster, Firm)
	}
}

local varlist "QrVaR_1 QrVaR_5 QrVaR_10 QrVaR_99 QrVaR_95 QrVaR_90"


foreach ivar of local varlist {
	replace `ivar' = `ivar'*100
}
	
//local varlist "Return"

// 外层循环：遍历所有风险度量变量
foreach ivar of local varlist {
	
	
	// 内层循环：遍历所有滞后期
	foreach lag in 5 21 63 84 105 126 {

		// 定义不同滞后期的 FRM 指标
		local FRM_1  L`lag'_FRM_1_High_Low
		local FRM_5  L`lag'_FRM_5_High_Low
		local FRM_10 L`lag'_FRM_10_High_Low

		// 定义输出路径（确保 $save/$date 已定义）
		local outreg_file0 "$date/`ivar'_1015_`lag'.xls"

		// ---------- 回归 1 ----------
		reg `ivar' `FRM_1'
		outreg2 using "`outreg_file0'", ///
			stat(coef se) bdec(3) sdec(3) addstat(Adj. R-squared, e(r2_a)) ///
			title("Lag `lag' - FRM_1") replace drop(_I** _est* o.**) ///
			addtext(Controls, YES, Firm FE, YES, Y-var, `ivar', Cluster, Firm)

		// ---------- 回归 2 ----------
		reg `ivar' `FRM_5'
		outreg2 using "`outreg_file0'", ///
			stat(coef se) bdec(3) sdec(3) addstat(Adj. R-squared, e(r2_a)) ///
			title("Lag `lag' - FRM_5") append drop(_I** _est* o.**) ///
			addtext(Controls, YES, Firm FE, YES, Y-var, `ivar', Cluster, Firm)

		// ---------- 回归 3 ----------
		reg `ivar' `FRM_10' 
		outreg2 using "`outreg_file0'", ///
			stat(coef se) bdec(3) sdec(3) addstat(Adj. R-squared, e(r2_a)) ///
			title("Lag `lag' - FRM_10") append drop(_I** _est* o.**) ///
			addtext(Controls, YES, Firm FE, YES, Y-var, `ivar', Cluster, Firm)
	}
}


*------------------------------------------------------------*
* 1. 数据准备
*------------------------------------------------------------*
use InsampleBeta_Premia_1014_month.dta, clear
tsset Date

/********************
  2. 生成滞后变量
*********************/
local laglist 5 21 63 84 105 126
foreach lag of local laglist {
    gen L`lag'_FRM_5_High_Low = L`lag'.FRM_5_High_Low
    gen L`lag'_FRM_1_High_Low = L`lag'.FRM_1_High_Low
    gen L`lag'_FRM_10_High_Low = L`lag'.FRM_10_High_Low
}






global control "Change_TY3M Slope TED RealEstate_excess MKreturn MKvol CarbonVol_Shenzhen CarbonVol_Guangdong CarbonVol_Hubei"


local laglist 5 21 63 84 105 126
foreach lag of local laglist {
    gen L`lag'_FRM_5_High_Low = L`lag'.FRM_5_High_Low
    gen L`lag'_FRM_1_High_Low = L`lag'.FRM_1_High_Low
    gen L`lag'_FRM_10_High_Low = L`lag'.FRM_10_High_Low
	
}

local varlist L63_FRM_5_High_Low $control

gen dropflag = 0
foreach var of local varlist {
    replace dropflag = 1 if missing(`var')
}

drop if dropflag==1
drop dropflag


gen t = _n
tsset t


program define myreg, rclass
    regress Return L63_FRM_5_High_Low $control
    return scalar beta_CARDI = _b[L63_FRM_5_High_Low]
end

bootstrap beta=r(beta_CARDI), reps(1000) seed(12345): myreg


* Newey–West HAC 标准误（滞后 12）
newey Return L84_FRM_5_High_Low , lag(9)

newey Return L63_FRM_5_High_Low  , lag(12)




newey Return L63_FRM_5_High_Low $control, lag(12)


reg Return L63_FRM_5_High_Low $control

gen Date_td = date(Date, "DMY")  

* 2. 设置日期显示格式
format Date_td %td  

* 3. tsset
tsset Date

* 补齐缺失日期
reg Return L63_FRM_5_High_Low $control, vce(hac 12)

tsset Date, daily





local varlist "Return"
// 外层循环：遍历所有风险度量变量
foreach ivar of local varlist {
	
	
	// 内层循环：遍历所有滞后期
	foreach lag in 5 21 63 84 105 126 {

		// 定义不同滞后期的 FRM 指标
		local FRM_1  L`lag'_FRM_1_High_Low
		local FRM_5  L`lag'_FRM_5_High_Low
		local FRM_10 L`lag'_FRM_10_High_Low

		// 定义输出路径（确保 $save/$date 已定义）
		local outreg_file0 "$date/`ivar'_1015_`lag'_month.xls"

		// ---------- 回归 1 ----------
		reg `ivar' `FRM_1' $control
		outreg2 using "`outreg_file0'", ///
			stat(coef se) bdec(3) sdec(3) addstat(Adj. R-squared, e(r2_a)) ///
			title("Lag `lag' - FRM_1") replace drop(_I** _est* o.**) ///
			addtext(Controls, YES, Firm FE, YES, Y-var, `ivar', Cluster, Firm)

		// ---------- 回归 2 ----------
		reg `ivar' `FRM_5' $control
		outreg2 using "`outreg_file0'", ///
			stat(coef se) bdec(3) sdec(3) addstat(Adj. R-squared, e(r2_a)) ///
			title("Lag `lag' - FRM_5") append drop(_I** _est* o.**) ///
			addtext(Controls, YES, Firm FE, YES, Y-var, `ivar', Cluster, Firm)

		// ---------- 回归 3 ----------
		reg `ivar' `FRM_10' $control
		outreg2 using "`outreg_file0'", ///
			stat(coef se) bdec(3) sdec(3) addstat(Adj. R-squared, e(r2_a)) ///
			title("Lag `lag' - FRM_10") append drop(_I** _est* o.**) ///
			addtext(Controls, YES, Firm FE, YES, Y-var, `ivar', Cluster, Firm)
	}
}
	
	