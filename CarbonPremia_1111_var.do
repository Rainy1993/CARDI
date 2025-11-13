

clear
clear matrix
clear mata
global date = "20251015"
global root = "/Users/ruting/Documents/macbook/PcBack/23.Topic_wkh_cooperation/FRM/Output/EfficiencyTest/"
cd $root

set maxvar  120000

/**************** 1. data generating ****************/
use InsampleBeta_Premia_1014_daily.dta, clear
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

// 定义变量组
local varlist "VaRRMAY_95 VaREMAY_95 VaREMAY_99 VaRRMAY_99 VaRRMAY_90 VaREMAY_90 VaRRMAY_5 VaREMAY_5 VaRRMAY_1 VaREMAY_1 VaRRMAY_10 VaREMAY_10 CResidual"

foreach ivar of local varlist {
	replace `ivar' = `ivar'*100
}

save InsampleBeta_Premia_1014_daily_lag.dta, replace

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
    quietly regress VaREMAY_5 L`lag'_FRM_`frm'_High_Low $control
    return scalar beta_CARDI = _b[L`lag'_FRM_`frm'_High_Low]
end

* Step 3. 双循环
foreach lag in 63 84 105 {
    foreach frm in 5 {
        use InsampleBeta_Premia_1014_daily_lag.dta, clear
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
	    local filename = "boot_lag`lag'_frm`frm'_var.dta"
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
export excel using "Bootstrap_FRM_var_results.xlsx", replace firstrow(variables)



