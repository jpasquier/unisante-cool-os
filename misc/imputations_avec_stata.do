use "\\tsclient\jerome\coolos.dta", clear

mi set M = 20
mi set wide
mi register imputed L1L4s_BMD BS_Age PO_BMI DX_age DX_BMI Time_postop deltaWeight
mi register imputed femurTot_DMO col_DMO OneY_TWL OneY_EBMIL ALMI FMI VitD1 PTH P1NP Bcrosslaps glyc HOMA HbA1c
mi impute chained (pmm, knn(5)) femurTot_DMO col_DMO OneY_TWL OneY_EBMIL ALMI FMI VitD1 PTH P1NP Bcrosslaps glyc HOMA HbA1c L1L4s_BMD BS_Age PO_BMI DX_age DX_BMI Time_postop deltaWeight, replace rseed(666)

mi estimate: regress L1L4s_BMD FMI