# CERFIT 0.1.1

* Added Out-of-Bag (OOB) predictions
* Fixed Spelling Errors

# CERFIT 0.2.0

* Added survival analysis support
* Added Cox proportional hazards residuals (martingale and deviance)
* Added support for time-dependent covariates in counting-process format
* Added response parameter for explicit response type specification
* Added surv_resid parameter to select residual type for survival models
* Added automatic type coercion in predict() to align newdata column classes with training data
* Fixed treatment variable extraction to use vector instead of one-column data frame
* Fixed continuous treatment detection to require is.numeric() check
* Fixed continuous ITE column naming to use gridval instead of utrt
* Fixed dimension dropping in single-contrast ITE predictions with drop = FALSE
* Fixed spelling errors
* Changed unique-treatment grid threshold from 20 to 15
* Expanded documentation