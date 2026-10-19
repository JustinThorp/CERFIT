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

# CERFIT 0.2.1

* Added URL and BugReports fields to DESCRIPTION
* Fixed an off-by-one error in the default gridval that dropped the maximum treatment value
* Added guard to reject continuous treatment with a gridval of 10 or fewer values 
* Fixed `randomForest` method IPTW bug
* Updated gridval documentation
* Added unit tests
* ITE prediction by iter now works with more than two arms and always returns
  a named list with one matrix per contrast against the reference arm (see
  `?predict.CERFIT`). For binary treatments this is a list of length one
  rather than a matrix, so existing code should use `ite[[1]]`.
* Added baseline model offset to prediction by iter when using residuals 
* Prediction by iter names list elements by gridval value for continuous treatments
* Added documentation for reproducibility with continuous treatments
* Predict function errors when trying to predict optimal treatment by iter instead of returning null
