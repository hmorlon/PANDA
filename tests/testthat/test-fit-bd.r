# Load data and get tot_time
data(Cetacea)
tot_time<-max(node.age(Cetacea)$ages)
options(digits=16)

# Parameters for validation
precision_lh <- 1e-3
precision_aicc <- 5e-2


context("B constant")
# B constant
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){0}
lamb_par <- c(0.09)
mu_par <- c()
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, cst.lamb = TRUE, fix.mu=TRUE, cond="stem")

reference_lh <- -279.0280
reference_aicc <- 560.1031
# reference_aicc <- 560.0796
test_that("B constant",{
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
})

# BD constant
context("B & D constant")
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){y[1]}
lamb_par <- c(0.09)
mu_par <- c(0.3)
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, cst.lamb=TRUE, cst.mu=TRUE, cond="stem")
reference_lh <- -279.0280
reference_aicc <- 562.1989
# reference_aicc <- 562.1270

test_that("BD constant",{
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
})


# 3) B variable E
context("B exponential")
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.05, 0.01)
mu_par<-c()
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, expo.lamb=TRUE, fix.mu=TRUE, cond="stem")
reference_lh <- -278.9887
reference_aicc <- 562.0485
# reference_aicc <- 562.1213

test_that("B variable exponential",{
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
})

# 4) B variable L
## TODO fix numerical results
# context("B linear")
# f.lamb <-function(t,y){y[1] + y[2] * t}
# f.mu<-function(t,y){0}
# lamb_par<-c(0.05, 0.01)
# mu_par<-c()
# res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, fix.mu=TRUE, cond="stem")
# reference_lh <- -278.9896
# reference_aicc <- 562.0502
# # reference_aicc <- 562.0485
#
# test_that("B variable linear",{
#   expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
#   expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
# })


# 5) B variable E, D constant
context("B exponential, D constant")
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
mu_par <- c(0.5)
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.05, 0.01)
mu_par<-c(0.1)
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, expo.lamb=TRUE, cst.mu=TRUE, cond="stem")
reference_lh <- -278.9887
reference_aicc <- 564.1204
# reference_aicc <- 564.2676

test_that("B variable exponential, D constant",{
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
})


# 6) B variable L, D constant
## TODO fix numerical results
# context("B linear, D constant")
# f.lamb <-function(t,y){y[1] + y[2] * t}
# f.mu<-function(t,y){y[1]}
# lamb_par<-c(0.5, 0.01)
# mu_par <- c(0.5)
# res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, cst.mu= TRUE, cond="stem")
# reference_lh <- -278.9896
# reference_aicc <- 564.1221
# # reference_aicc <- 564.2676
#
# test_that("B variable linear, D constant",{
#   expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
#   expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
# })


# 7) B constant, D variable E
context("B constant, D exponential")
f.lamb<-function(t,y){y[1]}
f.mu <-function(t,y){y[1] * exp(y[2] * t)}
lamb_par <- c(0.5)
mu_par <-c(0.5, 0.01)
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, cst.lamb=TRUE, expo.mu=TRUE, cond="stem")
reference_lh <- -279.0280
reference_aicc <- 564.1989
# reference_aicc <- 564.2676

test_that("B constant, D exponential",{
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
})


# 8) B constant, D variable L
## TODO fix numerical results
# context("B constant, D linear")
# f.lamb<-function(t,y){y[1]}
# f.mu <-function(t,y){y[1] + y[2] * t}
# lamb_par <- c(0.5)
# mu_par <-c(0.5, 0.01)
# res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, cst.lamb=TRUE, cond="stem")
# reference_lh <- -279.0280
# reference_aicc <- 564.1989
# # reference_aicc <- 564.2676
#
# test_that("B constant, D linear",{
#   expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
#   expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
# })
