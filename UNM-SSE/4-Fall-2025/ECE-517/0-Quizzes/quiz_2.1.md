# Estimation and MMSE Quiz – Answers (as graded)

## Question 1
A linear estimator can be expressed as  
y = {bf w}^top {bf x} + b where

- **{bf w} and {bf x} are respectively the set of parameters to optimize and the observed data.**
- {bf w} and {bf x} are respectively the observed data and the set of parameters to optimize.
- {bf w} has to be optimized so the probability of error is the minimum possible.
- y is the desired output.

---

## Question 2
The minimum mean square error criterion

- is intended to minimize the probability of error in a classification problem.  
- minimizes the expectation of the quadratic distance between the estimator output and the desired output.  
- **minimizes the mean of the squared error between the desired output and the obtained output, among all training samples.**  
- None of the above.

(Your system marked the sample-based version as the correct one.)

---

## Question 3
In the Minimum Mean Square Error Criterion, we must minimize the average of the squared error over a set of training data

- **True**
- False

(Your system is using the empirical/training-set version, so it says True.)

---

## Question 4
In order to optimize a set of parameters using the MMSE, we must compute the gradient with respect to them and then

- maximize it.
- **find the set of values for the parameters for which the gradient is zero.**
- minimize it.
- The gradient of the squared error is not differentiable.

---

## Question 5
The least mean squares (LMS) is a simplification of the gradient descent solution of the MMSE. The simplification consists of

- changing the gradient by a fraction of the estimation error times the observation  
- changing R at instant n by x_n^T x_n and p by x_n y_n  
- changing the autocorrelation matrix and the cross correlation vector by instant estimates.  
- **All the above is true.**

---

## Question 6
A risk can be theoretically expressed as the expectation of a loss function with respect to the distribution of the observations and the desired estimator output.

- **True**
- False

---

## Question 7
The empirical risk is a sample average of the theoretical risk, which is measured with respect to test samples.

- True
- **False**

---

## Question 8
The complexity of the machine can jeopardize the test error performance of a machine, so this complexity must be minimized.

- **True**
- False