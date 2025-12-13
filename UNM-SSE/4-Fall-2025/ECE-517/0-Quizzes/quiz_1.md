# Machine Learning Quiz 1 – Answers

## Question 1
**How do machines learn?**
- They need to be feed with human-produced knowledge.  
- They learn using deductive principles.  
- Machines cannot learn, only humans can.  
- **They use inductive learning from data.**

---

## Question 2
**Machine learning can be seen as the process of:**
- Taking data, extracting the knowledge and produce information.  
- **Taking data, extracting the information and inferring knowledge.**  
- Learning heuristic rules in order to take decisions.  
- None of the previous is valid.  

---

## Question 3
**Which one of these sentences is true?**
- Data needs to be fed to the algorithm in form of numbers.  
- The data fed into the machine consist of vectors of features.  
- Each vector used for training belongs to a sample of data.  
- **All of these are true.**

---

## Question 4
**A linear classifier consists of:**
- **One or several hyperplanes that separate the patterns into different classes.**  
- Just one hyperplane that separates the patterns in two classes.  
- Linear classifiers are useless in practice.  
- Linear classifiers work only in spaces of very low dimension.  

---

## Question 5
**A linear binary classification:**
- **is obtained by observing the sign of the estimator output.**  
- is obtained by computing a linear function over the observations.  
- has as many parameters as observed features, plus a scalar bias.  
- All of the above is true.

Note: the last option is not the best answer here because not every linear binary classifier is described only as "features + bias" in the same way (it depends on the model).

======================================================================

# Estimation and MMSE Quiz – Answers

## Question 1
**A linear estimator can be expressed as**  
\(y = \mathbf{\hat{w}}^{\text{top}} \mathbf{x} + b\) **where**
- **\(\mathbf{w}\) and \(\mathbf{x}\) are respectively the set of parameters to optimize and the observed data.**  
- \(\mathbf{w}\) and \(\mathbf{x}\) are respectively the observed data and the set of parameters to optimize.  
- \(\mathbf{w}\) has to be optimized so the probability of error is the minimum possible.  
- y is the desired output.

---

## Question 2
**The minimum mean square error (MMSE) criterion**
- is intended to minimize the probability of error in a classification problem.  
- **minimizes the expectation of the quadratic distance between the estimator output and the desired output.**  
- minimizes the mean of the squared error between the desired output and the obtained output, among all training samples.  
- None of the above.

---

## Question 3
**In the MMSE criterion, we must minimize the average of the squared error over a set of training data.**
- True  
- **False**

Explanation (short): MMSE is defined with respect to the expectation over the underlying data/target distribution; minimizing over only the training set is empirical MSE, not theoretical MMSE.

---

## Question 4
**In order to optimize a set of parameters using the MMSE, we must compute the gradient with respect to them and then**
- maximize it.  
- **find the set of values for the parameters for which the gradient is zero.**  
- minimize it.  
- The gradient of the squared error is not differentiable.

---

## Question 5
**The least mean squares (LMS) is a simplification of the gradient descent solution of the MMSE. The simplification consists of**
- changing the gradient by a fraction of the estimation error times the observation  
- changing R at instant n by x\_n^T x\_n and p by x\_n y\_n  
- changing the autocorrelation matrix and the cross correlation vector by instant estimates.  
- **All the above is true.**

---

## Question 6
**A risk can be theoretically expressed as the expectation of a loss function with respect to the distribution of the observations and the desired estimator output.**
- **True**  
- False

---

## Question 7
**The empirical risk is a sample average of the theoretical risk, which is measured with respect to test samples.**
- True  
- **False**

(Reason: empirical risk is computed over the available training samples, not the test set.)

---

## Question 8
**The complexity of the machine can jeopardize the test error performance of a machine, so this complexity must be minimized.**
- True  
- **False**

(Reason: we do not blindly minimize complexity; we control or regularize it to balance bias and variance.)
