# Regression & SVR Quiz – Answers

## Question 1  
Regression can be defined as the estimation of one or more discrete variables from a set of observations.

- True  
- **False**

---

## Question 2  
In order to perform regression, we always need to hide the bias in the weight vector **w**.

- True  
- **False**

---

## Question 3  
The MMSE criterion

- Is intended to minimize the expectation of the error, but since the distribution of the observed data is unknown, we need to approximate it by the minimization of a sample average of the quadratic error.  
- Is an optimization criterion that always has a solution and it is unique.  
- Is an optimization criterion that can be achieved by several algorithms in block or iteratively.  
- **All of them are true.**

---

## Question 4  
The Least Mean Squares (LMS) algorithm

- Is an approximation of the MMSE.  
- Is an algorithm that approximates an optimal criterion by gradient descent.  
- Is an algorithm that approximates the square-error gradient by the product of the error times the observation.  
- **All of the above is right.**

---

## Question 5  
Support Vector Regression

- **considers separately positive and negative errors in order to reduce the solution bias.**  
- uses positive slack variables because only the positive errors are taken into account.  
- **uses positive slack variables that represent the absolute value of the error minus a small quantity which is taken as an error tolerance.**  
- **is an optimization criterion that minimizes a structural term represented by the norm of the weights plus an empirical term that is a linear function of the error.**

---

## Question 6  
The dual variables of the support vectors can be expressed as \( \alpha_n - \alpha_n^\* \)

- **where \( \alpha_n \) is positive if the error is positive and larger than ε (with \( \alpha_n^\* = 0 \)), and \( \alpha_n^\* \) is positive if the error is negative and lower than −ε (with \( \alpha_n = 0 \)).**  
- If the sample is outside the margin, the corresponding dual variable is less than C, and if the sample is inside the margin, the dual variable is zero.  
- **If the sample has an error which is less th**
