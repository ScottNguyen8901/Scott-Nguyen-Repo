# Machine Learning Quiz 2.2 – Answers

## Question 1
Mark the WRONG sentence

- The VC dimension of a linear machine is equal to the dimension of the space plus one.  
- **The VC dimension only works in linear machines.**  
- The VC dimension is a measure of the complexity of a classification machine.  
- If the VC dimension of a machine is higher than the number of training data, the machine will probably overfit.

---

## Question 2
The VC dimension is a way to determine the maximum complexity of a learning machine.

- **True**
- False

---

## Question 3
If the VC dimension of a machine is higher than the number of samples in the training set, then the machine is guaranteed to overfit:

- **if we do not limit the complexity, because the machine will be able to shatter all the points.**  
- because it will always correctly classify all points regardless of the control for the complexity.  
- Not true. The higher the VC dimension is the better the machine will work.

---

## Question 4
The empirical risk is

- **the expectation of the error over the training samples.**  
- the sample average of the error over the training samples.  
- the sample average of the error over the test samples.  
- the expectation of the error over the test samples.

---

## Question 5
The VC dimension is restricted to linear machines.

- True  
- **False**

---

## Question 6
The actual risk is the expectation of the error over the test samples, and it is bounded by the empirical risk plus the structural risk, where

- the structural risk is a term that does not depend on anything related to the data, but only on the VC dimension.  
- the structural risk is a term that decreases when the VC dimension increases, and it is only true with probability 1 - η.  
- the structural risk increases with the VC dimension and the number of data.  
- **the structural risk is a term that increases with the VC dimension of the classifier, decreases with the number of training data. The actual risk is bounded with a probability less than 1.**

---

## Question 7
In Support Vector Machines, the complexity of the machine is limited by minimizing h, which is equivalent to minimize ||w||².

- **True**
- False

---

## Question 8
In SVM the minimization of the complexity is done by minimizing h, which is equivalent to minimize the distance d of the margin.

- True  
- **False**
