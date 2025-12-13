# Introduction to RKHS – Answers

## Question 1  
Why do we need nonlinear extensions to the linear algorithms  
Group of answer choices

- **To be able to learn the nonlinear behaviour of observable phenomena.**  
- To provide the user wit an alternative formulation.  
- To make classification algorithms faster.  
- To make regression algorithms more powerful. 

---

## Question 2  
What is a general approach top provide an algorithm with nonlinear properties?  
Group of answer choices

- Transforming the input pattern ${\bf x}$ using a nonlinear transformation.  
- **Transforming the input pattern ${\bf x}$ into a higher dimensional space using a nonlinear transformation of its components.**  
- **Using a Volterra expansion of the data.**  
- Applying a polynomial to the input pattern $\bf x$

---

## Question 3  
What is the curse of dimensionality?  
Group of answer choices

- **If we want to add nonlinear properties to the data using a nonlinear transformation of the input space, the corresponding feature space has a dimension that increases with the degree of nonlinearity.**  
- If the dimension of the input and the output of our learning machine are too high, then the machine has high complexity.  
- The curse of dimensionality does not really exist unless the transformation is polynomial with a high order, in which case, the dimensionality of the corresponding space increases polynomially.  
- The curse of dimensionality is an unavoidable increase of the complexity of nonlinear learning machines.

---

## Question 4  
Which one is an interpretation of the Representer theorem.  
Group of answer choices

- **All answers are correct.**  
- Under some conditions, the estimator is a linear combination of dot products between the training and the test data.  
- If the criterion for the optimization of the data contains a convex loss and a nonnegative function of the weight vector norm, then the weight vector norm is a linear function of the training data.  
- The estimator can be constructed as a dot product between the weight vector and the test data into the Hilbert space, or equivalently as a linear combination of dot products between training and test data in the Hilbert space, if certain conditions are satisfied.

---

## Question 5  
The curse of dimensionality can be avoided  
Group of answer choices

- if the estimator is constructed under the conditions of the Representer Theorem. In that case, the estimator can be constructed in a dual space whose dimension is constant.  
- **if the estimator is constructed under the conditions of the Representer Theorem. In that case, the estimator can be constructed in a dual space whose dimension is equal to or less than the total number of training data.**  
- if the estimator is constructed under the conditions of the Representer Theorem. In that case, the curse of dimensionality does not exist.  
- if the estimator is constructed under the conditions of the Representer Theorem, which consists of using a nonlinear representation that transforms the data into a dual space, and a convex cost function.
