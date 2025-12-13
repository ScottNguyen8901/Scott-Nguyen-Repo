# Machine Learning Quiz 3.1 – Answers

## Question 1
The SVM criterion takes into account the structural risk and the empirical risk. The optimization is intended to minimize both at the same time.

- **True**  
- False

---

## Question 2
The SVM optimization is performed through a Lagrange analysis. It leads to

- a dual functional that depends on the corresponding Lagrange multipliers. The solution is always guaranteed but it is not unique.  
- **The solution of the dual is always guaranteed because the matrix of the functional is non negative definite.**  
- The solution is guaranteed because the functional is the sum of a convex and a concave functions.  
- The solution is never guaranteed.

---

## Question 3
In a SVM classifier, the weight vector is a linear combination of the training data.

- **True**  
- False

---

## Question 4
In a Support Vector Machine, the solution depends on a subset of the training data.

- **True**  
- False

---

## Question 5
The KKT conditions

- are obtained from nulling the gradient of the lagrange functional with respect to the primal variables (except the complementary ones).  
- **The KKT conditions are obtained by nulling the gradient of the Lagrange functional wrt the primal and dual variables (except the complementary ones).**   
- **The complementary KKT conditions say that either the constraints or the corresponding Lagrange multiplier is zero.**  
- The KKT assure that the constraints vanish from the Lagrangian at the optimal point.

---

## Question 6
The dual variables or Lagrange multipliers

- **of the samples that are properly classified and outside the origin are zero.**  
- of the samples that are properly classified and inside the margin have a value strictly less than C.  
- **of the samples that are inside the margin have a value equal to C.**  
- of the samples that are misclassified and outside the margin are C.

---

## Question 7
All the dual variables are either 0 or C.

- True  
- **False**
