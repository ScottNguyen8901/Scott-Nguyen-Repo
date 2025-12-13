# Here you can import the data using the code provided by UCI.
# As it can be seen, this is the most reliable and easy way to import the data.
# You can also use scikit learn or pandas.


from ucimlrepo import fetch_ucirepo

# fetch dataset
breast_cancer = fetch_ucirepo(id=14) # If the repository has standarized CSV file for pandas to parse, just use the data ID number to download. The ID is in the dataset url. FOr example: https://archive.ics.uci.edu/dataset/53/iris

# data (as pandas dataframes)
X = breast_cancer.data.features
y = breast_cancer.data.targets

# metadata
print(breast_cancer.metadata)

# variable information
print(breast_cancer.variables)

from ucimlrepo import fetch_ucirepo

# fetch dataset
adult = fetch_ucirepo(id=53)

# data (as pandas dataframes)
X = adult.data.features
y = adult.data.targets

# metadata
print(adult.metadata)

# variable information
print(adult.variables)

from ucimlrepo import fetch_ucirepo

# fetch dataset
wine = fetch_ucirepo(id=109)

# data (as pandas dataframes)
X = wine.data.features
y = wine.data.targets

# metadata
print(wine.metadata)

# variable information
print(wine.variables)

import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs
from sklearn.svm import LinearSVC

# Generate synthetic data with 4 classes
X, y = make_blobs(n_samples=100, centers=4, random_state=42, cluster_std=1.5)

# Train one-vs-all linear SVM classifiers
n_classes = len(np.unique(y))
svm_classifiers = []
for i in range(n_classes):
    # Create a binary target where class i is positive and others are negative
    y_binary = (y == i).astype(int)
    svm = LinearSVC(random_state=42)
    svm.fit(X, y_binary)
    svm_classifiers.append(svm)

# Plot the data and decision boundaries
plt.figure(figsize=(10, 8))
scatter = plt.scatter(X[:, 0], X[:, 1], c=y, cmap='viridis', s=50)

# Plot decision boundaries
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.02),
                     np.arange(y_min, y_max, 0.02))

# Predict the class for each point in the meshgrid
Z = np.argmax([clf.decision_function(np.c_[xx.ravel(), yy.ravel()]) for clf in svm_classifiers], axis=0)
Z = Z.reshape(xx.shape)

plt.contourf(xx, yy, Z, alpha=0.3, cmap='viridis')
plt.title('One-vs-All Linear SVM Classification')
plt.xlabel('Feature 1')
plt.ylabel('Feature 2')
plt.colorbar(scatter, label='Class')
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs
from sklearn.svm import LinearSVC

# Generate synthetic data with 4 classes (using the same data as before for comparison)
X, y = make_blobs(n_samples=100, centers=4, random_state=42, cluster_std=1.5)

# Train one-vs-one linear SVM classifiers
n_classes = len(np.unique(y))
svm_classifiers_ovo = []
for i in range(n_classes):
    for j in range(i + 1, n_classes):
        # Select data for classes i and j
        X_ij = X[(y == i) | (y == j)]
        y_ij = y[(y == i) | (y == j)]

        # Create a binary target where class i is positive (1) and class j is negative (0)
        y_binary_ij = (y_ij == i).astype(int)

        svm = LinearSVC(random_state=42)
        svm.fit(X_ij, y_binary_ij)
        svm_classifiers_ovo.append((i, j, svm))

# Plot the data and decision boundaries
plt.figure(figsize=(10, 8))
scatter = plt.scatter(X[:, 0], X[:, 1], c=y, cmap='viridis', s=50)

# Plot decision boundaries
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.02),
                     np.arange(y_min, y_max, 0.02))

# Predict the class for each point in the meshgrid using one-vs-one
Z_ovo = np.zeros(xx.shape, dtype=int)
for i in range(xx.shape[0]):
    for j in range(xx.shape[1]):
        x_point = np.array([[xx[i, j], yy[i, j]]])
        votes = np.zeros(n_classes)
        for class1, class2, svm in svm_classifiers_ovo:
            decision = svm.decision_function(x_point)
            if decision > 0:
                votes[class1] += 1
            else:
                votes[class2] += 1
        Z_ovo[i, j] = np.argmax(votes)

plt.contourf(xx, yy, Z_ovo, alpha=0.3, cmap='viridis')
plt.title('One-vs-One Linear SVM Classification with Synthetic Data')
plt.xlabel('Feature 1')
plt.ylabel('Feature 2')
plt.colorbar(scatter, label='Class')
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs
from sklearn.svm import SVC # Import SVC

# Generate synthetic data with 4 classes (using the same data as before for comparison)
X, y = make_blobs(n_samples=300, centers=4, random_state=42, cluster_std=1.5)

# Train a standard multiclass SVM classifier using SVC with an RBF kernel
svm_multiclass_rbf = SVC(kernel='rbf', random_state=42) # Use SVC with RBF kernel
svm_multiclass_rbf.fit(X, y)

# Plot the data and decision boundaries
plt.figure(figsize=(10, 8))
scatter = plt.scatter(X[:, 0], X[:, 1], c=y, cmap='viridis', s=50)

# Plot decision boundaries
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.02),
                     np.arange(y_min, y_max, 0.02))

# Predict the class for each point in the meshgrid using the multiclass SVM with RBF kernel
Z_multiclass_rbf = svm_multiclass_rbf.predict(np.c_[xx.ravel(), yy.ravel()])
Z_multiclass_rbf = Z_multiclass_rbf.reshape(xx.shape)

plt.contourf(xx, yy, Z_multiclass_rbf, alpha=0.3, cmap='viridis')
plt.title('Standard Multiclass SVM Classification with RBF Kernel') # Updated title
plt.xlabel('Feature 1')
plt.ylabel('Feature 2')
plt.colorbar(scatter, label='Class')
plt.show()

