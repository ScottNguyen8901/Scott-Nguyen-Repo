import os
import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    roc_curve,
    auc,
    confusion_matrix,
    f1_score,
    accuracy_score,
    precision_score,
    recall_score
)

# ---------------------------------------------------------
# Problem 5: Multiobjective Optimization
# ---------------------------------------------------------

os.makedirs("plots", exist_ok=True)

# load dataset
digits = load_digits()
X = digits.data
y = digits.target

# binary labels
# positive (True): 2, 3, 5, 7
# negative (False): 0, 1, 4, 6, 8, 9
y_binary = np.isin(y, [2, 3, 5, 7]).astype(int)

# train/test split
X_train, X_test, y_train_bin, y_test_bin = train_test_split(
    X,
    y_binary,
    test_size=0.2,
    random_state=42,
    stratify=y_binary
)

# train classifier
rf = RandomForestClassifier(
    n_estimators=500,
    max_depth=None,
    random_state=42
)
rf.fit(X_train, y_train_bin)

# prediction scores
y_score = rf.predict_proba(X_test)[:, 1]

# ROC curve
fpr, tpr, thresholds = roc_curve(y_test_bin, y_score)
roc_auc = auc(fpr, tpr)

# ---------------------------------------------------------
# helper functions
# ---------------------------------------------------------
def compute_single_point_metrics(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

    acc = accuracy_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

    return {
        "Confusion Matrix": np.array([[tn, fp], [fn, tp]]),
        "Accuracy": acc,
        "F1 Score": f1,
        "Precision": precision,
        "Recall": recall,
        "Sensitivity": sensitivity,
        "Specificity": specificity
    }

def print_solution(name, idx):
    th = thresholds[idx]
    chosen_fpr = fpr[idx]
    chosen_tpr = tpr[idx]

    y_pred = (y_score >= th).astype(int)
    metrics = compute_single_point_metrics(y_test_bin, y_pred)

    print("==================================================")
    print(name)
    print("threshold:", th)
    print("fpr:", chosen_fpr)
    print("tpr:", chosen_tpr)
    for key, value in metrics.items():
        print(f"{key}:")
        print(value)
    print()

# ---------------------------------------------------------
# (i) maximize tpr subject to fpr <= 0.1
# ---------------------------------------------------------
valid_i = np.where(fpr <= 0.1)[0]
idx_i = valid_i[np.argmax(tpr[valid_i])]

# ---------------------------------------------------------
# (ii) minimize fpr subject to tpr >= 0.9
# ---------------------------------------------------------
valid_ii = np.where(tpr >= 0.9)[0]
idx_ii = valid_ii[np.argmin(fpr[valid_ii])]

# ---------------------------------------------------------
# (iii) balanced approach
# minimize (-tpr + 2*fpr)
# subject to (tpr >= 0.8) and (fpr <= 0.1)
# ---------------------------------------------------------
valid_iii = np.where((tpr >= 0.8) & (fpr <= 0.1))[0]
objective_iii = -tpr[valid_iii] + 2.0 * fpr[valid_iii]
idx_iii = valid_iii[np.argmin(objective_iii)]

# ---------------------------------------------------------
# print results
# ---------------------------------------------------------
print("ROC AUC:", roc_auc)
print()

print_solution("(i) Maximize tpr subject to fpr <= 0.1", idx_i)
print_solution("(ii) Minimize fpr subject to tpr >= 0.9", idx_ii)
print_solution("(iii) Balanced approach: minimize (-tpr + 2*fpr)", idx_iii)

# ---------------------------------------------------------
# (iv) explanation of negative sign
# ---------------------------------------------------------
print("==================================================")
print("(iv) Why is there a negative sign in front of tpr?")
print("Because we want to maximize tpr, but the optimization is written as a minimization problem.")
print("Using -tpr converts maximizing tpr into minimizing -tpr.")
print()

# ---------------------------------------------------------
# (v) slope relationship
# For objective J = -tpr + 2*fpr
# dJ/d(threshold) = -dtpr/dtau + 2*dfpr/dtau = 0
# so dtpr/dfpr = 2
# ---------------------------------------------------------
print("==================================================")
print("(v) Slope relationship for the balanced approach")
print("For J = -tpr + 2*fpr, setting dJ/d(threshold) = 0 gives:")
print("    -d(tpr)/d(threshold) + 2*d(fpr)/d(threshold) = 0")
print("So:")
print("    d(tpr)/d(fpr) = 2")
print("Thus, the optimal point is related to a ROC slope of 2.")
print()

# ---------------------------------------------------------
# plot ROC curve and selected operating points
# ---------------------------------------------------------
plt.figure()
plt.plot(fpr, tpr, label=f"ROC curve (AUC = {roc_auc:.4f})")
plt.scatter(fpr[idx_i], tpr[idx_i], label="(i) max tpr, fpr<=0.1")
plt.scatter(fpr[idx_ii], tpr[idx_ii], label="(ii) min fpr, tpr>=0.9")
plt.scatter(fpr[idx_iii], tpr[idx_iii], label="(iii) balanced")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Problem 5 ROC Operating Points")
plt.legend()
plt.grid()
plt.savefig("plots/problem5_roc.png")
plt.close()