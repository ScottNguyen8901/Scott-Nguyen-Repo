import matplotlib.pyplot as plt
import numpy as np
import os

from sklearn.datasets import load_digits

# create plots directory
os.makedirs("plots", exist_ok=True)

# ---------------- (1a iv) Ideal ROC ----------------
plt.figure()
fpr_ideal = [0, 0, 1]
tpr_ideal = [0, 1, 1]
plt.plot(fpr_ideal, tpr_ideal)
plt.scatter([0], [1])

plt.xlabel("False Positive Rate (1 - Specificity)")
plt.ylabel("True Positive Rate (Sensitivity)")
plt.title("Problem 1(a)(iv): Ideal ROC Curve")

plt.xlim(0, 1)
plt.ylim(0, 1)
plt.grid()

plt.savefig("plots/1a_iv_ideal_roc.png")
plt.close()


# ---------------- (1a v) Random Classifier ----------------
plt.figure()
fpr_rand = np.linspace(0, 1, 100)
tpr_rand = fpr_rand
plt.plot(fpr_rand, tpr_rand, linestyle='--')

plt.xlabel("False Positive Rate (1 - Specificity)")
plt.ylabel("True Positive Rate (Sensitivity)")
plt.title("Problem 1(a)(v): Random Classifier ROC")

plt.xlim(0, 1)
plt.ylim(0, 1)
plt.grid()

plt.savefig("plots/1a_v_random_roc.png")
plt.close()


# ---------------- (1a vi) Realistic ROC ----------------
plt.figure()
fpr_real = np.linspace(0, 1, 100)
tpr_real = 1 - (1 - fpr_real)**3
plt.plot(fpr_real, tpr_real)

# operating point (high sensitivity, small FPR)
plt.scatter([0.1], [1.0])

plt.xlabel("False Positive Rate (1 - Specificity)")
plt.ylabel("True Positive Rate (Sensitivity)")
plt.title("Problem 1(a)(vi): Realistic ROC with Operating Point")

plt.xlim(0, 1)
plt.ylim(0, 1)
plt.grid()

plt.savefig("plots/1a_vi_realistic_roc.png")
plt.close()

epochs = np.arange(0, 50)

# (i) Improving
plt.figure()
loss_improving = np.exp(-epochs/10)
plt.plot(epochs, loss_improving)
plt.title("Problem 2(a)(i): Improving Training Loss")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.grid()
plt.savefig("plots/2a_i_improving.png")
plt.close()

# (ii) Not improving
plt.figure()
loss_flat = np.ones_like(epochs) * 0.8
plt.plot(epochs, loss_flat)
plt.title("Problem 2(a)(ii): Not Improving")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.grid()
plt.savefig("plots/2a_ii_flat.png")
plt.close()

# (iii) Ideal (converged)
plt.figure()
loss_converged = np.exp(-epochs/10) + 0.1
plt.plot(epochs, loss_converged)
plt.title("Problem 2(a)(iii): Ideal Convergence")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.grid()
plt.savefig("plots/2a_iii_converged.png")
plt.close()

# ---------------- (i) Good Generalization ----------------
plt.figure()

train_loss = np.exp(-epochs/10)
val_loss = np.exp(-epochs/10) + 0.05  # slightly higher, close

plt.plot(epochs, train_loss, label="Training Loss")
plt.plot(epochs, val_loss, linestyle='--', label="Validation Loss")

plt.title("Problem 2(b)(i): Good Generalization")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2b_i_good_generalization.png")
plt.close()


# ---------------- (ii) Ideal Trend ----------------
plt.figure()

train_loss = np.exp(-epochs/10)
val_loss = 0.1 + 0.9 * np.exp(-epochs/8)  # converges slightly above train

plt.plot(epochs, train_loss, label="Training Loss")
plt.plot(epochs, val_loss, linestyle='--', label="Validation Loss")

plt.title("Problem 2(b)(ii): Ideal Trend")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2b_ii_ideal.png")
plt.close()


# ---------------- (iii) Worst Trend (Overfitting) ----------------
plt.figure()

train_loss = np.exp(-epochs/10)
val_loss = np.exp(-epochs/10) + 0.02 * epochs  # starts rising

plt.plot(epochs, train_loss, label="Training Loss")
plt.plot(epochs, val_loss, linestyle='--', label="Validation Loss")

plt.title("Problem 2(b)(iii): Overfitting (Worst Case)")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2b_iii_overfitting.png")
plt.close()

# ---------------- (i) Ideal ----------------
plt.figure()

train = np.exp(-epochs/10)
val = train + 0.02

plt.plot(epochs, train, label="Training Loss")
plt.plot(epochs, val, linestyle='--', label="Validation Loss")

plt.title("Problem 2(c)(i): Ideal")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2c_i_ideal.png")
plt.close()


# ---------------- (ii) Overfitting ----------------
plt.figure()

train = np.exp(-epochs/10)
val = np.exp(-epochs/10) + 0.03 * epochs

plt.plot(epochs, train, label="Training Loss")
plt.plot(epochs, val, linestyle='--', label="Validation Loss")

plt.title("Problem 2(c)(ii): Overfitting")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2c_ii_overfitting.png")
plt.close()


# ---------------- (iii) Underfitting ----------------
plt.figure()

train = 0.8 * np.ones_like(epochs)
val = 0.85 * np.ones_like(epochs)

plt.plot(epochs, train, label="Training Loss")
plt.plot(epochs, val, linestyle='--', label="Validation Loss")

plt.title("Problem 2(c)(iii): Underfitting")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2c_iii_underfitting.png")
plt.close()


# ---------------- (iv) Good Fit ----------------
plt.figure()

train = np.exp(-epochs/10)
val = train + 0.08

plt.plot(epochs, train, label="Training Loss")
plt.plot(epochs, val, linestyle='--', label="Validation Loss")

plt.title("Problem 2(c)(iv): Good Fit")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2c_iv_goodfit.png")
plt.close()

# ---------------- (i) No learning → Fix ----------------
plt.figure()

before = np.ones_like(epochs) * 0.9
after = np.exp(-epochs/10)

plt.plot(epochs, before, linestyle='--', label="Before (No Learning)")
plt.plot(epochs, after, label="After (Fixed)")

plt.title("Problem 2(d)(i): Fix No Learning")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2d_i_fix_no_learning.png")
plt.close()


# ---------------- (ii) Poor learning → Fix ----------------
plt.figure()

before = np.exp(-epochs/40)  # slow learning
after = np.exp(-epochs/10)   # improved learning

plt.plot(epochs, before, linestyle='--', label="Before (Poor Learning)")
plt.plot(epochs, after, label="After (Fixed)")

plt.title("Problem 2(d)(ii): Fix Poor Learning")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2d_ii_fix_poor_learning.png")
plt.close()


# ---------------- (iii) Underfitting → Fix ----------------
plt.figure()

before = 0.8 * np.ones_like(epochs)
after = np.exp(-epochs/10)

plt.plot(epochs, before, linestyle='--', label="Before (Underfitting)")
plt.plot(epochs, after, label="After (Fixed)")

plt.title("Problem 2(d)(iii): Fix Underfitting")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2d_iii_fix_underfitting.png")
plt.close()


# ---------------- (iv) Overfitting → Fix ----------------
plt.figure()

train = np.exp(-epochs/10)
val_before = np.exp(-epochs/10) + 0.03 * epochs
val_after = np.exp(-epochs/10) + 0.05

plt.plot(epochs, train, label="Training Loss")
plt.plot(epochs, val_before, linestyle='--', label="Before (Overfitting)")
plt.plot(epochs, val_after, linestyle=':', label="After (Fixed)")

plt.title("Problem 2(d)(iv): Fix Overfitting")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2d_iv_fix_overfitting.png")
plt.close()


# ---------------- (v) Improve good model ----------------
plt.figure()

before = np.exp(-epochs/10) + 0.1
after = np.exp(-epochs/10) + 0.05

plt.plot(epochs, before, linestyle='--', label="Before (Good)")
plt.plot(epochs, after, label="After (Improved)")

plt.title("Problem 2(d)(v): Improve Good Model")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.grid()

plt.savefig("plots/2d_v_improve_model.png")
plt.close()

np.random.seed(0)

# load digits
digits = load_digits()
X = digits.data
y = digits.target

# true: 2,3,5,7
true_mask = np.isin(y, [2, 3, 5, 7])
false_mask = np.logical_not(true_mask)

true_idx = np.where(true_mask)[0]
false_idx = np.where(false_mask)[0]

# bad balancing approach: sample with replacement BEFORE train/validation split
count_true = len(true_idx)
count_false = len(false_idx)

target_count = max(count_true, count_false)

if count_true < count_false:
    true_sample = np.random.choice(true_idx, target_count, replace=True)
    false_sample = np.random.choice(false_idx, target_count, replace=False)
else:
    true_sample = np.random.choice(true_idx, target_count, replace=False)
    false_sample = np.random.choice(false_idx, target_count, replace=True)

balanced_idx = np.concatenate([true_sample, false_sample])
np.random.shuffle(balanced_idx)

# split into training and validation
split = int(0.8 * len(balanced_idx))
train_idx = balanced_idx[:split]
val_idx = balanced_idx[split:]

# detect leakage: same original images appearing in both sets
overlap_idx = np.intersect1d(train_idx, val_idx)

print("number of training samples:", len(train_idx))
print("number of validation samples:", len(val_idx))
print("number of leaked original images:", len(overlap_idx))
print("example leaked indices:", overlap_idx[:20])

# optional: show repeated source images directly
print("leakage exists:", len(overlap_idx) > 0)

# load data
digits = load_digits()
X = digits.data
y = digits.target

# define classes
true_mask = np.isin(y, [2, 3, 5, 7])
false_mask = np.logical_not(true_mask)

true_idx = np.where(true_mask)[0]
false_idx = np.where(false_mask)[0]

# get counts
count_true = len(true_idx)
count_false = len(false_idx)

# downsample larger class to match smaller class
target_count = min(count_true, count_false)

true_sample = np.random.choice(true_idx, target_count, replace=False)
false_sample = np.random.choice(false_idx, target_count, replace=False)

# combine + shuffle
balanced_idx = np.concatenate([true_sample, false_sample])
np.random.shuffle(balanced_idx)

# -------- required outputs --------
print("first 10 elements:", balanced_idx[:10])
print("number of true samples:", len(true_sample))
print("number of false samples:", len(false_sample))

# ---------------- 3(c) Train/Test Split ----------------

# split 80/20
split_idx = int(0.8 * len(balanced_idx))

train_idx = balanced_idx[:split_idx]
test_idx = balanced_idx[split_idx:]

# get data + labels
X_train = X[train_idx]
y_train = y[train_idx]

X_test = X[test_idx]
y_test = y[test_idx]

# -------- required outputs --------
print("number of training samples:", len(train_idx))
print("number of testing samples:", len(test_idx))

print("X_train shape:", X_train.shape)
print("y_train shape:", y_train.shape)
print("X_test shape:", X_test.shape)
print("y_test shape:", y_test.shape)

import numpy as np
from sklearn.metrics import (
    roc_auc_score,
    f1_score,
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score
)
from sklearn.linear_model import LogisticRegression


def compute_metrics(y_true, y_pred, y_score):
    """
    Inputs:
        y_true  : ground-truth binary labels
        y_pred  : predicted binary labels at a chosen operating point
        y_score : predicted scores/probabilities for the positive class

    Returns:
        Dictionary containing:
            AUC
            F1
            Confusion Matrix
            Accuracy
            Precision
            Recall
            Sensitivity
            Specificity
    """

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

    auc = roc_auc_score(y_true, y_score)
    f1 = f1_score(y_true, y_pred)
    acc = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred)
    recall = recall_score(y_true, y_pred)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

    metrics = {
        "AUC": auc,
        "F1 Score": f1,
        "Confusion Matrix": np.array([[tn, fp], [fn, tp]]),
        "Accuracy": acc,
        "Precision": precision,
        "Recall": recall,
        "Sensitivity": sensitivity,
        "Specificity": specificity
    }

    return metrics


# ---------------- Binary labels for Problem 3 ----------------
# True digits: 2,3,5,7
# False digits: 0,1,4,6,8,9
y_binary = np.isin(y, [2, 3, 5, 7]).astype(int)

y_train_bin = y_binary[train_idx]
y_test_bin = y_binary[test_idx]

# ---------------- Train a classifier ----------------
model = LogisticRegression(max_iter=1000)
model.fit(X_train, y_train_bin)

# ---------------- Predictions ----------------
# scores/probabilities for positive class
y_score = model.predict_proba(X_test)[:, 1]

# predicted labels at a single operating point
y_pred = (y_score >= 0.5).astype(int)

# ---------------- Compute metrics ----------------
metrics = compute_metrics(y_test_bin, y_pred, y_score)

# ---------------- Print results ----------------
print("Inputs to compute_metrics:")
print("y_true shape:", y_test_bin.shape)
print("y_pred shape:", y_pred.shape)
print("y_score shape:", y_score.shape)
print()

for key, value in metrics.items():
    print(f"{key}:")
    print(value)
    print()

import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.base import clone
from sklearn.metrics import (
    roc_auc_score,
    f1_score,
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score
)

def compute_metrics(y_true, y_pred, y_score):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

    auc = roc_auc_score(y_true, y_score)
    f1 = f1_score(y_true, y_pred)
    acc = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred)
    recall = recall_score(y_true, y_pred)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

    return {
        "AUC": auc,
        "F1 Score": f1,
        "Confusion Matrix": np.array([[tn, fp], [fn, tp]]),
        "Accuracy": acc,
        "Precision": precision,
        "Recall": recall,
        "Sensitivity": sensitivity,
        "Specificity": specificity
    }


# ---------------------------------------------------------
# Problem 3(e): 5-fold cross-validation on balanced training set
# ---------------------------------------------------------

# Binary labels:
# True digits: 2, 3, 5, 7
# False digits: 0, 1, 4, 6, 8, 9
y_binary = np.isin(y, [2, 3, 5, 7]).astype(int)
y_train_bin = y_binary[train_idx]

# Random Forest classifier
rf_clf = RandomForestClassifier(
    n_estimators=100,
    random_state=42
)

# 5-fold stratified CV
skfolds = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

fold_results = []

for fold_num, (cv_train_index, cv_test_index) in enumerate(skfolds.split(X_train, y_train_bin), start=1):
    print(f"========== Fold {fold_num} ==========")

    # Split fold data
    X_train_folds = X_train[cv_train_index]
    y_train_folds = y_train_bin[cv_train_index]

    X_test_fold = X_train[cv_test_index]
    y_test_fold = y_train_bin[cv_test_index]

    # Verify balance in each fold
    train_true_count = np.sum(y_train_folds == 1)
    train_false_count = np.sum(y_train_folds == 0)
    test_true_count = np.sum(y_test_fold == 1)
    test_false_count = np.sum(y_test_fold == 0)

    print("Training fold class counts:")
    print("  True :", train_true_count)
    print("  False:", train_false_count)

    print("Validation fold class counts:")
    print("  True :", test_true_count)
    print("  False:", test_false_count)

    # Train classifier
    clone_clf = clone(rf_clf)
    clone_clf.fit(X_train_folds, y_train_folds)

    # Predictions
    y_score_fold = clone_clf.predict_proba(X_test_fold)[:, 1]
    y_pred_fold = (y_score_fold >= 0.5).astype(int)

    # Metrics
    metrics = compute_metrics(y_test_fold, y_pred_fold, y_score_fold)

    print("Metrics:")
    for key, value in metrics.items():
        print(f"{key}:")
        print(value)
    print()

    fold_results.append({
        "fold": fold_num,
        "train_true": train_true_count,
        "train_false": train_false_count,
        "test_true": test_true_count,
        "test_false": test_false_count,
        "metrics": metrics
    })

# ---------------------------------------------------------
# Highest AUC fold
# ---------------------------------------------------------
best_fold = max(fold_results, key=lambda r: r["metrics"]["AUC"])

print("========== Highest AUC Fold ==========")
print("Fold:", best_fold["fold"])
for key, value in best_fold["metrics"].items():
    print(f"{key}:")
    print(value)
print()

# ---------------------------------------------------------
# Lowest AUC fold
# ---------------------------------------------------------
worst_fold = min(fold_results, key=lambda r: r["metrics"]["AUC"])

print("========== Lowest AUC Fold ==========")
print("Fold:", worst_fold["fold"])
for key, value in worst_fold["metrics"].items():
    print(f"{key}:")
    print(value)
print()

# ---------------------------------------------------------
# Percentage variation between min and max values
# ---------------------------------------------------------
metric_names = [
    "AUC",
    "F1 Score",
    "Accuracy",
    "Precision",
    "Recall",
    "Sensitivity",
    "Specificity"
]

print("========== Percentage Variation ==========")
for metric_name in metric_names:
    values = np.array([r["metrics"][metric_name] for r in fold_results], dtype=float)
    min_val = np.min(values)
    max_val = np.max(values)

    if min_val != 0:
        percent_variation = 100.0 * (max_val - min_val) / min_val
    else:
        percent_variation = np.nan

    print(metric_name)
    print("  min:", min_val)
    print("  max:", max_val)
    print("  percent variation:", percent_variation)
print()

import os
import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import confusion_matrix

# ---------------------------------------------------------
# Problem 4: Multiclass classifier for classes A, B, C
# A: 0,1,2
# B: 3,4,5
# C: 6,7,8,9
# ---------------------------------------------------------

os.makedirs("plots", exist_ok=True)
np.random.seed(42)

# load dataset
digits = load_digits()
X = digits.data
y = digits.target

# create multiclass labels:
# 0 = Class A
# 1 = Class B
# 2 = Class C
y_abc = np.zeros_like(y)

mask_b = np.logical_or(np.logical_or(y == 3, y == 4), y == 5)
mask_c = np.logical_or(np.logical_or(np.logical_or(y == 6, y == 7), y == 8), y == 9)

y_abc[mask_b] = 1
y_abc[mask_c] = 2

# optional print to verify class counts
print("Class counts:")
print("Class A:", np.sum(y_abc == 0))
print("Class B:", np.sum(y_abc == 1))
print("Class C:", np.sum(y_abc == 2))
print()

# stratified train/test split
X_train, X_test, y_train_abc, y_test_abc = train_test_split(
    X,
    y_abc,
    test_size=0.2,
    random_state=42,
    stratify=y_abc
)

# train SGD classifier
sgd_clf = SGDClassifier(random_state=42, max_iter=1000, tol=1e-3)
sgd_clf.fit(X_train, y_train_abc)

# predictions
y_pred_abc = sgd_clf.predict(X_test)

# confusion matrix
cm = confusion_matrix(y_test_abc, y_pred_abc)

print("Confusion Matrix:")
print(cm)
print()

# find one correct example
correct_idx = np.where(y_pred_abc == y_test_abc)[0][0]

# find one incorrect example
incorrect_candidates = np.where(y_pred_abc != y_test_abc)[0]
if len(incorrect_candidates) > 0:
    incorrect_idx = incorrect_candidates[0]
else:
    incorrect_idx = None

class_names = {
    0: "Class A",
    1: "Class B",
    2: "Class C"
}

# print correct example info
print("Correct example:")
print("Test index:", correct_idx)
print("True label:", y_test_abc[correct_idx], "-", class_names[y_test_abc[correct_idx]])
print("Predicted label:", y_pred_abc[correct_idx], "-", class_names[y_pred_abc[correct_idx]])
print()

# print incorrect example info
if incorrect_idx is not None:
    print("Incorrect example:")
    print("Test index:", incorrect_idx)
    print("True label:", y_test_abc[incorrect_idx], "-", class_names[y_test_abc[incorrect_idx]])
    print("Predicted label:", y_pred_abc[incorrect_idx], "-", class_names[y_pred_abc[incorrect_idx]])
    print()
else:
    print("No incorrect examples found.")
    print()

# save correct example image
plt.figure()
plt.imshow(X_test[correct_idx].reshape(8, 8), cmap="gray")
plt.title(
    f"Correct Example\nTrue: {class_names[y_test_abc[correct_idx]]}, "
    f"Pred: {class_names[y_pred_abc[correct_idx]]}"
)
plt.axis("off")
plt.savefig("plots/problem4_correct_example.png")
plt.close()

# save incorrect example image
if incorrect_idx is not None:
    plt.figure()
    plt.imshow(X_test[incorrect_idx].reshape(8, 8), cmap="gray")
    plt.title(
        f"Incorrect Example\nTrue: {class_names[y_test_abc[incorrect_idx]]}, "
        f"Pred: {class_names[y_pred_abc[incorrect_idx]]}"
    )
    plt.axis("off")
    plt.savefig("plots/problem4_incorrect_example.png")
    plt.close()