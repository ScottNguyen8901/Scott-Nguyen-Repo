import numpy as np
import pandas as pd

from ucimlrepo import fetch_ucirepo
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.impute import SimpleImputer

v_folds = 5
cv = StratifiedKFold(n_splits=v_folds, shuffle=True, random_state=42)

param_grid = [
    {
        "kernel": ["linear"],
        "C": [0.1, 1, 10, 100],
    },
    {
        "kernel": ["rbf"],
        "C": [0.1, 1, 10, 100],
        "gamma": [0.01, 0.1, 1],
    },
    {
        "kernel": ["poly"],
        "C": [0.1, 1, 10],
        "degree": [2, 3, 4],
        "gamma": ["scale", "auto"],
    },
]


def prepare_features(X_train, X_test):

    if hasattr(X_train, "select_dtypes"):
        X_train_enc = pd.get_dummies(X_train)
        X_test_enc = pd.get_dummies(X_test)

        X_train_enc, X_test_enc = X_train_enc.align(
            X_test_enc, join="left", axis=1, fill_value=0
        )
    else:
        X_train_enc, X_test_enc = X_train, X_test

    X_train_enc = X_train_enc.astype(float)
    X_test_enc = X_test_enc.astype(float)

    imputer = SimpleImputer(strategy="most_frequent")
    X_train_imp = imputer.fit_transform(X_train_enc)
    X_test_imp = imputer.transform(X_test_enc)

    return X_train_imp, X_test_imp

dataset_ids = {
    "car_evaluation": 19,
    "dermatology": 33,
    "glass": 42,
    "wine": 109,
}

datasets = {}
splits = {}

for name, uci_id in dataset_ids.items():
    print(f"\n=== {name} (UCI id={uci_id}) ===")
    ds = fetch_ucirepo(id=uci_id)

    X = ds.data.features
    y = ds.data.targets.squeeze()

    if getattr(y, "ndim", 1) > 1:
        y = y.iloc[:, 0] if hasattr(y, "iloc") else y[:, 0]

    n_samples, n_features = X.shape
    classes = np.unique(y)
    n_classes = len(classes)

    print(f"samples: {n_samples}, features: {n_features}, classes: {n_classes}")
    print("class labels:", classes)
    assert n_classes > 2, f"{name} is not multiclass"

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.2,
        stratify=y,
        random_state=42,
    )

    datasets[name] = ds
    splits[name] = {
        "X_train": X_train,
        "X_test": X_test,
        "y_train": y_train,
        "y_test": y_test,
    }

    print(f"train shape: {X_train.shape}, test shape: {X_test.shape}")

easy_datasets = ["car_evaluation", "wine"]
difficult_datasets = ["dermatology", "glass"]

print("\nEasy datasets:", easy_datasets)
print("Difficult datasets:", difficult_datasets)

results = {}

for name, split in splits.items():
    print("\n" + "=" * 70)
    print(f"SVM experiment for dataset: {name}")
    print("=" * 70)

    X_train = split["X_train"]
    X_test = split["X_test"]
    y_train = split["y_train"]
    y_test = split["y_test"]

    X_train_pre, X_test_pre = prepare_features(X_train, X_test)

    grid = GridSearchCV(
        estimator=SVC(),
        param_grid=param_grid,
        cv=cv,
        scoring="accuracy",
        n_jobs=-1,
        refit=True, 
        return_train_score=False,
    )

    grid.fit(X_train_pre, y_train)

    best_model = grid.best_estimator_
    best_params = grid.best_params_
    best_cv_acc = grid.best_score_

    print("Best parameters (v-fold CV):", best_params)
    print(f"Best mean CV accuracy: {best_cv_acc:.4f}")

    cv_results = grid.cv_results_
    kernel_best_cv = {}
    for mean_score, params in zip(cv_results["mean_test_score"], cv_results["params"]):
        k = params["kernel"]
        if k not in kernel_best_cv or mean_score > kernel_best_cv[k]:
            kernel_best_cv[k] = mean_score

    print("\nBest CV accuracy per kernel:")
    for k, acc in kernel_best_cv.items():
        print(f"  {k:>6}: {acc:.4f}")

    y_pred_test = best_model.predict(X_test_pre)
    test_acc = accuracy_score(y_test, y_pred_test)
    cm = confusion_matrix(y_test, y_pred_test)

    print(f"\nTest accuracy (on 20% hold-out set): {test_acc:.4f}")
    print("Confusion matrix (rows = true, cols = predicted):")
    print(cm)

    results[name] = {
        "best_params": best_params,
        "best_cv_acc": best_cv_acc,
        "kernel_best_cv": kernel_best_cv,
        "test_acc": test_acc,
        "confusion_matrix": cm,
    }

print("\n\n=== Summary of results across all datasets ===")
for name, res in results.items():
    print("\nDataset:", name)
    print("  Best params:", res["best_params"])
    print(f"  Best CV accuracy: {res['best_cv_acc']:.4f}")
    print(f"  Test accuracy:    {res['test_acc']:.4f}")
    print("  Best CV accuracy per kernel:")
    for k, acc in res["kernel_best_cv"].items():
        print(f"    {k:>6}: {acc:.4f}")