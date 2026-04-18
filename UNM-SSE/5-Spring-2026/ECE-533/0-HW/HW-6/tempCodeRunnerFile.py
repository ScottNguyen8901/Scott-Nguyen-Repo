
# save incorrect example image
if incorrect_idx is not None:
    plt.figure()
    plt.imshow(X_test[incorrect_idx].reshape(8, 8), cmap="gray")
    plt.title(
        f"Incorrect Example\nTrue: {class_names[y_test_abc[incorrect_idx]]}, "
        f"Pred: {class_names[y_pred_abc[incorrect_idx]]}"