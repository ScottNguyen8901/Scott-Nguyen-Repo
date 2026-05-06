from pathlib import Path

dataset = Path("C:/Users/snguye17/Documents/folder/martianlunar-crater-detection-dataset/craters")

for split in ["train", "valid", "test"]:
    image_dir = dataset / split / "images"
    label_dir = dataset / split / "labels"

    images = list(image_dir.glob("*.jpg"))
    labels = list(label_dir.glob("*.txt"))

    box_count = 0
    empty_labels = 0

    for label in labels:
        lines = label.read_text().strip().splitlines()
        if not lines:
            empty_labels += 1
        box_count += len(lines)

    print(f"{split}:")
    print(f"  images: {len(images)}")
    print(f"  label files: {len(labels)}")
    print(f"  crater boxes: {box_count}")
    print(f"  empty label files: {empty_labels}")
    print()