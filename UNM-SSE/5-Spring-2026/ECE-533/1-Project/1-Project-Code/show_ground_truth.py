import cv2
from pathlib import Path

dataset = Path("C:/Users/snguye17/Documents/folder/martianlunar-crater-detection-dataset/craters")
image_dir = dataset / "test/images"
label_dir = dataset / "test/labels"
output_dir = Path("ground_truth_examples")
output_dir.mkdir(exist_ok=True)

for image_path in list(image_dir.glob("*.jpg"))[:5]:
    img = cv2.imread(str(image_path))
    h, w = img.shape[:2]

    label_path = label_dir / f"{image_path.stem}.txt"

    if label_path.exists():
        for line in label_path.read_text().strip().splitlines():
            cls, xc, yc, bw, bh = map(float, line.split())

            x1 = int((xc - bw / 2) * w)
            y1 = int((yc - bh / 2) * h)
            x2 = int((xc + bw / 2) * w)
            y2 = int((yc + bh / 2) * h)

            cv2.rectangle(img, (x1, y1), (x2, y2), (0, 255, 0), 2)
            cv2.putText(img, "crater", (x1, max(y1 - 5, 15)),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 0), 1)

    cv2.imwrite(str(output_dir / image_path.name), img)

print(f"saved ground truth examples to {output_dir.resolve()}")