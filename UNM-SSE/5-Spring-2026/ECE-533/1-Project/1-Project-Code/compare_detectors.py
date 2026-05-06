from pathlib import Path

import pandas as pd
import torch
from PIL import Image
from torch.utils.data import Dataset, DataLoader
from torchmetrics.detection.mean_ap import MeanAveragePrecision
from torchvision.transforms import functional as F
from torchvision.models.detection import fasterrcnn_resnet50_fpn, retinanet_resnet50_fpn
from torchvision.models.detection.faster_rcnn import FastRCNNPredictor
from torchvision.models.detection.retinanet import RetinaNetClassificationHead
from ultralytics import YOLO

import matplotlib.pyplot as plt


DATASET = Path("C:/Users/snguye17/Documents/folder/martianlunar-crater-detection-dataset/craters")
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class CraterDataset(Dataset):
    def __init__(self, split="test"):
        self.image_dir = DATASET / split / "images"
        self.label_dir = DATASET / split / "labels"
        self.images = sorted(self.image_dir.glob("*.jpg"))

    def __len__(self):
        return len(self.images)

    def __getitem__(self, idx):
        image_path = self.images[idx]
        image = Image.open(image_path).convert("RGB")
        w, h = image.size

        label_path = self.label_dir / f"{image_path.stem}.txt"

        boxes = []
        labels = []

        if label_path.exists():
            for line in label_path.read_text().strip().splitlines():
                cls, xc, yc, bw, bh = map(float, line.split())

                x1 = (xc - bw / 2) * w
                y1 = (yc - bh / 2) * h
                x2 = (xc + bw / 2) * w
                y2 = (yc + bh / 2) * h

                boxes.append([x1, y1, x2, y2])
                labels.append(1)

        if len(boxes) == 0:
            boxes = torch.zeros((0, 4), dtype=torch.float32)
            labels = torch.zeros((0,), dtype=torch.int64)
        else:
            boxes = torch.as_tensor(boxes, dtype=torch.float32)
            labels = torch.as_tensor(labels, dtype=torch.int64)

        return F.to_tensor(image), {"boxes": boxes, "labels": labels}


def collate_fn(batch):
    return tuple(zip(*batch))


def get_faster_rcnn():
    model = fasterrcnn_resnet50_fpn(weights=None, weights_backbone=None)
    num_classes = 2

    in_features = model.roi_heads.box_predictor.cls_score.in_features
    model.roi_heads.box_predictor = FastRCNNPredictor(in_features, num_classes)

    return model


def get_retinanet():
    model = retinanet_resnet50_fpn(weights=None, weights_backbone=None)
    num_classes = 2

    num_anchors = model.head.classification_head.num_anchors
    in_channels = model.head.classification_head.conv[0].out_channels

    model.head.classification_head = RetinaNetClassificationHead(
        in_channels=in_channels,
        num_anchors=num_anchors,
        num_classes=num_classes,
    )

    return model


def evaluate_torchvision_model(model, weights_path, model_name):
    model.load_state_dict(torch.load(weights_path, map_location=DEVICE))
    model.to(DEVICE)
    model.eval()

    test_ds = CraterDataset("test")
    test_loader = DataLoader(
        test_ds,
        batch_size=1,
        shuffle=False,
        collate_fn=collate_fn,
    )

    metric = MeanAveragePrecision(box_format="xyxy", iou_type="bbox")

    with torch.no_grad():
        for images, targets in test_loader:
            images = [img.to(DEVICE) for img in images]
            outputs = model(images)

            preds = []
            gts = []

            for output, target in zip(outputs, targets):
                preds.append(
                    {
                        "boxes": output["boxes"].cpu(),
                        "scores": output["scores"].cpu(),
                        "labels": output["labels"].cpu(),
                    }
                )

                gts.append(
                    {
                        "boxes": target["boxes"].cpu(),
                        "labels": target["labels"].cpu(),
                    }
                )

            metric.update(preds, gts)

    results = metric.compute()

    return {
        "model": model_name,
        "mAP50": float(results["map_50"]),
        "mAP50_95": float(results["map"]),
        "mAP75": float(results["map_75"]),
    }


def evaluate_yolo():
    model = YOLO("runs/detect/crater_detector/weights/best.pt")

    metrics = model.val(
        data=str(DATASET / "data.yaml"),
        split="test",
        imgsz=640,
        verbose=False,
    )

    return {
        "model": "YOLOv8",
        "mAP50": float(metrics.box.map50),
        "mAP50_95": float(metrics.box.map),
        "mAP75": float(metrics.box.map75),
    }


def make_plot(df):
    ax = df.plot(
        x="model",
        y=["mAP50", "mAP50_95", "mAP75"],
        kind="bar",
        figsize=(9, 5),
    )

    ax.set_title("Crater Detection Model Comparison")
    ax.set_xlabel("Model")
    ax.set_ylabel("Score")
    ax.set_ylim(0, 1)
    ax.legend(["mAP@0.5", "mAP@0.5:0.95", "mAP@0.75"])
    plt.xticks(rotation=0)
    plt.tight_layout()

    plt.savefig("model_comparison_plot.png", dpi=300)
    print("Saved plot to model_comparison_plot.png")


if __name__ == "__main__":
    print(f"Using device: {DEVICE}")

    results = []

    print("Evaluating YOLOv8...")
    results.append(evaluate_yolo())

    print("Evaluating Faster R-CNN...")
    faster_rcnn = get_faster_rcnn()
    results.append(
        evaluate_torchvision_model(
            faster_rcnn,
            "faster_rcnn_crater.pth",
            "Faster R-CNN",
        )
    )

    print("Evaluating RetinaNet...")
    retinanet = get_retinanet()
    results.append(
        evaluate_torchvision_model(
            retinanet,
            "retinanet_crater.pth",
            "RetinaNet",
        )
    )

    df = pd.DataFrame(results)
    print()
    print(df)

    df.to_csv("model_comparison_results.csv", index=False)
    print("Saved results to model_comparison_results.csv")

    make_plot(df)