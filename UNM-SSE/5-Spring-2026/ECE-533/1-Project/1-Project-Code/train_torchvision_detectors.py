from pathlib import Path

import torch
from torch.utils.data import Dataset, DataLoader
from PIL import Image
from tqdm import tqdm

from torchvision.transforms import functional as F
from torchvision.models.detection import (
    fasterrcnn_resnet50_fpn,
    retinanet_resnet50_fpn,
)
from torchvision.models.detection.faster_rcnn import FastRCNNPredictor
from torchvision.models.detection.retinanet import RetinaNetClassificationHead


DATASET = Path(
    "C:/Users/snguye17/Documents/folder/martianlunar-crater-detection-dataset/craters"
)
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class CraterDataset(Dataset):
    def __init__(self, split):
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
            lines = label_path.read_text().strip().splitlines()

            for line in lines:
                cls, xc, yc, bw, bh = map(float, line.split())

                x1 = (xc - bw / 2) * w
                y1 = (yc - bh / 2) * h
                x2 = (xc + bw / 2) * w
                y2 = (yc + bh / 2) * h

                boxes.append([x1, y1, x2, y2])
                labels.append(1)  # 1 = crater, 0 = background

        if len(boxes) == 0:
            boxes = torch.zeros((0, 4), dtype=torch.float32)
            labels = torch.zeros((0,), dtype=torch.int64)
        else:
            boxes = torch.as_tensor(boxes, dtype=torch.float32)
            labels = torch.as_tensor(labels, dtype=torch.int64)

        target = {
            "boxes": boxes,
            "labels": labels,
            "image_id": torch.tensor([idx]),
        }

        image = F.to_tensor(image)
        return image, target


def collate_fn(batch):
    return tuple(zip(*batch))


def get_faster_rcnn():
    model = fasterrcnn_resnet50_fpn(weights="DEFAULT")

    num_classes = 2  # background + crater
    in_features = model.roi_heads.box_predictor.cls_score.in_features
    model.roi_heads.box_predictor = FastRCNNPredictor(in_features, num_classes)

    return model


def get_retinanet():
    model = retinanet_resnet50_fpn(weights="DEFAULT")

    num_classes = 2  # background + crater
    num_anchors = model.head.classification_head.num_anchors

    # torchvision version compatibility fix
    in_channels = model.head.classification_head.conv[0].out_channels

    model.head.classification_head = RetinaNetClassificationHead(
        in_channels=in_channels,
        num_anchors=num_anchors,
        num_classes=num_classes,
    )

    return model

def train_model(model, model_name, epochs=3, batch_size=2, lr=1e-4):
    model.to(DEVICE)

    train_ds = CraterDataset("train")
    train_loader = DataLoader(
        train_ds,
        batch_size=batch_size,
        shuffle=True,
        collate_fn=collate_fn,
    )

    optimizer = torch.optim.AdamW(model.parameters(), lr=lr)

    model.train()

    for epoch in range(epochs):
        total_loss = 0.0

        progress_bar = tqdm(
            train_loader,
            desc=f"{model_name} epoch {epoch + 1}/{epochs}",
            unit="batch",
        )

        for images, targets in progress_bar:
            images = [img.to(DEVICE) for img in images]
            targets = [{k: v.to(DEVICE) for k, v in t.items()} for t in targets]

            loss_dict = model(images, targets)
            loss = sum(loss for loss in loss_dict.values())

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            progress_bar.set_postfix(loss=f"{loss.item():.4f}")

        avg_loss = total_loss / len(train_loader)
        print(
            f"{model_name} epoch {epoch + 1}/{epochs}, "
            f"total loss: {total_loss:.4f}, avg loss: {avg_loss:.4f}"
        )

    output_path = f"{model_name}.pth"
    torch.save(model.state_dict(), output_path)
    print(f"saved {output_path}")


if __name__ == "__main__":
    print(f"Using device: {DEVICE}")

    faster_rcnn = get_faster_rcnn()
    train_model(faster_rcnn, "faster_rcnn_crater", epochs=3)

    retinanet = get_retinanet()
    train_model(retinanet, "retinanet_crater", epochs=3)