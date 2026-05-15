from pathlib import Path

import torch
from PIL import Image, ImageDraw, ImageFont
from torchvision.transforms import functional as F
from torchvision.models.detection import fasterrcnn_resnet50_fpn, retinanet_resnet50_fpn
from torchvision.models.detection.faster_rcnn import FastRCNNPredictor
from torchvision.models.detection.retinanet import RetinaNetClassificationHead


DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

IMAGE_PATH = Path(
    "martianlunar-crater-detection-dataset/craters/test/images/"
    "015_png.rf.7d5b2091b6339c9480a171a59c52c3b9.jpg"
)

OUT_DIR = Path("torchvision_predictions")
OUT_DIR.mkdir(exist_ok=True)


def get_faster_rcnn():
    model = fasterrcnn_resnet50_fpn(weights=None, weights_backbone=None)

    num_classes = 2  # background + crater
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


def draw_predictions(model, weights_path, output_name, score_thresh=0.3):
    model.load_state_dict(torch.load(weights_path, map_location=DEVICE))
    model.to(DEVICE)
    model.eval()

    image = Image.open(IMAGE_PATH).convert("RGB")
    image_tensor = F.to_tensor(image).to(DEVICE)

    with torch.no_grad():
        output = model([image_tensor])[0]

    draw = ImageDraw.Draw(image)

    boxes = output["boxes"].cpu()
    scores = output["scores"].cpu()

    count = 0

    for box, score in zip(boxes, scores):
        if score < score_thresh:
            continue

        x1, y1, x2, y2 = box.tolist()

        draw.rectangle(
            [x1, y1, x2, y2],
            outline="red",
            width=3,
        )

        draw.text(
            (x1, max(0, y1 - 12)),
            f"crater {score:.2f}",
            fill="red",
        )

        count += 1

    out_path = OUT_DIR / output_name
    image.save(out_path)

    print(f"Saved {out_path}")
    print(f"{output_name}: {count} detections above threshold {score_thresh}")


if __name__ == "__main__":
    print(f"Using device: {DEVICE}")

    faster_rcnn = get_faster_rcnn()
    draw_predictions(
        faster_rcnn,
        "faster_rcnn_crater.pth",
        "faster_rcnn_prediction_015.jpg",
        score_thresh=0.3,
    )

    retinanet = get_retinanet()
    draw_predictions(
        retinanet,
        "retinanet_crater.pth",
        "retinanet_prediction_015.jpg",
        score_thresh=0.3,
    )