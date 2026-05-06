from ultralytics import YOLO

model = YOLO("runs/detect/crater_detector/weights/best.pt")

metrics = model.val(
    data="C:/Users/snguye17/Documents/folder/martianlunar-crater-detection-dataset/craters/data.yaml",
    split="test",
    imgsz=640
)