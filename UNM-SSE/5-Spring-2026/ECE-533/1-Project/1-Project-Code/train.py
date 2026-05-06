from ultralytics import YOLO

model = YOLO("yolov8n.pt")

model.train(
    data="C:/Users/snguye17/Documents/folder/martianlunar-crater-detection-dataset/craters/data.yaml",
    epochs=50,
    imgsz=640,
    batch=8,
    name="crater_detector"
)