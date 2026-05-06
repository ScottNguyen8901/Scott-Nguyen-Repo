from ultralytics import YOLO

model = YOLO("runs/detect/crater_detector/weights/best.pt")

model.predict(
    source="C:/Users/snguye17/Documents/folder/martianlunar-crater-detection-dataset/craters/test/images",
    imgsz=640,
    conf=0.25,
    save=True,
    name="crater_test_predictions"
)