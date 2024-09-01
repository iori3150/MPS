import json
import os
import sys
import time

import cv2


def load_config(config_file):
    if not os.path.exists(config_file):
        print(f"Error: Config file '{config_file}' not found.")
        sys.exit(1)

    with open(config_file, "r") as f:
        config = json.load(f)

    return config


def images_to_video(image_folder, output_video, frame_rate):
    # Get a list of image files and sort them by filename
    images = [
        img
        for img in os.listdir(image_folder)
        if img.endswith((".png", ".jpg", ".jpeg"))
    ]
    images.sort()

    if len(images) == 0:
        print("No images found.")
        return

    # Get the size of the first image to set the video format
    first_image_path = os.path.join(image_folder, images[0])
    frame = cv2.imread(first_image_path)
    height, width, layers = frame.shape

    # Set up the output video
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    video = cv2.VideoWriter(output_video, fourcc, frame_rate, (width, height))

    # Start time for progress tracking
    start_time = time.time()

    # Write images to the video
    total_images = len(images)
    for idx, image in enumerate(images):
        image_path = os.path.join(image_folder, image)
        frame = cv2.imread(image_path)
        video.write(frame)

        # Progress calculation
        current_time = time.time()
        elapsed_time = current_time - start_time
        percentage_complete = (idx + 1) / total_images * 100
        estimated_total_time = elapsed_time / (idx + 1) * total_images
        estimated_time_remaining = estimated_total_time - elapsed_time

        # Log progress
        print(
            f"Processing file: {image} | "
            f"Progress: {percentage_complete:.2f}% | "
            f"Elapsed time: {elapsed_time:.2f}s | "
            f"Estimated time remaining: {estimated_time_remaining:.2f}s"
        )

    # Release the video file
    video.release()
    print(f"Video file created: {output_video}")


config_file = "config.json"
config = load_config(config_file)

image_folder = config["image_folder"]
output_video = config["output_video"]
frame_rate = config["frame_rate"]

base_dir = os.path.dirname(os.path.abspath(__file__))
image_folder = os.path.join(base_dir, config["image_folder"])
output_video = os.path.join(base_dir, config["output_video"])
frame_rate = config["frame_rate"]

# Create video from images
images_to_video(image_folder, output_video, frame_rate)
