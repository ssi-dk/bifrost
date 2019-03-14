import os
import glob


#DEV images
image_directory = "./assets/img/"
list_of_images = [os.path.basename(x) for x in glob.glob(
    "./assets/img/*.svg")]
static_image_route = "/assets/"