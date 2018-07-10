import os
import glob
from components.mongo_interface import get_species_colors


#DEV images
image_directory = '/Users/mbas/Documents/SerumQC-private/reporter/resources/img/'
list_of_images = [os.path.basename(x) for x in glob.glob(
    '{}*.svg'.format(image_directory))]
static_image_route = '/static/'

COLOR_DICT = get_species_colors()
def get_species_color(species):
    color = COLOR_DICT.get(species, "#b3ccc1")
    if color == '':
        color = "#b3ccc1"  # Default
    return color
