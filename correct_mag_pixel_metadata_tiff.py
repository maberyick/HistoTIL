import tifffile as tf
import numpy as np

def update_metadata(tiff_path):
    # Load the TIFF image
    with tf.TiffFile(tiff_path) as tif:
        # Read the image and metadata
        image = tif.asarray()
        metadata = tif.pages[0].tags
        
        # Check and update 'Magnification' metadata
        if 'Magnification' in metadata:
            magnification = metadata['Magnification'].value
            if magnification == 'Unknown' or not isinstance(magnification, (int, float)):
                metadata['Magnification'] = 40.0
        
        # Check and update 'Pixel width' metadata
        if 'Pixel width' in metadata:
            pixel_width = metadata['Pixel width'].value
            if pixel_width == 'Unknown' or not isinstance(pixel_width, (int, float)):
                metadata['Pixel width'] = 0.25
        
        # Check and update 'Pixel height' metadata
        if 'Pixel height' in metadata:
            pixel_height = metadata['Pixel height'].value
            if pixel_height == 'Unknown' or not isinstance(pixel_height, (int, float)):
                metadata['Pixel height'] = 0.25
        
        # Save the updated metadata to the TIFF image
        new_tiff_path = tiff_path.rstrip('.tif') + '_updated.tif'
        with tf.TiffWriter(new_tiff_path) as tif_out:
            tif_out.save(image, description=tif.pages[0].description, metadata=metadata)
        
        print("Updated TIFF image with new metadata saved successfully.")

# Usage example
tiff_path = '/mnt/datas1/ccf_io/wsi/161.tiff'
update_metadata(tiff_path)